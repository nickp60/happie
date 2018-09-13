#!/usr/bin/env python3
#-*- coding: utf-8 -*-
import re
import argparse
import sys
import os
import unittest
import itertools
import multiprocessing
import subprocess

from Bio.Seq import Seq
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

from . import __version__


def get_args():  # pragma: no cover
    parser = argparse.ArgumentParser(
        description="extract the mobile elements and run pangenome analysis on",
        add_help=False)
    parser.add_argument("-c", "--contigs", action="store",
                        help="FASTA formatted genome or set of contigs",
                        required="-p" not in sys.argv)
    parser.add_argument("-p", "--prokka_dir", action="store",
                        help="prokka file",
                        required="-c" not in sys.argv)
    parser.add_argument("-o", "--output", action="store",
                        help="destination dir", required=True)

    optional = parser.add_argument_group('optional arguments')
    optional.add_argument("--cores", dest='cores', default=4,
                          help="cores to use" +
                          "without extension.")
    optional.add_argument("-n", "--name", dest='name',
                          help="name of experiment; defaults to file name " +
                          "without extension.", default="test")
    optional.add_argument("-h", "--help",
                          action="help", default=argparse.SUPPRESS,
                          help="Displays this help message")
    optional.add_argument('--version',
                          action='version',
                          version='%(prog)s ' + __version__)
    args = parser.parse_args()
    return args


def main(args=None):
    if args is None:
        args = get_args()
    output_root = os.path.abspath(os.path.expanduser(args.output))
    os.makedirs(output_root, exist_ok=False)
    # make output dir names
    prophet_dir = os.path.join(output_root, "ProphET")
    prophet_results = os.path.join(prophet_dir, "phages_coords")
    mobsuite_dir = os.path.join(output_root, "mobsuite")
    island_dir = os.path.join(output_root, "dimob")
    island_results = os.path.join(output_root, "dimob", "results.txt")
    mlplasmids_dir = os.path.join(output_root, "mlplasmids")
    mlplasmids_results = os.path.join(output_root, "mlplasmids", "results.txt")
    for path in [prophet_dir, mobsuite_dir, mlplasmids_dir]:
        os.makedirs(path, exist_ok=False)

    if args.cores is None:
        args.cores = multiprocessing.cpu_count()

    if args.prokka_dir is None:
        args.prokka_dir = os.path.join(output_root, "prokka")
        print("running prokka")
        prokka_cmd = "{exe} {file} -o {outdir} --prefix {name} --fast --cpus {cpus}".format(
            exe="prokka", file=args.contigs, outdir=args.prokka_dir, name=args.name, cpus=args.cores)
        print(prokka_cmd)
        subprocess.run([prokka_cmd],
                       shell=sys.platform != "win32",
                       stdout=subprocess.PIPE,
                       stderr=subprocess.PIPE,
                       check=True)
    else:
        args.prokka_dir = os.path.abspath(os.path.expanduser(args.prokka_dir))
        args.name = os.path.splitext(os.path.basename(args.prokka_dir))[0]

    # set some names of shtuff
    prokka_fna = os.path.join(args.prokka_dir, args.name + ".fna")
    prokka_gbk = os.path.join(args.prokka_dir, args.name + ".gbk")
    prokka_gff = os.path.join(args.prokka_dir, args.name + ".gff")
    prokka_new_gff=os.path.join(args.prokka_dir, args.name + "_new.gbk")

    print( "Reformatting gff")
    print(os.getcwd())

    ###########################################################################
    # try not to vomit
    os.chdir("./submodules/ProphET/UTILS.dir/GFFLib/")
    print(os.getcwd())
    new_gff_cmd = "{exe} --input {file} --output {out} --add_missing_features".format(
        exe="./gff_rewrite.pl", file=prokka_gff, out=prokka_new_gff)
    subprocess.run([new_gff_cmd],
                   shell=sys.platform != "win32",
                   stdout=subprocess.PIPE,
                   stderr=subprocess.PIPE,
                   check=True)
                   #cwd=os.getcwd())
    os.chdir("../../../../")
    print(os.getcwd())
    ###########################################################################

    prophet_cmd = "{exe} --fasta_in  {file} --gff_in {gff} --outdir {out}".format(
        exe="./submodules/ProphET/ProphET_standalone.pl", file=prokka_fna, gff=prokka_new_gff, out=prophet_dir)
    print(prophet_cmd)
    subprocess.run([prophet_cmd],
                   shell=sys.platform != "win32",
                   stdout=subprocess.PIPE,
                   stderr=subprocess.PIPE,
                   check=True)

    # mobsuite_cmd = "{exe} --infile {file} --outdir {out} --run_typer --keep_tmp".format(
    #     exe="mob_recon", file=prokka_fna, out=mobsuite_dir)
    # print(mobsuite_cmd)
    # subprocess.run([mobsuite_cmd],
    #                shell=sys.platform != "win32",
    #                stdout=subprocess.PIPE,
    #                stderr=subprocess.PIPE,
    #                check=True)



    mlplasmids_cmd = "{exe} {file} {out} .8 'Escherichia coli'".format(
        exe="Rscript scripts/run_mlplasmids.R", file=prokka_fna, out=mlplasmids_results)
    print(mlplasmids_cmd)
    subprocess.run([mlplasmids_cmd],
                   shell=sys.platform != "win32",
                   stdout=subprocess.PIPE,
                   stderr=subprocess.PIPE,
                   check=True)



    ###########################################################################
    # program type sequence start end
    all_results = []
    prophet_results_text = []
    mlplasmids_results_text = []
    with open(prophet_results) as inf:
        for line in inf:
            results = ["prophet", "prophage"]
            results.extend(line.strip().split("\t"))
            prophet_results_text.append(results)
    with open(mlplasmids_results) as inf:
        for line in inf:
            # here we add 1 as a start index, and we will use the contig length as the end
            results = ["mlplasmids", "plasmid", "0"]
            results.extend(line.strip().split("\t"))
            mlplasmids_results_text.append(results)
    print(prophet_results_text)
    print(mlplasmids_results_text)
    for line in prophet_results_text:
        subresults = []
        for i in [0, 1, 2, 4, 5]:
            subresults.append(line[i])
        all_results.append(subresults)
    for line in mlplasmids_results_text:
        if line[5] == '"Plasmid"':
            subresults = []
            for i in [0, 1, 6, 2, 7]:
                subresults.append(line[i])
            all_results.append(subresults)

    print("sjjvajsjvajsjvavjjas")
    print(all_results)


if __name__ == "__main__":
    args = get_args()
    main(args)
