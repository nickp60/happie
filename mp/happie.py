#!/usr/bin/env python3
#-*- coding: utf-8 -*-
import re
import argparse
import sys
import shutil
import os
import unittest
import itertools
import multiprocessing
import subprocess
import glob

from Bio.Seq import Seq
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

from . import __version__


def get_args():  # pragma: no cover
    parser = argparse.ArgumentParser(
        description="extract the mobile elements for pangenome analysis",
        add_help=False)
    parser.add_argument("--contigs", action="store",
                        help="FASTA formatted genome or set of contigs",
                        required="-p" not in sys.argv)
    parser.add_argument("-p", "--prokka_dir", action="store",
                        help="prokka file",
                        required="--contigs" not in sys.argv)
    parser.add_argument("-o", "--output", action="store",
                        help="destination dir", required=True)

    optional = parser.add_argument_group('optional arguments')
    optional.add_argument("--cores", dest='cores', default=4,
                          help="cores to use" +
                          "without extension.")
    optional.add_argument("-n", "--name", dest='name',
                          help="name of experiment; defaults to file name " +
                          "without extension.", default="test")
    optional.add_argument("--skip_dimob", dest='skip_dimob',
                          action="store_true",
                          help="skip island finding.")
    optional.add_argument("--skip_rename", dest='skip_rename',
                          action="store_true",
                          help="skip initial contig renaming.")
    optional.add_argument("--elements", dest='elements',
                          action="store", nargs="*",
                          default=["plasmids", "islands", "prophages"],
                          choices=["plasmids", "islands", "prophages"],
                          help="skip initial contig renaming.")
    optional.add_argument("-s", "--restart_stage", dest='restart_stage',
                          choices=[1, 2, 3, 4, 5],
                          help="stage, if restarting:" +
                          "1 - from the begining, run prokka |" +
                          "2 - run prophage finders |" +
                          "3 - run plasmid finders | " +
                          "4 - run genomic island finders | " +
                          "5 - analyze the results",
                          type=int,
                          default=1)
    optional.add_argument("-h", "--help",
                          action="help", default=argparse.SUPPRESS,
                          help="Displays this help message")
    optional.add_argument('--version',
                          action='version',
                          version='%(prog)s ' + __version__)
    args = parser.parse_args()
    return args


def test_exes(exes):
    for exe in exes:
        if shutil.which(exe) is None:
            raise ValueError("%s executable not found")

def parse_prophet_results(prophet_results):
    prophet_results_text = []
    templated = []
    with open(prophet_results) as inf:
        for line in inf:
            results = ["prophet", "prophage"]
            results.extend(line.strip().split("\t"))
            prophet_results_text.append(results)
    for line in prophet_results_text:
        subresults = []
        for i in [0, 1, 2, 4, 5]:
            subresults.append(line[i])
        templated.append(subresults)
    return templated


def parse_mlplasmids_results(mlplasmids_results):
    mlplasmids_results_text = []
    templated = []
    with open(mlplasmids_results) as inf:
        for line in inf:
            # here we add 1 as a start index, and we will use the
            #  contig length as the end.  Its not great, but its what we have
            results = ["mlplasmids", "plasmid", "0"]
            results.extend(line.strip().split("\t"))
            mlplasmids_results_text.append(results)
    for line in mlplasmids_results_text:
        if line[5] == '"Plasmid"':
            subresults = []
            for i in [0, 1, 6, 2, 7]:
                subresults.append(line[i])
            templated.append(subresults)
    return templated


def parse_cafe_results(cafe_results):
    mlplasmids_results_text = []
    templated = []
    with open(mlplasmids_results) as inf:
        for line in inf:
            # here we add 1 as a start index, and we will use the
            #  contig length as the end.  Its not great, but its what we have
            results = ["mlplasmids", "plasmid", "0"]
            results.extend(line.strip().split("\t"))
            mlplasmids_results_text.append(results)
    for line in mlplasmids_results_text:
        if line[5] == '"Plasmid"':
            subresults = []
            for i in [0, 1, 6, 2, 7]:
                subresults.append(line[i])
            templated.append(subresults)
    return templated


def write_sequence_regions_of_interest(contigs, output_path,  all_results):
    """write out all regions of interest to a single fasta file
    """
    total_length = 0
    with open(contigs, "r") as inf, open(output_path, "w") as outf:
        for rec in SeqIO.parse(inf, "fasta"):
            for program, type, recid, start, stop in all_results:
                if rec.id == recid:
                    start, stop = int(start), int(stop)
                    header = \
                        "lcl|{program}|{type}|{recid}:{start}-{stop}".format(
                            **locals())
                    total_length = total_length + (stop - start)
                    SeqIO.write(
                        SeqRecord(
                            # account for 1 index in ext tools
                            rec.seq[start - 1: stop - 1],
                            id=header,
                            description="",
                            name=""
                        ),
                        outf, "fasta"
                    )
    print("wrote out %i bases" % total_length)


def main(args=None):
    if args is None:
        args = get_args()
    output_root = os.path.abspath(os.path.expanduser(args.output))
    if args.restart_stage == 1:
        os.makedirs(output_root, exist_ok=False)
    # make output dir names
    if args.prokka_dir is None:
        args.prokka_dir = os.path.join(output_root, "prokka")
    prophet_dir = os.path.join(output_root, "ProphET", "")
    prophet_results = os.path.join(prophet_dir, "phages_coords")
    mobsuite_dir = os.path.join(output_root, "mobsuite", "")
    island_dir = os.path.join(output_root, "dimob", "")
    island_results = os.path.join(output_root, "dimob", "results.txt")
    mlplasmids_dir = os.path.join(output_root, "mlplasmids", "")
    mlplasmids_results = os.path.join(output_root, "mlplasmids", "results.txt")
    cafe_dir = os.path.join(output_root, "cafe", "")
    cafe_results = os.path.join(output_root, "cafe", "results.txt")
    test_exes(exes=["prokka"])
    if args.restart_stage < 2:
        for path in [prophet_dir, mobsuite_dir, mlplasmids_dir]:
            os.makedirs(path, exist_ok=False)

    if args.cores is None:
        args.cores = multiprocessing.cpu_count()
    if args.restart_stage < 2:
        if not args.skip_rename:
            dest_fasta = os.path.join(output_root, "new_fasta.fasta")
            with open(args.contigs, "r") as inf, open(dest_fasta, "w") as outf:
                for i, rec in enumerate(SeqIO.parse(inf, "fasta")):
                    header = "lcl_" + str(i+1)
                    SeqIO.write(
                        SeqRecord(
                            rec.seq,
                            id=header,
                            description="",
                            name=""
                        ),
                        outf, "fasta"
                    )
            args.contigs = dest_fasta
        print("running prokka")
        prokka_cmd = "{exe} {file} -o {outdir} --prefix {name} --fast --cpus {cpus} 2> {log}".format(
            exe=shutil.which("prokka"), file=args.contigs, outdir=args.prokka_dir, name=args.name, cpus=args.cores,
            log=os.path.join(output_root, "log.txt"))
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
    try:
        prokka_fna = glob.glob(os.path.join(args.prokka_dir, "*.fna"))[0]
        prokka_gbk = glob.glob(os.path.join(args.prokka_dir, "*.gbk"))[0]
        prokka_gff = glob.glob(os.path.join(args.prokka_dir, "*.gff"))[0]
    except IndexError:
        raise IndexError("File not found - something went wrong in step 1 with prokak")
    # for path in [prokka_fna, prokka_gbk, prokka_gff]:
    #     if not os.path.exists(path):
    #         raise ValueError("File not found - something went wrong in step 1 with prokka")


    prokka_new_gff = os.path.splitext(prokka_gbk)[0] + "_new.gff"
    if args.restart_stage < 3 and any([x=="prophages" for x in args.elements]):
        print( "Reformatting gff")
        print(os.getcwd())

        #######################################################################
        # try not to vomit
        os.chdir("./submodules/ProphET/UTILS.dir/GFFLib/")
        print(os.getcwd())
        new_gff_cmd = \
            "{exe} --input {file} --output {o} --add_missing_features 2> {log}".format(
                exe="./gff_rewrite.pl",
                file=prokka_gff,
                o=prokka_new_gff,
                log=os.path.join(output_root, "gff_rewrite_log.txt"))
        subprocess.run([new_gff_cmd],
                       shell=sys.platform != "win32",
                       stdout=subprocess.PIPE,
                       stderr=subprocess.PIPE,
                       check=True)
        #cwd=os.getcwd())
        os.chdir("../../../../")
        print(os.getcwd())
        #######################################################################
        shutil.rmtree(prophet_dir)
        prophet_cmd = "{exe} --fasta_in  {file} --gff_in {gff} --outdir {out} 2> {log}".format(
            exe="perl ./submodules/ProphET/ProphET_standalone.pl",
            file=prokka_fna, gff=prokka_new_gff, out=prophet_dir,
            log=os.path.join(output_root, "ProphET_log.txt"))
        print(prophet_cmd)
        subprocess.run([prophet_cmd],
                       shell=sys.platform != "win32",
                       stdout=subprocess.PIPE,
                       stderr=subprocess.PIPE,
                       check=True)
    if not os.path.exists(prokka_new_gff):
        raise ValueError("Issue with recreating gff %s" % prokka_new_gff)

    if args.restart_stage < 4 and any([x=="plasmids" for x in args.elements]):
        # mobsuite_cmd = \
        #     "{exe} --infile {file} --outdir {out} --run_typer --keep_tmp".format(
        #         exe="mob_recon", file=prokka_fna, out=mobsuite_dir)
        # print(mobsuite_cmd)
        # shutil.rmtree(mobsuite_dir)
        # subprocess.run([mobsuite_cmd],
        #                shell=sys.platform != "win32",
        #                stdout=subprocess.PIPE,
        #                stderr=subprocess.PIPE,
        #                check=True)
        mlplasmids_cmd = "{exe} {file} {out} .8 'Escherichia coli'".format(
            exe="Rscript scripts/run_mlplasmids.R",
            file=prokka_fna, out=mlplasmids_results)
        print(mlplasmids_cmd)
        try: # to get rid of old results if we are rerunning
            os.remove(mlplasmids_results)
        except FileNotFoundError:
            pass
        subprocess.run([mlplasmids_cmd],
                       shell=sys.platform != "win32",
                       stdout=subprocess.PIPE,
                       stderr=subprocess.PIPE,
                       check=True)
    if args.restart_stage < 5 and any([x=="islands" for x in args.elements]):
        #######################################################################
        # try not to vomit again
        # os.chdir("./submodules/CAFE")
        # print(os.getcwd())
        cafe_cmd = \
            "{exe} -gbk {file} --out {o} --verbose 2> {log}".format(
                exe="perl cafe",
                file=prokka_gbk,
                o=cafe_results,
                log=os.path.join(output_root, "cafe_log.txt"))
        subprocess.run([cafe_cmd],
                       shell=sys.platform != "win32",
                       stdout=subprocess.PIPE,
                       stderr=subprocess.PIPE,
                       check=True)
        #cwd=os.getcwd())
        os.chdir("../../")
        print(os.getcwd())
        #######################################################################

    for path in [mlplasmids_results]:
        if not os.path.exists(prokka_new_gff):
            raise ValueError("Issue running mlplasmids gff")
    if args.restart_stage < 5 and not args.skip_dimob:
        island_cmd = "{exe} {file} {out}".format(
            exe="perl ./submodules/islandpath/Dimob.pl",
            file=prokka_gbk, out=island_results)
        print(island_cmd)
        subprocess.run([island_cmd],
                       shell=sys.platform != "win32",
                       stdout=subprocess.PIPE,
                       stderr=subprocess.PIPE,
                       check=True)


    ###########################################################################
    # program type sequence start end
    all_results = []
    prophet_parsed_results = parse_prophet_results(prophet_results)
    mlplasmids_parsed_results = parse_mlplasmids_results(mlplasmids_results)

    for i in [prophet_parsed_results, mlplasmids_parsed_results]:
        all_results.extend(i)

    print(all_results)
    output_path = os.path.join(output_root, "mobile_genome.fasta")
    write_sequence_regions_of_interest(
        contigs=args.contigs,
        output_path=output_path,
        all_results=all_results)


if __name__ == "__main__":
    args = get_args()
    main(args)
