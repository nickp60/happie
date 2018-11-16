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
import pkg_resources

from argparse import Namespace
from Bio.Seq import Seq
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

from . import __version__
from . import shared_methods as sm

def get_args():  # pragma: no cover
    parser = argparse.ArgumentParser(
        description="extract the mobile elements for pangenome analysis",
        add_help=False)
    parser.add_argument("--contigs", action="store",
                        help="FASTA formatted genome or set of contigs",
                        required="-p" not in sys.argv)
    parser.add_argument("-o", "--output", action="store",
                        help="destination dir", required=True)

    parser.add_argument("--virtualization", action="store",
                        help="Whether this will be run with Docker or Singularity",
                        choices=["docker", "singularity"],
                        default="docker")
    optional = parser.add_argument_group('optional arguments')
    optional.add_argument("--cores", dest='cores',
                          default=4,
                          help="cores to use" +
                          "without extension.")
    optional.add_argument("--memory", dest='memory',
                          default=4,
                          help="memory limit to use" +
                          "without extension.")
    optional.add_argument("-n", "--name", dest='name',
                          help="name of experiment; defaults to file name " +
                          "without extension.", default="test")
    optional.add_argument("--skip_rename", dest='skip_rename',
                          action="store_true",
                          help="skip initial contig renaming.")
    optional.add_argument("--elements", dest='elements',
                          action="store", nargs="*",
                          default=["plasmids", "islands", "prophages", "is"],
                          choices=["plasmids", "islands", "prophages", "is"],
                          help="which regions to look for.")
    optional.add_argument("--analyses", dest='analyses',
                          action="store", nargs="*",
                          default=["resfinder", "vfdb"],
                          choices=["ncbi", "card", "resfinder",
                                   "argannot",  "vfdb", "ecoli_vf"],
                          help="which analyses to perform on the " +
                          "extracted mobile genome.")
    optional.add_argument("-s", "--restart_stage", dest='restart_stage',
                          choices=[1, 2, 3, 4, 5, 6, 7],
                          help="stage, if restarting:" +
                          "1 - from the begining, run prokka |" +
                          "2 - run prophage finders |" +
                          "3 - run plasmid finders | " +
                          "4 - run genomic island finders | " +
                          "5 - run insertion sequences finder" +
                          "6 - analyze the results",
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



def make_containerized_cmd(args, command, indir=None, outdir=None):
    """
    note that memory must be provided in gigabytes
    """
    if indir is None:
        indir = os.getcwd()
    if outdir is None:
        outdir = os.getcwd()
    if args.virtualization == "docker":
        cmd = str(
            "docker run --rm " +
            "--memory={args.memory}g " +
            "--cpus={args.cores} " +
            "-v {indir}:/input -v {outdir}:/output {command}"
        ).format(**locals())
    else:
        cmd = str(
            "singularity exec --rm -B " +
        "{indir}:/input -v {outdir}:/output {command}").format(**locals())
    return cmd

def parse_prophet_results(results):
    results_text = []
    templated = []
    with open(results) as inf:
        for line in inf:
            results = ["prophet", "prophage"]
            results.extend(line.strip().split("\t"))
            results_text.append(results)
    for line in results_text:
        subresults = []
        for i in [0, 1, 2, 4, 5]:
            subresults.append(line[i])
        templated.append(subresults)
    return templated


def parse_dimob_results(results):
    results_text = []
    templated = []
    with open(results) as inf:
        for line in inf:
            results = ["dimob", "island"]
            results.extend(line.strip().split("\t"))
            results_text.append(results)
    for line in results_text:
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
            results = ["mlplasmids", "plasmid", "1"]
            results.extend(line.strip().split("\t"))
            mlplasmids_results_text.append(results)
    for line in mlplasmids_results_text:
        if line[5] == '"Plasmid"':
            subresults = []
            for i in [0, 1, 6, 2, 7]:
                # if i == 2:
                #     subresults.append(line[i])
                # else:
                # mlplasmids quotes contig name
                #  I should really talk to whoever wrote that run script.
                subresults.append(line[i].replace('"', ""))
            templated.append(subresults)
    return templated

def condensce_regions(all_results):
    merged_labeled = []
    for seq in set([x[2] for x in all_results]):
        intervals = [(int(x[3]), int(x[4])) for x in all_results if x[2] == seq]
        sorted_by_lower_bound = sorted(intervals, key=lambda tup: tup[0])
        for i in intervals:
            print(i)
        merged = []
        for higher in sorted_by_lower_bound:
            if not merged:
                merged.append(higher)
            else:
                lower = merged[-1]
                # test for intersection between lower and higher:
                # we know via sorting that lower[0] <= higher[0]
                if higher[0] <= lower[1]:
                    upper_bound = max(lower[1], higher[1])
                    merged[-1] = (lower[0], upper_bound)  # replace by merged interval
                else:
                    merged.append(higher)
        for i in merged:
            merged_labeled.append(["combined", "combined", seq, i[0], i[1]])
    return merged_labeled


def annotate_overlaps():
    for prog, feat, seq, start, end in all_results:
        pass


def write_sequence_regions_of_interest(contigs, output_path,  all_results):
    """write out all regions of interest to a single fasta file
    """
    start_stop_dict = [] # orig
    total_length = 0
    with open(contigs, "r") as inf, open(output_path, "w") as outf:
        for rec in SeqIO.parse(inf, "fasta"):
            start_stop_dict[rec.id] = [len(rec.seq)]
            this_length = 0
            for program, type, recid, start, stop in all_results:
                these_features = []
                if rec.id == recid:
                    start, stop = int(start), int(stop)
                    header = \
                        "lcl|{recid}-{start}-{stop}|{program}-{type}".format(
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
                    these_features.append(
                        this_length += (stop-start)
    print("wrote out %i bases" % total_length)
    return total_length


def write_out_names_key(inA, inB, outfile):
    contig_key = []
    inA_names = []
    inB_names = []
    with open(inA, "r") as inAf:
        for rec in SeqIO.parse(inAf, "fasta"):
            inA_names.append([len(rec.seq), rec.id])
    with open(inB, "r") as inBf:
        for rec in SeqIO.parse(inBf, "fasta"):
            inB_names.append([len(rec.seq), rec.id])
    # I hope this never happens
    assert len(inA_names) == len(inB_names), \
        "length of fasta before and after prokka is different!"
    with open(outfile, "w") as outf:
        outf.write("original_name\toriginal_length\tnew_name\tnew_length\n")
        for a, b  in zip(sorted(inA_names, reverse=True), \
                         sorted(inB_names, reverse=True)):
            outf.write("{}\t{}\t{}\t{}\n".format(
                a[0], a[1], b[0], b[1]))


def run_annotation(args, images_dict):
    if not args.skip_rename:
        dest_fasta = os.path.join(args.output, "new_fasta.fasta")
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
    prokka_cmd = make_containerized_cmd(
        args=args,
        command=str(
            "{image} {exe} " +
            "/input/{infile} -o /output/{outdir} --prefix {name} " +
            "--fast --cpus {cpus} 2> {log}").format(
                image=images_dict['prokka']["image"],
                exe=images_dict['prokka']["exe"],
                infile=os.path.relpath(args.contigs),
                outdir=os.path.relpath(prokka_dir),
                name=args.name,
                cpus=args.cores,
                log=os.path.join(args.output, "log.txt"))
    )
    print(prokka_cmd)
    subprocess.run([prokka_cmd],
                   shell=sys.platform != "win32",
                   stdout=subprocess.PIPE,
                   stderr=subprocess.PIPE,
                   check=True)


def run_prophet(args, prokka, prophet_dir, images_dict):
    if os.path.exists(prophet_dir):
        shutil.rmtree(prophet_dir)
    prophet_cmd = make_containerized_cmd(
        args=args,
        command=str(
            "{image} {exe} " +
            "--fasta_in /input/{infilefasta} --gff_in /input/{infilegff} " +
            "--outdir /output/{outdir}/  2> {log}").format(
                image=images_dict['prophet']['image'],
                exe=images_dict['prophet']['exe'],
                infilefasta=os.path.relpath(prokka.fna),
                infilegff=os.path.relpath(prokka.gff),
                outdir=os.path.relpath(prophet_dir),
                log=os.path.join(args.output, "PropheET_log.txt"))
    )
    print(prophet_cmd)
    subprocess.run([prophet_cmd],
                   shell=sys.platform != "win32",
                   stdout=subprocess.PIPE,
                   stderr=subprocess.PIPE,
                   check=True)

def run_mlplasmids(args, prokka, mlplasmids_results, images_dict):
    if os.path.exists(mlplasmids_results):
        os.remove(mlplasmids_results)
    mlplasmids_cmd = make_containerized_cmd(
        args=args,
        command=str(
            "{image} {exe} " +
            "/input/{infilefasta} /output/{outdir}  " +
            ".8 'Escherichia coli' 2> {log}").format(
                image=images_dict['mlplasmids']['image'],
                exe=images_dict['mlplasmids']['exe'],
                infilefasta=os.path.relpath(prokka.fna),
                outdir=os.path.relpath(mlplasmids_results),
                log=os.path.join(args.output, "mlplasmids_log.txt"))
    )
    print(mlplasmids_cmd)
    subprocess.run([mlplasmids_cmd],
                   shell=sys.platform != "win32",
                   stdout=subprocess.PIPE,
                   stderr=subprocess.PIPE,
                   check=True)

def run_dimob(args, prokka, island_results, images_dict):
    island_dir = os.path.dirname(island_results)
    island_seqs_dir = os.path.join(island_dir, "tmp")
    island_results_dir = os.path.join(island_dir, "results")
    if os.path.exists(island_dir):
        shutil.rmtree(island_dir)
    for path in [island_dir, island_seqs_dir, island_results_dir]:
        os.makedirs(path)
    # dimob doesnt process assemblies; we write out each sequence to a tmp file prior to processing
    # {id: {gbk: path_to_gbk, results: path_to_results}}
    island_path_results = {}
    with open (prokka.gbk, "r") as inf:
        for i, rec in enumerate(SeqIO.parse(inf, "genbank")):
            gbk = os.path.join(island_seqs_dir, "%i.gbk" %i)
            result = os.path.join(island_results_dir, "results_%s" %i)
            # dimob wont work if number of CDS's is 0
            if rec.features:
                ncds = sum([1 for feat in rec.features if feat.type == "CDS"])
                if ncds > 0:
                    island_path_results[rec.id] = {"gbk": gbk,
                                                   "result": result}
                    SeqIO.write(rec, gbk, "genbank")
                else:
                    print(str(
                        "Note: contig {0} does not have any CDSs, " +
                        "and will not be analyzed with Dimob").format(rec.id))

    for k, v  in island_path_results.items():
        # Note that dimob does not have an exe
        island_cmd = make_containerized_cmd(
            args=args,
            command=str(
                "{image} /input/{infile} /output/{outdir}" ).format(
                    image=images_dict['dimob']['image'],
                    infile=os.path.relpath(v["gbk"]),
                    outdir=os.path.relpath(v['result']),
                    log=os.path.join(args.output, "dimob_log.txt"))
        )
        print(island_cmd)
        subprocess.run([island_cmd],
                       shell=sys.platform != "win32",
                       stdout=subprocess.PIPE,
                       stderr=subprocess.PIPE,
                       check=True)
        # modify output to include sequence name,
        #   append to final results file
        with open(v['result'], "r") as inf, open(island_results, "a") as outf:
            for line in inf:
                outf.write(k + "\t" + line)


def run_abricate(args, abricate_dir, mobile_fasta, images_dict):
    # remove old results
    if os.path.exists(abricate_dir):
        shutil.rmtree(abricate_dir)
    os.makedirs(abricate_dir)
    choices=["ncbi", "card", "resfinder",
             "argannot",  "vfdb", "ecoli_vf"],
    cmds = []
    for db in args.analyses:
        if db not in choices:
            next
        # write to output stdout
        this_cmd = make_containerized_cmd(
            args=args,
            command=str(
                "{image} {exe} " +
                "--db {db}  /input/{infilefasta} > {outdir}/{db}.tab  2> {log}_{db}_log.txt").format(
                    image=images_dict['abricate']['image'],
                    db=db,
                    exe=images_dict['abricate']['exe'],
                    infilefasta=os.path.relpath(mobile_fasta),
                    outdir=abricate_dir,
                    log=os.path.join(args.output, "abricate"))
        )
        cmds.append(this_cmd)
    for cmd in cmds:
        print(cmd)
        subprocess.run([cmd],
                   shell=sys.platform != "win32",
                   stdout=subprocess.PIPE,
                   stderr=subprocess.PIPE,
                   check=True)

def make_cgview_tab_file(args, seqlen, prokka_fna, reference_fasta,  results_list):
    """
    http://wishart.biology.ualberta.ca/cgview/tab_input.html
    """
    results = []
    header=[
        "#{}".format(args.name),
        "%{}".format(seqlen),
        "!strand\tslot\tstart\tstop\ttype\tlabel\tmouseover\thyperlink"
    ]
    entries = []
    start_stop_dict = {} # seq=[orig_len, new_start, new_end]
    with open(reference_fasta, "r") as inf:
        global_start, global_end = 0, 0
        for rec in SeqIO.parse(inf, "fasta"):
            start_stop_dict[rec.id] = [len(rec.seq)]

    with open(reference_fasta, "r") as inf:
        global_start, global_end = 0, 0
        for rec in SeqIO.parse(inf, "fasta"):
            global_start = global_end + 1
            global_end = global_start + len(rec.seq)
            name = rec.id.split("|")[1].split("-")[0]
            entries.append(
                "\t".join(
                    ["forward", "1", str(global_start), str(global_end), "open_reading_frame",
                     name, rec.id, "github.com/nickp60/happie"]))
            try:
                start_stop_dict[name][2] += len(rec.seq)
            except KeyError:
                    print("error matching names!")
                    sys.exit()
            except IndexError:
                    start_stop_dict[name].extend([global_start, global_end])
            except Exception as e:
                    raise(e)

    for i, l in enumerate(results_list):
        seqN = i + 1
        for program, region, seq, start, end in l:
            this_start = start_stop_dict[seq][1] + int(start)
            this_end = start_stop_dict[seq][0] + int(start) + int(end)
            entries.append(
                "\t".join(
                    ["forward", str(seqN + 1),
                     str(this_start),
                     str(this_end),
                     "predicted_gene",
                     program, region, "github.com/nickp60/happie"]))
    results.extend(header)
    results.extend(entries)
    return results


def main(args=None):
    if args is None:
        args = get_args()
    args.output = os.path.abspath(os.path.expanduser(args.output))
    try:
        if args.restart_stage == 1:
            os.makedirs(args.output, exist_ok=False)
    except FileExistsError:
        print("Output directory already exsists! "+
              "Please chose other desination or, " +
              "if restarting previous analysis, set --stages 2 or above")
        sys.exit(1)
    # make output dir names
    config_data = sm.get_containers_manifest()
    images_dict = sm.parse_docker_images(config_data)
    prokka_dir = os.path.join(args.output, "prokka")
    prophet_dir = os.path.join(args.output, "ProphET", "")
    prophet_results = os.path.join(prophet_dir, "phages_coords")
    mobsuite_dir = os.path.join(args.output, "mobsuite", "")
    island_dir = os.path.join(args.output, "dimob", "")
    island_results = os.path.join(args.output, "dimob", "dimob_results")
    mlplasmids_dir = os.path.join(args.output, "mlplasmids", "")
    mlplasmids_results = os.path.join(args.output, "mlplasmids", "results.txt")
    abricate_dir = os.path.join(args.output, "abricate", "")
    cafe_results = os.path.join(args.output, "cafe", "results.txt")
    test_exes(exes=[args.virtualization])
    # make sub directories.  We don't care if they already exist;
    #  cause we clobber them later if they will cause problems for reexecution
    for path in [prophet_dir, mobsuite_dir, mlplasmids_dir, island_dir,
                 abricate_dir]:
        os.makedirs(path, exist_ok=True)

    # This causes problems on HPCs when I forget to set cores:(
    #    I set a default in the args
    # if args.cores is None:
    #     args.cores = multiprocessing.cpu_count()
    if args.restart_stage < 2:
        run_annotation(args, images_dict),
    else:
        prokka_dir = os.path.abspath(os.path.expanduser(prokka_dir))
        example_prokka_for_name_parsing = \
            glob.glob(os.path.join(prokka_dir, "*.fna"))[0]
        args.name = os.path.splitext(
            os.path.basename(example_prokka_for_name_parsing))[0]

    # set some names of shtuff
    try:
        prokka_fna = glob.glob(os.path.join(prokka_dir, "*.fna"))[0]
        prokka_gbk = glob.glob(os.path.join(prokka_dir, "*.gbk"))[0]
        prokka_gff = glob.glob(os.path.join(prokka_dir, "*.gff"))[0]
    except IndexError:
        raise FileNotFoundError("File not found - something went wrong in step 1 with prokka")

    prokka = Namespace(
        fna=prokka_fna,
        gbk=prokka_gbk,
        gff=prokka_gff,)
    ##########################  finding mobile element ################################3
    if args.restart_stage < 3 and any([x=="prophages" for x in args.elements]):
        run_prophet(args, prokka, prophet_dir, images_dict)
    if args.restart_stage < 4 and any([x=="plasmids" for x in args.elements]):
        run_mlplasmids(args, prokka, mlplasmids_results, images_dict)
    if args.restart_stage < 5 and any([x=="islands" for x in args.elements]):
        run_dimob(args, prokka, island_results, images_dict)
    if args.restart_stage < 6 and any([x=="is" for x in args.elements]):
        pass
        #run_dimob(args, prokka, island_results, images_dict)


    ###########################################################################
    # program type sequence start end
    all_results = []
    results_list = []
    if any([x=="prophages" for x in args.elements]):
        prophet_parsed_result = parse_prophet_results(prophet_results)
        all_results.extend(prophet_parsed_result)
        results_list.append(prophet_parsed_result)
    if any([x=="plasmids" for x in args.elements]):
        mlplasmids_parsed_result = parse_mlplasmids_results(mlplasmids_results)
        all_results.extend(mlplasmids_parsed_result)
        results_list.append(mlplasmids_parsed_result)
    if any([x=="islands" for x in args.elements]):
        dimob_parsed_result = parse_dimob_results(island_results)
        all_results.extend(dimob_parsed_result)
        results_list.append(dimob_parsed_result)

    print(all_results)
    non_overlapping_results = condensce_regions(all_results)
    reference_mobile_genome_path = os.path.join(args.output, "reference_mobile_genome.fasta")
    mobile_genome_path = os.path.join(args.output, "total_mobile_genome.fasta")
    output_regions = os.path.join(args.output, "mobile_genome_coords")
    output_key_path = os.path.join(args.output, "names_key")
    redundant_length = write_sequence_regions_of_interest(
        contigs=prokka_fna,
        output_path=mobile_genome_path,
        all_results=all_results)
    non_redundant_length = write_sequence_regions_of_interest(
        contigs=prokka_fna,
        output_path=reference_mobile_genome_path,
        all_results=non_overlapping_results)
    write_out_names_key(inA=args.contigs, inB=prokka_fna,
                        outfile=output_key_path)
    with open(output_regions, "w") as outf:
        for line in all_results:
            outf.write("\t".join(line) + "\n")

    tab_data = make_cgview_tab_file(args, results_list=results_list,
                                    prokka_fna=prokka.fna,
                                    reference_fasta=reference_mobile_genome_path,
                                    seqlen=non_redundant_length)
    with open("tmp.tab", "w") as outf:
        for line in tab_data:
            outf.write(line + "\n")
    #########################################
    # run_abricate(args, abricate_dir, mobile_fasta=reference_mobile_genome_path,
    #              images_dict=images_dict)

if __name__ == "__main__":
    args = get_args()
    main(args)
