
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
import random
import string
import subprocess
import glob
import pkg_resources
import yaml
from argparse import Namespace
from Bio.Seq import Seq
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
from Bio.SeqFeature import SeqFeature, FeatureLocation

from . import __version__
from . import shared_methods as sm

def get_args():  # pragma: no cover
    parser = argparse.ArgumentParser(
        description="extract the mobile elements for pangenome analysis;" +
        "if running QC, check default values, as those are geared towards " +
        "E. coli <https://enterobase.readthedocs.io/en/latest/pipelines/backend-pipeline-qaevaluation.html>",
        add_help=False)
    # contigs not needed if just re-analyzing the results
    parser.add_argument("--contigs", action="store",
                        help="FASTA or Genbank formatted genome or set of " +
                        "contigs"
    )
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
    optional.add_argument("--images_dir", dest='images_dir',
                          help="if using singularity, where to store your "
                          +"singularity images",
                          default=sm.get_happie_dir())
    optional.add_argument("-n", "--name", dest='name',
                          help="name of experiment; defaults to file name " +
                          "without extension.", default="test")
    optional.add_argument("--skip_QC", dest='skip_QC',
                          action="store_true",
                          help="skip initial contig QC.")
    optional.add_argument("--skip_rename", dest='skip_rename',
                          action="store_true",
                          help="skip initial contig renaming.")
    optional.add_argument("--skip_reannotate", dest='skip_reannotate',
                          action="store_true",
                          help="skip re-annotation of mobile genome with prokka.")
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
                          "6 - analyze the results " +
                          "7 - run annofilt",
                          type=int,
                          default=1)
    optional.add_argument("--skip_annofilt", dest="skip_annofilt",
                          action="store_true",
                          help="skip annotation filtering mobile genome with annofilt.")
    optional.add_argument("--annofilt_reference", dest='annofilt_reference',
                          help="if running annofilt, set this to the " +
                          "species-specific pangnome for annotation " +
                          "filtering")
    optional.add_argument("--QC_min_assembly", dest='QC_min_assembly',
                          default=3700000,
                          type=int,
                          help="if running QC, minimum total assembly length; default: %(default)s")
    optional.add_argument("--QC_max_assembly", dest='QC_max_assembly',
                          default=6400000,
                          type=int,
                          help="if running QC, maximum total assembly length; default: %(default)s")
    optional.add_argument("--QC_min_contig", dest='QC_min_contig',
                          default=800,
                          type=int,
                          help="if running QC, minimum contig length; default: %(default)s")
    optional.add_argument("--QC_min_cov", dest='QC_min_cov',
                          default=0.2,
                          type=float,
                          help="if running QC, minimum coverage fraction of mean coverage; default: %(default)s")
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
            raise ValueError("%s executable not found" % exe)



def make_containerized_cmd(args, image, dcommand, scommand, indir=None, outdir=None, sing=None):
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
            "-v {indir}:/input " +
            "-v {outdir}:/output " +
            "{image} {dcommand}"
        ).format(**locals())
    else:
        cmd = str(
            "{args.images_dir}{sing} {scommand}").format(**locals())
    return cmd

def parse_prophet_results(results):
    results_text = []
    templated = []
    with open(results) as inf:
        for line in inf:
            results = ["prophet", "prophages"]
            results.extend(line.strip().split("\t"))
            results_text.append(results)
    for line in results_text:
        subresults = []
        for i in [0, 1, 2, 4, 5]:
            if i in [4,5]:
                subresults.append(int(line[i]))
            else:
                subresults.append(line[i])
        templated.append(subresults)
    return templated


def parse_dimob_results(results):
    results_text = []
    templated = []
    with open(results) as inf:
        for line in inf:
            results = ["dimob", "islands"]
            results.extend(line.strip().split("\t"))
            results_text.append(results)
    for line in results_text:
        subresults = []
        for i in [0, 1, 2, 4, 5]:
            if i in [4,5]:
                subresults.append(int(line[i]))
            else:
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
            results = ["mlplasmids", "plasmids", "1"]
            results.extend(line.strip().split("\t"))
            mlplasmids_results_text.append(results)
    for line in mlplasmids_results_text:
        line = [x.replace('"', '') for x in line]
        if line[5] == 'Plasmid':
            subresults = []
            for i in [0, 1, 6, 2, 7]:
                # if i == 2:
                #     subresults.append(line[i])
                # else:
                # mlplasmids quotes contig name
                #  I should really talk to whoever wrote that run script.
                if i in [2,7]:
                    subresults.append(int(line[i]))
                else:
                    subresults.append(line[i])
                # subresults.append(line[i].replace('"', ""))
            templated.append(subresults)
            # print(subresults)
    return templated

def condensce_regions(all_results):
    merged_labeled = []
    for seq in set([x[2] for x in all_results]):
        intervals = [(int(x[3]), int(x[4])) for x in all_results if x[2] == seq]
        sorted_by_lower_bound = sorted(intervals, key=lambda tup: tup[0])
        # for i in intervals:
        #     print(i)
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


def write_annotated_mobile_genome(contigs,seed, output_path,  all_results, non_overlapping_results):
    """write out all regions of interest to a single fasta file
    """
    # we want to make our hit list along side this for cgview
    # ring 1 is the "genome"
    # additional rings are for possible element/annoations
    # {seq: {length, global_start, mobile_start, }}
                    # entries.append(
                    #     "\t".join(
                    #         ["forward",
                    #          str(i+2),
                    #          str(rel_start),
                    #          str(rel_end ),
                    #          "gene",
                    #          typ,
                    #          program,
                    #          "github.com/nickp60/happie"]))
    cgview_entries= []
    total_length = 0
    programs = set([x[0] for x in all_results ])
    # start at ring 2; righ 1 is the base
    program_rings = dict(zip(programs, ([x + 2 for x in range(len(programs))])))
    # make prefix based on seed
    random.seed(seed)
    locus_prefix = ''.join(random.choice(string.ascii_uppercase) for x in range(6))
    locus_start, locus_increment = 5, 5
    with open(contigs, "r") as inf, open(output_path + ".fasta", "w") as outfasta, \
         open(output_path + ".gbk", "w") as outgbk:
        for rec in SeqIO.parse(inf, "fasta"):
            ring = 1
            entry_start = total_length + 1
            regions_subset = [(x[3], x[4]) for x in non_overlapping_results if x[2] == rec.id]
            if not regions_subset:
                continue
            seq = ""
            for region_start, region_end in regions_subset:
                seq = seq + str(rec.seq[region_start - 1: region_end - 1])
            cgview_entries.append(
                "\t".join(
                    ["forward",
                     str(ring),
                     str(entry_start),
                     str(entry_start + len(seq)),
                     "open_reading_frame",
                     rec.id,
                     "happie",
                     "github.com/nickp60/happie"]))

            seqrec = SeqRecord(
                Seq(str(seq), IUPAC.unambiguous_dna),
                id=rec.id,  # should we differentiate?
                name='happie_elements',
                annotations={
                    "version": "0.1",
                    "source": "test"
                },
                description='mobile elements from assembly')
            total_length += len(seqrec.seq)
            previous_end = 1
            for program, typ, recid, start, stop in all_results:
                if (rec.id == recid):
                    this_len = stop - start
                    cgview_entries.append(
                        "\t".join(
                            ["forward",
                             str(program_rings[program]),
                             str(previous_end),
                             str(previous_end + this_len),
                             "open_reading_frame",
                             program,
                             typ,
                             "github.com/nickp60/happie"]))

                    feature = SeqFeature(
                        FeatureLocation(
                            # start=4, end=199),
                            start=previous_end + 1,
                            end=previous_end + 1 + this_len),
                        type='gene',
                        qualifiers={
                            "locus_tag": "{}_{}".format(
                                locus_prefix, locus_start),
                            "program": program,
                            "product": typ}
                    )
                    locus_start += locus_increment
                    previous_end = previous_end + this_len
                    seqrec.features.append(feature)
            SeqIO.write(seqrec, outfasta, "fasta")
            SeqIO.write(seqrec, outgbk, "genbank")
    print("wrote out %i bases" % total_length)
    return(total_length, cgview_entries)


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


def run_annotation(args, contigs, prokka_dir, images_dict, skip_rename=True,
                   new_name="new_fasta.fasta", log_dir=None):
    if os.path.exists(prokka_dir):
        shutil.rmtree(prokka_dir)
    if not skip_rename:
        dest_fasta = os.path.join(args.output, new_name)
        with open(contigs, "r") as inf, open(dest_fasta, "w") as outf:
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
        contigs = dest_fasta
    print("running prokka")
    prokka_cmd = make_containerized_cmd(
        args=args,
        image=images_dict['prokka']["image"],
        sing=images_dict['prokka']["sing"],
        dcommand=str(
            "/input/{infile} -o /output/{outdir} --prefix {name} " +
            "--fast --cpus {cpus} 2> {log}").format(
                infile=os.path.relpath(contigs),
                outdir=os.path.relpath(prokka_dir),
                name=args.name,
                cpus=args.cores,
                log=os.path.join(log_dir, "prokka.log")),
        scommand=str(
            "{infile} -o {outdir} --prefix {name} " +
            "--fast --cpus {cpus} 2> {log}").format(
                infile=contigs,
                outdir=prokka_dir,
                name=args.name,
                cpus=args.cores,
                log=os.path.join(log_dir, "prokka.log"))

    )
    print(prokka_cmd)
    subprocess.run([prokka_cmd],
                   shell=sys.platform != "win32",
                   stdout=subprocess.PIPE,
                   stderr=subprocess.PIPE,
                   check=True)


def run_annofilt(args, annofilt_dir, prokka_dir, images_dict, log_dir=None):
    if os.path.exists(annofilt_dir):
        print("removing old annofilt dir")
        shutil.rmtree(annofilt_dir)
    print("running annofilt")
    annofilt_cmd = make_containerized_cmd(
        args=args,
        image=images_dict['annofilt']["image"],
        sing=images_dict['annofilt']["sing"],
        dcommand=str(
            "/input/{ref} /input/{infile} -o /output/{outdir} " +
            " --threads {cpus} 2> {log}").format(
                ref=os.path.relpath(args.annofilt_reference),
                infile=os.path.relpath(prokka_dir),
                outdir=os.path.relpath(annofilt_dir),
                cpus=args.cores,
                log=os.path.join(log_dir, "annofilt.log")),
        scommand=str(
            "{ref} {infile} -o {outdir} " +
            "--threads {cpus} 2> {log}").format(
                ref=args.annofilt_reference,
                infile=prokka_dir,
                outdir=annofilt_dir,
                cpus=args.cores,
                log=os.path.join(log_dir, "annofilt.log"))

    )
    print(annofilt_cmd)
    subprocess.run([annofilt_cmd],
                   shell=sys.platform != "win32",
                   stdout=subprocess.PIPE,
                   stderr=subprocess.PIPE,
                   check=True)


def run_prophet(args, prokka, prophet_dir, images_dict, log_dir=None):
    if os.path.exists(prophet_dir):
        shutil.rmtree(prophet_dir)
    prophet_cmd = make_containerized_cmd(
        args=args,
        image=images_dict['prophet']['image'],
        sing=images_dict['prophet']["sing"],
        dcommand=str(
            "--fasta_in /input/{infilefasta} --gff_in /input/{infilegff} " +
            "--outdir /output/{outdir}/ --cores {cores} 2> {log}").format(
                infilefasta=os.path.relpath(prokka.fna),
                infilegff=os.path.relpath(prokka.gff),
                outdir=os.path.relpath(prophet_dir),
                cores=args.cores,
                log=os.path.join(log_dir, "PropheET.log")),
        scommand=str(
            "--fasta_in {infilefasta} --gff_in {infilegff} " +
            "--outdir {outdir}/ --cores {cores} 2> {log}").format(
                infilefasta=os.path.relpath(prokka.fna),
                infilegff=os.path.relpath(prokka.gff),
                outdir=os.path.relpath(prophet_dir),
                cores=args.cores,
                log=os.path.join(log_dir, "PropheET.log")),
    )
    print(prophet_cmd)
    subprocess.run([prophet_cmd],
                   shell=sys.platform != "win32",
                   stdout=subprocess.PIPE,
                   stderr=subprocess.PIPE,
                   check=True)


def run_mlplasmids(args, prokka, mlplasmids_results, images_dict, log_dir=None):
    if os.path.exists(mlplasmids_results):
        os.remove(mlplasmids_results)
    mlplasmids_cmd = make_containerized_cmd(
        args=args,
                image=images_dict['mlplasmids']['image'],
        sing=images_dict['mlplasmids']["sing"],
        dcommand=str(
            "/input/{infilefasta} /output/{outdir}  " +
            ".8 'Escherichia coli' 2> {log}").format(
                infilefasta=os.path.relpath(prokka.fna),
                outdir=os.path.relpath(mlplasmids_results),
                log=os.path.join(log_dir, "mlplasmids.log")),
        scommand=str(
            "{infilefasta} {outdir}  " +
            ".8 'Escherichia coli' 2> {log}").format(
                infilefasta=prokka.fna,
                outdir=mlplasmids_results,
                log=os.path.join(log_dir, "mlplasmids.log")),


    )
    print(mlplasmids_cmd)
    subprocess.run([mlplasmids_cmd],
                   shell=sys.platform != "win32",
                   stdout=subprocess.PIPE,
                   stderr=subprocess.PIPE,
                   check=True)


def run_dimob(args, prokka, island_results, images_dict, log_dir=None):
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
            sing=images_dict['dimob']["sing"],
            image=images_dict['dimob']['image'],
            dcommand=str(
                "/input/{infile} /output/{outdir}" ).format(
                    infile=os.path.relpath(v["gbk"]),
                    outdir=os.path.relpath(v['result']),
                    log=os.path.join(log_dir, "dimob.log")),
            scommand=str(
                "{infile} {outdir}" ).format(
                    infile=v["gbk"],
                    outdir=v['result'],
                    log=os.path.join(log_dir, "dimob.log")),

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


def run_abricate(args, abricate_dir, mobile_fasta, images_dict, all_results, log_dir=None):
    # remove old results
    if os.path.exists(abricate_dir):
        shutil.rmtree(abricate_dir)
    os.makedirs(abricate_dir)
    cmds = []
    outfiles = []
    for db in args.analyses:
        print(db)
        outfiles.append("{outdir}/{db}.tab".format(outdir=abricate_dir, db=db))
        this_cmd = make_containerized_cmd(
            args=args,
            image=images_dict['abricate']['image'],
            sing=images_dict['abricate']["sing"],
            dcommand=str(
                "--db {db}  /input/{infilefasta} > {outdir}/{db}.tab  " +
                "2> {log}_{db}_log.txt").format(
                    db=db,
                    infilefasta=os.path.relpath(mobile_fasta),
                    outdir=abricate_dir,
                    log=os.path.join(log_dir, "abricate")),
            scommand=str(
                "--db {db}  {infilefasta} > {outdir}/{db}.tab  " +
                "2> {log}_{db}_log.txt").format(
                    db=db,
                    infilefasta=mobile_fasta,
                    outdir=abricate_dir,
                    log=os.path.join(log_dir, "abricate")),


        )
        cmds.append(this_cmd)
    for cmd in cmds:
        print(cmd)
        subprocess.run([cmd],
                       shell=sys.platform != "win32",
                       stdout=subprocess.PIPE,
                       stderr=subprocess.PIPE,
                       check=True)

    with open(all_results, "w") as outf:
        for f in outfiles:
            for line in open(f, "r"):
                outf.write(line)




def make_cgview_tab_file(args, seqlen, cgview_entries):
    """
    http://wishart.biology.ualberta.ca/cgview/tab_input.html
    """
    results = []
    header=[
        "#{}".format(args.name),
        "%{}".format(seqlen + 1),
        "!strand\tslot\tstart\tstop\ttype\tlabel\tmouseover\thyperlink"
    ]
    results.extend(header)
    results.extend(cgview_entries)
    return results


def run_cgview(args, cgview_tab, cgview_dir, images_dict, log_dir=None):
    if os.path.exists(cgview_dir):
        shutil.rmtree(cgview_dir)
    os.makedirs(cgview_dir)
    cgview_cmd = make_containerized_cmd(
        args=args,
        image=images_dict['cgview']['image'],
        sing=images_dict['cgview']["sing"],
        dcommand=str(
            "-i /input/{infilefasta} -f svg -o /output/{outdir} -I T").format(
                exe=images_dict['cgview']['exe'],
                infilefasta=os.path.relpath(cgview_tab),
                outdir=os.path.join(os.path.relpath(cgview_dir), "cgview.svg"),
                log=os.path.join(log_dir, "cgview.log")),
        scommand=str(
            "-i {infilefasta} -f svg -o {outdir} -I T").format(
                exe=images_dict['cgview']['exe'],
                infilefasta=cgview_tab,
                outdir=os.path.join(cgview_dir, "cgview.svg"),
                log=os.path.join(log_dir, "cgview.log")),

    )
    print(cgview_cmd)
    subprocess.run([cgview_cmd],
                   shell=sys.platform != "win32",
                   stdout=subprocess.PIPE,
                   stderr=subprocess.PIPE,
                   check=True)


def make_circleator():
    """

    type - either the name of a predefined track type OR the keyword new
    name - a name by which the track may be referenced from elsewhere in the configuration file
    glyph - the Circleator “glyph” used to render this track
    heightf - the height of the track as a fraction of the circle’s radius (0-1)
    innerf - position of the innermost part of the track as a fraction of the circle’s radius (0-1)
    outerf - position of the outermost part of the track as a fraction of the circle’s radius (0-1)
    data - path to a data file, if one is required by the chosen track type and/or glyph
    feat_type - display only features of the specified type (e.g., “gene”, “tRNA”)
    feat_strand - display only features on the specified strand (e.g., “-“, “+”, “-1”, “1”)
    color1 - interpretation depends on the track type: usually the SVG fill color
    color2 - interpretation depends on the track type: usually the SVG stroke (outline) color
    opacity - opacity of the track between 0 and 1, where 0 = invisible/completely transparent and 1=completely opaque
    zindex - integer z-index of the track: tracks with higher z-indexes are drawn on top of those with lower z-indexes
    options - a comma-delimited list of track options in the format “key=value”


    """
    pass


def write_out_yaml_args(args, outpath):
    if os.path.exists(outpath):
        new_path = outpath + ".bak"
        shutil.move(outpath, new_path)
    with open(outpath, "w") as outf:
        yaml.dump(args, outf)

def read_in_yaml_args(outpath):
    if not os.path.exists(outpath):
        raise FileNotFoundError("cannot open file %s" % outpath)
    with open(outpath, "r") as outf:
        new_args = yaml.load( outf)
    return new_args


def recheck_required_args(args):
    if args.contigs is None:
        raise ValueError("no --contigs provided! see 'happie -h'")
    if not args.skip_annofilt:
        if args.annofilt_reference is None:
            raise ValueError("If using annofilt, you must provide a species specific pangenome from whole genomes for comparison. See the annofilt documentaion for more details. <https://github.com/nickp60/annofilt/>")


def coords_to_merged_gff(coords):
    all_coords = [x[3] for x in coords]
    all_coords.expand([x[4] for x in coords])
    mn, mx = min(all_coords), max(all_coords)
    pass


def QC_bug(args, min_length, max_length, cov_threshold=.2, min_contig_length=800):
    # check assembly size
    log_strings = []
    lengths = []
    with open(args.contigs) as inf:
        for rec in SeqIO.parse(inf, "fasta"):
            lengths.append([rec.id, len(rec.seq)])
    total_length = sum([x[1] for x in lengths])
    log_strings.append("QC criteria -- see happy_args.yaml")
    log_strings.append("QC'd " + str(len(lengths)) + " sequences")
    log_strings.append("Assembly is" + str(total_length) + "bases")
    if not min_length < total_length < max_length:
        with open(os.path.join(args.output, "sublogs", "QC.log"), "w") as logoutf:
            for s in log_strings:
                logoutf.write(s + "\n")
        raise ValueError(
            "Assembly length {} falls outside of QC range of {} to {}".format(
                total_length, min_length, max_length)
    )
    short_contigs = {x[0]: x[1] for x in lengths if
                       x[1] < min_contig_length}
    log_strings.append("Removed " + str(len(short_contigs)) + " short contigs ")

    ncontigs = len(lengths)
    # get SPAdes coverage
    low_cov_contigs = {}
    with open(args.contigs) as inf:
        header_info = {}
        for rec in SeqIO.parse(inf, "fasta"):
            if not rec.id.startswith("NODE"):
                log_strings.append("Warning: cannot QC by assembly coverage; " +
                      "happie only parses SPAdes headers")
                print(log_strings[-1])
            else:
                # NODE_1_length_10442_cov_5.92661
                p = re.compile(r'NODE_(?P<node>\d*?)_length_(?P<length>\d*?)_cov_(?P<cov>[\d|\.]*)')
                m = p.search(rec.id)
                header_info[rec.id] = {
                    "length": int(m.group("length")),
                    "cov": float(m.group("cov"))
                }
    if header_info:
        # this gets skipped if we didnt have spades headers
        mean_coverage = sum([y['cov'] for x, y in header_info.items()])/ncontigs
        low_cov_contigs = {x: y for x, y in header_info.items() if
                       y['cov'] < (mean_coverage * cov_threshold)}
    bad_contigs = {**short_contigs, **low_cov_contigs}
    log_strings.append("Removed " + str(len(low_cov_contigs)) +
                       " low covereage contigs ")
    if bad_contigs:
        print("the following contigs have too low loverage or are too short "
              "and will be removed from the analysis")
        print(bad_contigs)
        # make a filtered file
        outfile = os.path.join(
            args.output,
            os.path.basename(os.path.splitext(args.contigs)[0]) + "_filtered.fasta")
        with open(args.contigs, "r") as inf, open(outfile, "w") as outf:
            for rec in SeqIO.parse(inf, "fasta"):
                if not rec.id in bad_contigs.keys():
                    SeqIO.write(rec, outf, "fasta")
        args.contigs = outfile
    with open(os.path.join(args.output, "sublogs", "QC.log"), "w") as logoutf:
        for s in log_strings:
            logoutf.write(s + "\n")
    if len(bad_contigs) >= ncontigs:
        raise ValueError("All of the contigs are filtered out with the current criteria")


def main(args=None):
    if args is None:
        args = get_args()
    args.output = os.path.abspath(os.path.expanduser(args.output))
    test_exes(exes=[args.virtualization])
    if args.restart_stage == 1:
        recheck_required_args(args)
    recheck_required_args(args)
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
    log_dir = os.path.join(args.output, "sublogs")
    prokka_dir = os.path.join(args.output, "prokka")
    mobile_prokka_dir = os.path.join(args.output, "mobile_prokka")
    prophet_dir = os.path.join(args.output, "ProphET", "")
    prophet_results = os.path.join(prophet_dir, "phages_coords")
    mobsuite_dir = os.path.join(args.output, "mobsuite", "")
    island_dir = os.path.join(args.output, "dimob", "")
    island_results = os.path.join(args.output, "dimob", "dimob_results")
    mlplasmids_dir = os.path.join(args.output, "mlplasmids", "")
    mlplasmids_results = os.path.join(args.output, "mlplasmids", "results.txt")
    abricate_dir = os.path.join(args.output, "abricate", "")
    wgs_abricate_dir = os.path.join(args.output, "abricate_wgs", "")
    annofilt_dir = os.path.join(args.output, "annofilt", "")
    wgs_annofilt_dir = os.path.join(args.output, "annofilt_wgs", "")
    cgview_dir = os.path.join(args.output, "cgview", "")
    # make sub directories.  We don't care if they already exist;
    #  cause we clobber them later if they will cause problems for reexecution
    # except dont make prokka dirs
    for path in [prophet_dir, mobsuite_dir, mlplasmids_dir, island_dir,
                 abricate_dir, wgs_abricate_dir, log_dir]:
        os.makedirs(path, exist_ok=True)

    # write out args for easier re-running
    write_out_yaml_args(args, outpath=os.path.join(args.output, "happie_args.yaml"))

    # This causes problems on HPCs when I forget to set cores:(
    #    I set a default in the args
    # if args.cores is None:
    #     args.cores = multiprocessing.cpu_count()
    if args.restart_stage < 2:
        isfasta = False
        with open(args.contigs, "r") as inf:
            for line in inf:
                if  line.startswith(">"):
                    isfasta = True
                break
        if not isfasta:
            print("converting input to fasta")
            fasta_version = os.path.join(
                args.output,
                os.path.basename(os.path.splitext(args.contigs)[0] + ".fna")
            )
            with open(args.contigs, "r") as i2, open(fasta_version, "w") as o2:
                for rec in SeqIO.parse(i2, "genbank"):
                    SeqIO.write(rec, o2, "fasta")
            args.contigs = fasta_version
        else:
            print("input is fasta")
        if not args.skip_QC:
            QC_bug(
                args,
                min_length=args.QC_min_assembly,
                max_length=args.QC_max_assembly,
                cov_threshold=args.QC_min_cov,
                min_contig_length=args.QC_min_contig)
        print("running prokka")
        run_annotation(args, contigs=args.contigs, prokka_dir=prokka_dir,
                       images_dict=images_dict, skip_rename=args.skip_rename,
                       new_name="postQC_input.fasta", log_dir=log_dir)
    else:
        # read in old config file, if it exists. for now we just get the old
        # path to the contigs, so you dont have to remember how exacly you ran
        # the command to reprocess the results
        try:
            old_args = read_in_yaml_args(
                outpath=os.path.join(args.output, "happie_args.yaml")
            )
            args.contigs = old_args.contigs
        except ValueError:
            print("cannot parse old config!")
            try:
                recheck_required_args(args)
            except Exception as e:
                print(e)
                raise ValueError("When rerunning old analyses without a " +
                                 "'happie_args.yaml' file in outdir, all "+
                                 " the required args must be provided")
        prokka_dir = os.path.abspath(os.path.expanduser(prokka_dir))
        print(prokka_dir)
        examples_prokka_for_name_parsing = \
            glob.glob(os.path.join(prokka_dir, "*.fna"))
        print(examples_prokka_for_name_parsing)
        args.name = os.path.splitext(
            os.path.basename(examples_prokka_for_name_parsing[0]))[0]

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
        run_prophet(args, prokka, prophet_dir, images_dict, log_dir=log_dir)
    if args.restart_stage < 4 and any([x=="plasmids" for x in args.elements]):
        run_mlplasmids(args, prokka, mlplasmids_results, images_dict, log_dir=log_dir)
    if args.restart_stage < 5 and any([x=="islands" for x in args.elements]):
        run_dimob(args, prokka, island_results, images_dict, log_dir=log_dir)
    if args.restart_stage < 6 and any([x=="is" for x in args.elements]):
        pass
        #run_dimob(args, prokka, island_results, images_dict)

    if args.restart_stage < 7:
    ###########################################################################
    # program type sequence start end
        all_results = []
        results_list = []
        if any([x=="prophages" for x in args.elements]):
            if os.path.exists(prophet_results) and os.path.getsize(prophet_results) > 0:
                prophet_parsed_result = parse_prophet_results(prophet_results)
                all_results.extend(prophet_parsed_result)
                results_list.append(prophet_parsed_result)
        print(all_results)
        if any([x=="plasmids" for x in args.elements]):
            if os.path.exists(mlplasmids_results) and os.path.getsize(mlplasmids_results) > 0:
                mlplasmids_parsed_result = parse_mlplasmids_results(mlplasmids_results)
                all_results.extend(mlplasmids_parsed_result)
                results_list.append(mlplasmids_parsed_result)
        print(all_results)
        if any([x=="islands" for x in args.elements]):
            if os.path.exists(island_results) and os.path.getsize(island_results) > 0:
                dimob_parsed_result = parse_dimob_results(island_results)
                all_results.extend(dimob_parsed_result)
                results_list.append(dimob_parsed_result)
        # print(all_results)
        non_overlapping_results = condensce_regions(all_results)
        mobile_genome_path_prefix = os.path.join(args.output, "mobile_genome")
        output_regions = os.path.join(args.output, "mobile_genome_coords")
        with open(output_regions, "w") as outf:
            for line in all_results:
                outf.write(args.contigs + "\t" + "\t".join([str(x) for x in line]) + "\n")
        output_key_path = os.path.join(args.output, "names_key")
        write_out_names_key(inA=args.contigs, inB=prokka_fna,
                            outfile=output_key_path)
        seq_length, cgview_entries = write_annotated_mobile_genome(
            seed=os.path.basename(args.contigs),
            contigs=prokka_fna,
            output_path=mobile_genome_path_prefix,
            non_overlapping_results=non_overlapping_results,
            all_results=all_results)
        if seq_length == 0:
            print("WARNING: none of the genome was detected to be mobile")
            return 0

        tab_data = make_cgview_tab_file(args, cgview_entries=cgview_entries,
                                    seqlen=seq_length)
        cgview_data = os.path.join(args.output, "cgview.tab")
        with open(cgview_data, "w") as outf:
            for line in tab_data:
                outf.write(line + "\n")
        #########################################
        # run abricate on both the mobile genome, and the entire sequence, for enrichment comparison
        abricate_data = os.path.join(args.output, "abricate.tab")
        run_abricate(
            args, abricate_dir,
            mobile_fasta=mobile_genome_path_prefix + ".fasta",
            images_dict=images_dict,
            all_results=abricate_data, log_dir=log_dir)
        wgs_abricate_data = os.path.join(
            args.output, "wgs_abricate.tab")
        run_abricate(
            args, wgs_abricate_dir,
            mobile_fasta=args.contigs,
            images_dict=images_dict,
            all_results=wgs_abricate_data, log_dir=log_dir)
        if not args.skip_reannotate:
            run_annotation(
                args, contigs=mobile_genome_path_prefix + ".fasta",
                prokka_dir=mobile_prokka_dir, images_dict=images_dict,
                skip_rename=False,
                new_name="tmp_mobile.fasta", log_dir=log_dir)
    # Run Annofilt on mobile genome and annotated genome
    if not args.skip_annofilt:
        try:
            mobile_prokka_fna = glob.glob(
                os.path.join(mobile_prokka_dir, "*.fna"))[0]
        except IndexError:
            raise FileNotFoundError(
                "File not found - something went " +
                "wrong annotating the mobile genome.")
        # copy pan_genome temporarrity:
        tmp_pangenome = os.path.join(args.output, "pangenome.fasta")
        shutil.copyfile(args.annofilt_reference, tmp_pangenome)
        old_annofilt_reference = args.annofilt_reference
        args.annofilt_reference = tmp_pangenome
        print("Running annofilt to remove truncated CDSs")
        run_annofilt(
            args,
            annofilt_dir=annofilt_dir,
            prokka_dir=mobile_prokka_dir,
            images_dict=images_dict,
            log_dir=log_dir)
        run_annofilt(
            args,
            annofilt_dir=wgs_annofilt_dir,
            prokka_dir=prokka_dir,
            images_dict=images_dict,
            log_dir=log_dir)
        # run_cgview(args, cgview_tab=cgview_data, cgview_dir=cgview_dir, images_dict=images_dict)

if __name__ == "__main__":
    args = get_args()
    main(args)
