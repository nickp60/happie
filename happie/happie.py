#!/usr/bin/env python3
#-*- coding: utf-8 -*-

import re
import argparse
import sys
import shutil
import os
import random
import string
import glob
import yaml
from argparse import Namespace
from Bio.Seq import Seq
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
from Bio.SeqFeature import SeqFeature, FeatureLocation

from . import __version__
from . import shared_methods as sm
from . import runners as runners
from . import parsers as parsers


def get_args():  # pragma: no cover
    parser = argparse.ArgumentParser(
        description="extract the mobile elements for pangenome analysis;" +
        "if running QC, check default values, as those are geared towards " +
        "E. coli <https://enterobase.readthedocs.io/en/latest/pipelines/" +
        "backend-pipeline-qaevaluation.html>",
        add_help=False)
    # contigs not needed if just re-analyzing the results
    parser.add_argument("--contigs", action="store",
                        help="FASTA or Genbank formatted genome or set of " +
                        "contigs")
    parser.add_argument("-o", "--output", action="store",
                        help="destination dir", required=True)
    parser.add_argument("--virtualization", action="store",
                        help="Whether to run with Docker or Singularity",
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
                          help="if using singularity, where to store your " +
                          "singularity images",
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
                          help="skip re-annotation of mobile genome " +
                          "with prokka.")
    optional.add_argument("--elements", dest='elements',
                          action="store", nargs="*",
                          default=["plasmids", "islands", "prophages", "is"],
                          choices=["plasmids", "islands", "prophages", "is"],
                          help="which regions to look for.")
    optional.add_argument("--plasmid_tools", dest='plasmid_tools',
                          action="store", nargs='+',
                          default=["mlplasmids"],
                          choices=["mlplasmids", "mobsuite", "plasflow"],
                          help="Which plasmid finder tool(s) to use. " +
                          "Note: only use mlplasmids if you  are working on " +
                          "E. coli, E. faecium, or K. pneumoniae")
    optional.add_argument("--mlplasmidsdb", dest='mlplasmidsdb',
                          action="store", default="Escherichia coli",
                          choices=["Escherichia coli",
                                   "Klebsiella pneumoniae",
                                   "Enterococcus faecium"],
                          help="Use mob-suite tools for plasmid ID rather "+
                          "than mlplasmids.  Use this option if not " +
                          "working on E. coli, E. faecium, or K "+
                          " pneumoniae ")
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


def write_annotated_mobile_genome(contigs,seed, output_path, all_results, non_overlapping_results):
    """write out all regions of interest to a single fasta file
    we want to make our hit list along side this for cgview
    ring 1 is the "genome"
    additional rings are for possible element/annoations
    {seq: {length, global_start, mobile_start, }}
                    entries.append(
                        "\t".join(
                            ["forward",
                             str(i+2),
                             str(rel_start),
                             str(rel_end ),
                             "gene",
                             typ,
                             program,
                             "github.com/nickp60/happie"]))
    """
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


def QC_bug(args, QC_dir, min_length, max_length, cov_threshold=.2, min_contig_length=800):
    # check assembly size
    log_strings = []
    lengths = []
    with open(args.contigs) as inf:
        for rec in SeqIO.parse(inf, "fasta"):
            lengths.append([rec.id, len(rec.seq)])
    total_length = sum([x[1] for x in lengths])
    log_strings.append("#QC criteria -- see happy_args.yaml")
    log_strings.append("N sequences\t" + str(len(lengths)))
    log_strings.append("Assembly Length\t" + str(total_length))
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
    log_strings.append("N too short\t" + str(len(short_contigs)))

    ncontigs = len(lengths)
    # get SPAdes coverage
    low_cov_contigs = {}
    spades_headers = True
    with open(args.contigs) as inf:
        header_info = {}
        for rec in SeqIO.parse(inf, "fasta"):
            if not rec.id.startswith("NODE"):
                spades_headers = False
                # log_strings.append("Warning: cannot QC by assembly coverage; " +
                #       "happie only parses SPAdes headers")
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
        # this should prrobably get changed to median or something;  shortie contigs have super high cover
        mean_coverage = sum([y['cov'] for x, y in header_info.items()])/ncontigs
        low_cov_contigs = {x: y for x, y in header_info.items() if
                       y['cov'] < (mean_coverage * cov_threshold)}
    bad_contigs = {**short_contigs, **low_cov_contigs}
    if spades_headers:
        log_strings.append("N low coverage\t" + str(len(low_cov_contigs)))
    else:
        log_strings.append("N low coverage\tNA")
    retained_length = 0
    if bad_contigs:
        print("the following contigs have too low loverage or are too short "
              "and will be removed from the analysis")
        print(bad_contigs)
        # make a filtered file
        outfile = os.path.join(
            args.output, QC_dir,
            os.path.basename(os.path.splitext(args.contigs)[0]) + "_postQC.fasta")
        with open(args.contigs, "r") as inf, open(outfile, "w") as outf:
            for rec in SeqIO.parse(inf, "fasta"):
                if not rec.id in bad_contigs.keys():
                    retained_length += len(rec.seq)
                    SeqIO.write(rec, outf, "fasta")
        args.contigs = outfile
    log_strings.append("Filtered assembly length\t" + str(retained_length))
    with open(os.path.join(args.output, "sublogs", "QC.log"), "w") as logoutf:
        for s in log_strings:
            logoutf.write(s + "\n")
    if len(bad_contigs) >= ncontigs:
        raise ValueError("All of the contigs are filtered out with the current criteria")

    log_strings.append("Filtered assembly is " + str(retained_length) + " bases")



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
    QC_dir = os.path.join(args.output, "QC")
    prokka_dir = os.path.join(args.output, "wgs_prokka")
    mobile_prokka_dir = os.path.join(args.output, "mobile_prokka")
    prophet_dir = os.path.join(args.output, "ProphET", "")
    prophet_results = os.path.join(prophet_dir, "phages_coords")
    island_dir = os.path.join(args.output, "dimob", "")
    island_results = os.path.join(args.output, "dimob", "dimob_results")
    mobsuite_dir = os.path.join(args.output, "mobsuite", "")
    mobsuite_results = os.path.join(args.output, "mobsuite", "contig_report.txt")
    plasflow_dir = os.path.join(args.output, "plasflow", "")
    plasflow_results = os.path.join(args.output, "plasflow", "results.txt")
    mlplasmids_dir = os.path.join(args.output, "mlplasmids", "")
    mlplasmids_results = os.path.join(args.output, "mlplasmids", "results.txt")
    abricate_dir = os.path.join(args.output, "mobile_abricate", "")
    wgs_abricate_dir = os.path.join(args.output, "wgs_abricate", "")
    annofilt_dir = os.path.join(args.output, "mobile_annofilt", "")
    wgs_annofilt_dir = os.path.join(args.output, "wgs_annofilt", "")
    cgview_dir = os.path.join(args.output, "cgview", "")
    # make sub directories. We don't care if they already exist;
    #  cause we clobber them later if they will cause problems for reexecution
    # except dont make prokka dirs
    for path in [prophet_dir, mlplasmids_dir, island_dir, mobsuite_dir,
                 plasflow_dir,
                 abricate_dir, wgs_abricate_dir, log_dir, QC_dir]:
        os.makedirs(path, exist_ok=True)

    # write out args for easier re-running
    write_out_yaml_args(args, outpath=os.path.join(args.output, "happie_args.yaml"))

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
                QC_dir,
                min_length=args.QC_min_assembly,
                max_length=args.QC_max_assembly,
                cov_threshold=args.QC_min_cov,
                min_contig_length=args.QC_min_contig)
        print("running prokka")
        runners.run_annotation(args, contigs=args.contigs, prokka_dir=prokka_dir,
                       images_dict=images_dict, skip_rename=args.skip_rename,
                       new_name="renamed_preprokka_input.fasta",
                       subset="wgs", log_dir=log_dir)
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
        runners.run_prophet(args, prokka, prophet_dir, images_dict, subset="mobile",log_dir=log_dir)
    if args.restart_stage < 4 and any([x=="plasmids" for x in args.elements]):
        if "mobsuite" in args.plasmid_tools:
            runners.run_mobsuite(args, prokka, mobsuite_dir, images_dict, subset="mobile", log_dir=log_dir)
        if "mlplasmids" in args.plasmid_tools:
            runners.run_mlplasmids(args, prokka, mlplasmids_results, images_dict, subset="mobile", log_dir=log_dir)
        if "plasflow" in args.plasmid_tools:
            runners.run_plasflow(args, prokka, plasflow_results, images_dict, subset="mobile", log_dir=log_dir)
    if args.restart_stage < 5 and any([x=="islands" for x in args.elements]):
        runners.run_dimob(args, prokka, island_results, images_dict, subset="mobile", log_dir=log_dir)
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
                prophet_parsed_result = parsers.parse_prophet_results(prophet_results)
                all_results.extend(prophet_parsed_result)
                results_list.append(prophet_parsed_result)
        print(all_results)
        if any([x=="plasmids" for x in args.elements]):
            if "mobsuite" in args.plasmid_tools:
                if os.path.exists(mobsuite_results) and os.path.getsize(mobsuite_results) > 0:
                    mobsuite_parsed_result = \
                        parsers.parse_mobsuite_results(mobsuite_results)
                    all_results.extend(mobsuite_parsed_result)
                    results_list.append(mobsuite_parsed_result)
            if "plasflow" in args.plasmid_tools:
                if os.path.exists(plasflow_results) and os.path.getsize(plasflow_results) > 0:
                    plasflow_parsed_result = \
                        parsers.parse_plasflow_results(plasflow_results)
                    all_results.extend(plasflow_parsed_result)
                    results_list.append(plasflow_parsed_result)
            if "mlplasmids" in args.plasmid_tools:
                if os.path.exists(mlplasmids_results) and os.path.getsize(mlplasmids_results) > 0:
                    mlplasmids_parsed_result = parsers.parse_mlplasmids_results(mlplasmids_results)
                    all_results.extend(mlplasmids_parsed_result)
                    results_list.append(mlplasmids_parsed_result)
        print(all_results)
        if any([x=="islands" for x in args.elements]):
            if os.path.exists(island_results) and os.path.getsize(island_results) > 0:
                dimob_parsed_result = parsers.parse_dimob_results(island_results)
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

        tab_data = runners.make_cgview_tab_file(args, cgview_entries=cgview_entries,
                                    seqlen=seq_length)
        cgview_data = os.path.join(args.output, "cgview.tab")
        with open(cgview_data, "w") as outf:
            for line in tab_data:
                outf.write(line + "\n")
        #########################################
        # run abricate on both the mobile genome, and the entire sequence, for enrichment comparison
        abricate_data = os.path.join(args.output, "mobile_abricate.tab")
        runners.run_abricate(
            args, abricate_dir,
            mobile_fasta=mobile_genome_path_prefix + ".fasta",
            images_dict=images_dict,
            all_results=abricate_data,subset="mobile", log_dir=log_dir)
        wgs_abricate_data = os.path.join(
            args.output, "wgs_abricate.tab")
        runners.run_abricate(
            args, wgs_abricate_dir,
            mobile_fasta=args.contigs,
            images_dict=images_dict,
            all_results=wgs_abricate_data, subset="wgs", log_dir=log_dir)
        if not args.skip_reannotate:
            runners.run_annotation(
                args, contigs=mobile_genome_path_prefix + ".fasta",
                prokka_dir=mobile_prokka_dir, images_dict=images_dict,
                skip_rename=True,
                new_name="tmp_mobile.fasta", subset="mobile", log_dir=log_dir)
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
        tmp_pangenome = os.path.join(args.output, "tmp_pangenome.fasta")
        shutil.copyfile(args.annofilt_reference, tmp_pangenome)
        old_annofilt_reference = args.annofilt_reference
        args.annofilt_reference = tmp_pangenome
        print("Running annofilt to remove truncated CDSs")
        runners.run_annofilt(
            args,
            annofilt_dir=annofilt_dir,
            prokka_dir=mobile_prokka_dir,
            images_dict=images_dict,
            subset="mobile",
            log_dir=log_dir)
        runners.run_annofilt(
            args,
            annofilt_dir=wgs_annofilt_dir,
            prokka_dir=prokka_dir,
            images_dict=images_dict,
            subset="wgs",
            log_dir=log_dir)
        os.remove(tmp_pangenome)
        # runners.run_cgview(args, cgview_tab=cgview_data, cgview_dir=cgview_dir, images_dict=images_dict)


if __name__ == "__main__":
    args = get_args()
    main(args)
