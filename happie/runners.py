#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import sys
import shutil
import os
import urllib.request
import subprocess
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord


def make_containerized_cmd(args, image, dcommand, scommand,
                           indir=None, outdir=None, sing=None):
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


def run_annotation(args, contigs, prokka_dir, images_dict, skip_rename=True,
                   new_name="new_fasta.fasta", subset="wgs", log_dir=None):
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
    prokka_cmd = make_containerized_cmd(
        args=args,
        image=images_dict['prokka']["image"],
        sing=images_dict['prokka']["sing"],
        dcommand=str(
            "/input/{infile} -o /output/{outdir} --prefix {name} " +
            "--fast --cpus {cpus} 2> {log}").format(
                infile=os.path.relpath(contigs),
                outdir=os.path.relpath(prokka_dir),
                name=subset + "_" + args.name,
                cpus=args.cores,
                log=os.path.join(log_dir, subset + "_prokka.log")),
        scommand=str(
            "{infile} -o {outdir} --prefix {name} " +
            "--fast --cpus {cpus} 2> {log}").format(
                infile=contigs,
                outdir=prokka_dir,
                name=subset + "_" + args.name,
                cpus=args.cores,
                log=os.path.join(log_dir, subset + "_prokka.log"))

    )
    print(prokka_cmd)
    subprocess.run([prokka_cmd],
                   shell=sys.platform != "win32",
                   stdout=subprocess.PIPE,
                   stderr=subprocess.PIPE,
                   check=True)
    with open(os.path.join(args.output, "citing.txt"), "a") as outf:
        outf.write(images_dict['prokka']["url"] + "\n")


def run_annofilt(args, annofilt_dir, prokka_dir, images_dict,
                 subset="wgs", log_dir=None):
    if os.path.exists(annofilt_dir):
        print("removing old annofilt dir")
        shutil.rmtree(annofilt_dir)
    print("Running annofilt")
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
                log=os.path.join(log_dir, subset + "_annofilt.log")),
        scommand=str(
            "{ref} {infile} -o {outdir} " +
            "--threads {cpus} 2> {log}").format(
                ref=args.annofilt_reference,
                infile=prokka_dir,
                outdir=annofilt_dir,
                cpus=args.cores,
                log=os.path.join(log_dir, subset + "_annofilt.log"))

    )
    print(annofilt_cmd)
    subprocess.run([annofilt_cmd],
                   shell=sys.platform != "win32",
                   stdout=subprocess.PIPE,
                   stderr=subprocess.PIPE,
                   check=True)
    with open(os.path.join(args.output, "citing.txt"), "a") as outf:
        outf.write(images_dict['annofilt']["url"] + "\n")


def run_prophet(args, prokka, prophet_dir, images_dict,
                subset="wgs", log_dir=None):
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
                log=os.path.join(log_dir, subset + "_PropheET.log")),
        scommand=str(
            "--fasta_in {infilefasta} --gff_in {infilegff} " +
            "--outdir {outdir}/ --cores {cores} 2> {log}").format(
                infilefasta=os.path.relpath(prokka.fna),
                infilegff=os.path.relpath(prokka.gff),
                outdir=os.path.relpath(prophet_dir),
                cores=args.cores,
                log=os.path.join(log_dir, subset + "_PropheET.log")),
    )
    print(prophet_cmd)
    subprocess.run([prophet_cmd],
                   shell=sys.platform != "win32",
                   stdout=subprocess.PIPE,
                   stderr=subprocess.PIPE,
                   check=True)
    with open(os.path.join(args.output, "citing.txt"), "a") as outf:
        outf.write(images_dict['prophet']["url"] + "\n")


def run_mlplasmids(args, prokka, mlplasmids_results, images_dict,
                   subset="wgs", log_dir=None):
    if os.path.exists(mlplasmids_results):
        os.remove(mlplasmids_results)
    mlplasmids_cmd = make_containerized_cmd(
        args=args,
        image=images_dict['mlplasmids']['image'],
        sing=images_dict['mlplasmids']["sing"],
        dcommand=str(
            "/input/{infilefasta} /output/{outdir}  " +
            ".8 '{db}' 2> {log}").format(
                db=args.mlplasmidsdb,
                infilefasta=os.path.relpath(prokka.fna),
                outdir=os.path.relpath(mlplasmids_results),
                log=os.path.join(log_dir, subset + "_mlplasmids.log")),
        scommand=str(
            "{infilefasta} {outdir}  " +
            ".8 '{db}' 2> {log}").format(
                db=args.mlplasmidsdb,
                infilefasta=prokka.fna,
                outdir=mlplasmids_results,
                log=os.path.join(log_dir, subset + "_mlplasmids.log")),
    )
    print(mlplasmids_cmd)
    subprocess.run([mlplasmids_cmd],
                   shell=sys.platform != "win32",
                   stdout=subprocess.PIPE,
                   stderr=subprocess.PIPE,
                   check=True)
    with open(os.path.join(args.output, "citing.txt"), "a") as outf:
        outf.write(images_dict['mlplasmids']["url"] + "\n")


def run_plasflow(args, prokka, plasflow_results, images_dict,
                 subset="wgs", log_dir=None):
    if os.path.exists(plasflow_results):
        os.remove(plasflow_results)
    plasflow_cmd = make_containerized_cmd(
        args=args,
        image=images_dict['plasflow']['image'],
        sing=images_dict['plasflow']["sing"],
        dcommand=str(
            "--input /input/{infilefasta} --output /output/{outdir}  " +
            " 2> {log}").format(
                infilefasta=os.path.relpath(prokka.fna),
                outdir=os.path.relpath(plasflow_results),
                log=os.path.join(log_dir, subset + "_plasflow.log")),
        scommand=str(
            "--input {infilefasta} --output {outdir}  " +
            " 2> {log}").format(
                infilefasta=prokka.fna,
                outdir=plasflow_results,
                log=os.path.join(log_dir, subset + "_plasflow.log")),
    )
    print(plasflow_cmd)
    subprocess.run([plasflow_cmd],
                   shell=sys.platform != "win32",
                   stdout=subprocess.PIPE,
                   stderr=subprocess.PIPE,
                   check=True)
    with open(os.path.join(args.output, "citing.txt"), "a") as outf:
        outf.write(images_dict['plasflow']["url"] + "\n")


def run_mobsuite(args, prokka,  mobsuite_results, images_dict,
                 subset="wgs", log_dir=None):
    if os.path.exists(mobsuite_results):
        shutil.rmtree(mobsuite_results)
    os.makedirs(mobsuite_results)
    mobsuite_cmd = make_containerized_cmd(
        args=args,
        image=images_dict['mobsuite']['image'],
        sing=images_dict['mobsuite']["sing"],
        dcommand=str(
            "--infile /input/{infilefasta} --outdir /output/{outdir}  " +
            "--run_typer 2> {log}").format(
                infilefasta=os.path.relpath(prokka.fna),
                outdir=os.path.relpath(mobsuite_results),
                log=os.path.join(log_dir, subset + "_mobsuite.log")),
        scommand=str(
            "--infile {infilefasta} --outdir {outdir}  " +
            "--run_typer 2> {log}").format(
                infilefasta=prokka.fna,
                outdir=mobsuite_results,
                log=os.path.join(log_dir, subset + "_mobsuite.log")),
    )
    print(mobsuite_cmd)
    subprocess.run([mobsuite_cmd],
                   shell=sys.platform != "win32",
                   stdout=subprocess.PIPE,
                   stderr=subprocess.PIPE,
                   check=True)
    with open(os.path.join(args.output, "citing.txt"), "a") as outf:
        outf.write(images_dict['mobsuite']["url"] + "\n")


def run_dimob(args, prokka, island_results, images_dict,
              subset="wgs", log_dir=None):
    island_dir = os.path.dirname(island_results)
    island_seqs_dir = os.path.join(island_dir, "tmp")
    island_results_dir = os.path.join(island_dir, "results")
    if os.path.exists(island_dir):
        shutil.rmtree(island_dir)
    for path in [island_dir, island_seqs_dir, island_results_dir]:
        os.makedirs(path)
    # dimob doesnt process assemblies;
    # we write out each sequence to a tmp file prior to processing
    # {id: {gbk: path_to_gbk, results: path_to_results}}
    island_path_results = {}
    with open(prokka.gbk, "r") as inf:
        for i, rec in enumerate(SeqIO.parse(inf, "genbank")):
            gbk = os.path.join(island_seqs_dir, "%i.gbk" % i)
            result = os.path.join(island_results_dir, "results_%s" % i)
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

    for k, v in island_path_results.items():
        # Note that dimob does not have an exe
        island_cmd = make_containerized_cmd(
            args=args,
            sing=images_dict['dimob']["sing"],
            image=images_dict['dimob']['image'],
            dcommand=str(
                "/input/{infile} /output/{outdir}").format(
                    infile=os.path.relpath(v["gbk"]),
                    outdir=os.path.relpath(v['result'])),
            scommand=str(
                "{infile} {outdir}").format(
                    infile=v["gbk"],
                    outdir=v['result']),
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
    with open(os.path.join(args.output, "citing.txt"), "a") as outf:
        outf.write(images_dict['dimob']["url"] + "\n")


def run_abricate(args, abricate_dir, mobile_fasta, images_dict,
                 all_results, subset="wgs", log_dir=None):
    # remove old results
    if os.path.exists(abricate_dir):
        shutil.rmtree(abricate_dir)
    os.makedirs(abricate_dir)
    cmds = []
    outfiles = []
    for db in args.analyses:
        if db == "antismash":
            continue
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
                    log=os.path.join(log_dir, subset + "_abricate")),
            scommand=str(
                "--db {db}  {infilefasta} > {outdir}/{db}.tab  " +
                "2> {log}_{db}_log.txt").format(
                    db=db,
                    infilefasta=mobile_fasta,
                    outdir=abricate_dir,
                    log=os.path.join(log_dir, subset + "_abricate")),
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
    with open(os.path.join(args.output, "citing.txt"), "a") as outf:
        outf.write(images_dict['abricate']["url"] + "\n")


def make_cgview_tab_file(args, seqlen, cgview_entries):
    """
    http://wishart.biology.ualberta.ca/cgview/tab_input.html
    """
    results = []
    header = [
        "#{}".format(args.name),
        "%{}".format(seqlen + 1),
        "!strand\tslot\tstart\tstop\ttype\tlabel\tmouseover\thyperlink"
    ]
    results.extend(header)
    results.extend(cgview_entries)
    return results


def run_cgview(args, cgview_tab, cgview_dir, images_dict,
               subset="wgs", log_dir=None):
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
                log=os.path.join(log_dir, subset + "_cgview.log")),
        scommand=str(
            "-i {infilefasta} -f svg -o {outdir} -I T").format(
                exe=images_dict['cgview']['exe'],
                infilefasta=cgview_tab,
                outdir=os.path.join(cgview_dir, "cgview.svg"),
                log=os.path.join(log_dir, subset + "_cgview.log")),

    )
    print(cgview_cmd)
    subprocess.run([cgview_cmd],
                   shell=sys.platform != "win32",
                   stdout=subprocess.PIPE,
                   stderr=subprocess.PIPE,
                   check=True)
    with open(os.path.join(args.output, "citing.txt"), "a") as outf:
        outf.write(images_dict['cgview']["url"] + "\n")


def run_antismash(args, antismash_dir, gbk, images_dict,
                 subset="wgs", log_dir=None):
    # remove old results
    if os.path.exists(antismash_dir):
        shutil.rmtree(antismash_dir)
    # check if the rrunner script is in the PATH
    this_cmd = make_containerized_cmd(
        args=args,
        image=images_dict['antismash']['image'],
        sing=images_dict['antismash']["sing"],
        dcommand=str(
            "--cpus {cores} /input/{gbk} --output-dir /output/{outdir} " +
            "--genefinding-tool prodigal " +
            "2> {log}_log.txt").format(
                cores=args.cores,
                gbk=os.path.relpath(gbk),
                outdir=os.path.relpath(antismash_dir),
                log=os.path.join(log_dir, subset + "_antismash")),
        scommand=str(
            "--cpus {cores} {gbk} --output-dir {outdir} " +
            "--genefinding-tool prodigal " +
            "2> {log}_log.txt").format(
                cores=args.cores,
                gbk=os.path.relpath(gbk),
                outdir=antismash_dir,
                log=os.path.join(log_dir, subset + "_antismash")),
    )
    print(this_cmd)
    subprocess.run([this_cmd],
                   shell=sys.platform != "win32",
                   stdout=subprocess.PIPE,
                   stderr=subprocess.PIPE,
                   check=True)
