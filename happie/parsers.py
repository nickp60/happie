#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import os
import glob
from Bio import SeqIO

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
            if i in [4, 5]:
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
            if i in [4, 5]:
                subresults.append(int(line[i]))
            else:
                subresults.append(line[i])
        templated.append(subresults)
    return templated


def parse_plasflow_results(results):
    results_text = []
    templated = []
    with open(results) as inf:
        for i, line in enumerate(inf):
            if i == 0:
                pass
            else:
                results = ["plasflow", "plasmids", "1"]
                results.extend(line.strip().split("\t"))
                results_text.append(results)
    for line in results_text:
        subresults = []
        if not line[8].startswith("plasmid"):
            continue
        for i in [0, 1, 5, 2, 6]:
            if i in [2, 6]:
                subresults.append(int(line[i]))
            else:
                subresults.append(line[i])
        templated.append(subresults)
    return templated


def parse_mobsuite_results(results):
    results_text = []
    templated = []
    with open(results) as inf:
        for i, line in enumerate(inf):
            if i == 0:
                pass
            else:
                results = ["mobsuite", "plasmids", "1"]
                results.extend(line.strip().split("\t"))
                results_text.append(results)
    for line in results_text:
        subresults = []
        if line[4].startswith("chromosome"):
            continue
        for i in [0, 1, 5, 2, 6]:
            if i in [2, 6]:
                subresults.append(int(line[i]))
            elif i == 5.:
                subresults.append(line[i].split("|")[1])
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
                if i in [2, 7]:
                    subresults.append(int(line[i]))
                else:
                    subresults.append(line[i])
                # subresults.append(line[i].replace('"', ""))
            templated.append(subresults)
            # print(subresults)
    return templated



def parse_antismash_results(results_dir, region):
    results = []
    d = os.path.join(results_dir, "")
    for f in glob.glob(d + "*region*.gbk"):
        print("parsing antistmash result %s" %f)
        with open(f, "r") as inf:
            for record in SeqIO.parse(inf, "genbank"):
                for feat in record.features:
                    if feat.type == "protocluster":
                        start, stop = feat.qualifiers.get(
                            "core_location")[0].replace(
                                "[", "").split("]")[0].split(":")
                        results.append([
                            record.id,
                            region,
                            "antismash",
                            feat.qualifiers.get("product"),
                            int(start), int(stop)])
    return results
