# -*- coding: utf-8 -*-
"""
Created on Tue Aug 30 08:57:31 2016
@author: nicholas

"""
import sys
import logging
import shutil
# import subprocess
import glob
import os
import unittest
# import multiprocessing

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from argparse import Namespace
from unittest.mock import MagicMock, patch

from . import happie as hh
from . import shared_methods as sm
from . import runners as runners
from . import parsers as parsers

sys.dont_write_bytecode = True

logger = logging


@unittest.skipIf((sys.version_info[0] != 3) or (sys.version_info[1] < 5),
                 "Subprocess.call among other things wont run if tried " +
                 " with less than python 3.5")
class MpTestCase(unittest.TestCase):
    """ tests for happie
    """
    def setUp(self):
        self.test_dir = os.path.join(os.path.dirname(__file__),
                                     "output_mp_tests")
        self.ref_dir = os.path.join(
            os.path.dirname(__file__), "test_files")
        self.ref_contigs = os.path.join(
            self.ref_dir,
            'GCA_000225105.2_ASM22510v2_genomic_fragments.fna')
        self.ref_prophet = os.path.join(
            self.ref_dir,
            'phages_coords')
        self.ref_antismash = os.path.join(
            self.ref_dir, "lcl_1.region001.gbk")
        self.ref_mlplasmids = os.path.join(
            self.ref_dir,
            'results.txt')
        self.spades_assembly = os.path.join(
            self.ref_dir,
            'contigs.fasta')
        self.alt_assembly = os.path.join(
            self.ref_dir,
            'alt.fa')
        self.maxDiff = 2000
        self.to_be_removed = []
        if not os.path.exists(self.test_dir):
            os.makedirs(self.test_dir, exist_ok=True)

    def test_parse_prophet_results(self):
        ref = [
            ['prophet', 'prophages', 'AFDU01000019.1', 177, 24461],
            ['prophet', 'prophages', 'AFDU01000019.1', 808102, 848966],
            ['prophet', 'prophages', 'AFDU01000012.1', 216544, 249812],
            ['prophet', 'prophages', 'AFDU01000012.1', 352274, 368386],
            ['prophet', 'prophages', 'AFDU01000014.1', 159, 8432],
            ['prophet', 'prophages', 'AFDU01000015.1', 385, 39315],
            ['prophet', 'prophages', 'AFDU01000031.1', 57319, 73419],
            ['prophet', 'prophages', 'AFDU01000005.1', 317383, 360232],
            ['prophet', 'prophages', 'AFDU01000011.1', 2908, 15720]
        ]
        self.assertEqual(
            parsers.parse_prophet_results(self.ref_prophet),
            ref
        )

    def test_parse_antismash_results(self):
        test_results = parsers.parse_antismash_results(
            self.ref_dir, region="wgs")
        ref_results = [["lcl_1", "wgs", "antismash",
                        ["bacteriocin"], 525047, 525299]]
        print(test_results)
        print(ref_results)
        assert ref_results == test_results



    def test_condensce_regions(self):
        all_results = [
            ['prophet', 'prophages', 'A', '5', '10'],
            ['prophet', 'prophages', 'A', '1', '15'],
            ['prophet', 'prophages', 'B', '55', '66'],
            ['prophet', 'prophages', 'B', '34', '100'],
            ['mlplasmids', 'plasmids', '"A"', '0', '1214'],
        ]
        non_overlapping_results = hh.condensce_regions(all_results)
        print(non_overlapping_results)

    def test_parse_mlplasmids_results(self):
        ref = [
            ['mlplasmids', 'plasmids', 'AFDU01000034.1', 1, 1214],
        ]

        self.assertEqual(
            parsers.parse_mlplasmids_results(self.ref_mlplasmids),
            ref
        )

    def test_write_annotated_mobile_genome(self):
        outf = os.path.join(self.test_dir, "mobile_genome.fasta")
        all_results = [
            ['prophet', 'prophages', 'AFDU01000019.1', 5, 10],
            ['prophet', 'prophages', 'AFDU01000005.1', 1, 15],
            ['prophet', 'prophages', 'AFDU01000012.1', 55, 66],
            ['prophet', 'prophages', 'AFDU01000031.1', 34, 100],
        ]

        seq_length, cgview_entries = hh.write_annotated_mobile_genome(
            contigs=self.spades_assembly,
            output_path=self.test_dir,
            non_overlapping_results=all_results,
            seed=123,
            all_results=all_results)

    def test_QC_bug(self):
        args = Namespace(contigs=self.spades_assembly, output="./")
        if os.path.exists("sublogs"):
            shutil.rmtree("sublogs")
        os.makedirs("sublogs")
        hh.QC_bug(
            args,
            QC_dir=self.test_dir,
            min_length=10,
            max_length=100000000,
            cov_threshold=.2)

    def test_alt_cov_filter(self):
        args = Namespace(contigs=self.alt_assembly, output="./")
        if os.path.exists("sublogs"):
            shutil.rmtree("sublogs")
        os.makedirs("sublogs")
        hh.QC_bug(
            args,
            QC_dir=self.test_dir,
            min_length=10,
            min_contig_length=10,
            max_length=100000000,
            cov_threshold=.2)

    # def test_remove_bad_contig(self):
    #     outfile = os.path.join(self.test_dir, "contigs_minus_NODE_2.fasta")
    #     remove_bad_contig(infile=self.bad_fasta,
    #                       # outfile=self.good_fasta,
    #                       outfile=outfile,
    #                       # bad_name="NODE_3_length_222_cov_5.70526",
    #                       bad_name="NODE_3_",
    #                       logger=logger)
    #     with open(outfile, 'r') as of:
    #         recs = list(SeqIO.parse(of, 'fasta'))
    #         for i in recs:
    #             self.assertTrue("NODE_3_" not in i.id)
    #     self.to_be_removed.append(outfile)

    # def test_fail_remove_bad_contig(self):
    #     outfile = os.path.join(self.test_dir, "contigs_minus_NODE_3.fasta")
    #     with self.assertRaises(ValueError):
    #         remove_bad_contig(infile=self.bad_fasta,
    #                           outfile=outfile,
    #                           bad_name="NODE_",
    #                           logger=logger)

    # def test_append_replacement_contigs(self):
    #     outfile = os.path.join(self.test_dir, "contigs_minus_NODE_3.fasta")
    #     remove_bad_contig(infile=self.bad_fasta,
    #                       # outfile=self.good_fasta,
    #                       outfile=outfile,
    #                       # bad_name="NODE_3_length_222_cov_5.70526",
    #                       bad_name="NODE_3_",
    #                       logger=logger)
    #     append_replacement_contigs(infile=self.good_fasta, outfile=outfile,
    #                                name_list="NODE_4_:NODE_5_".split(":"),
    #                                logger=logger)
    #     self.assertEqual(md5(outfile),
    #                      md5(os.path.join(self.ref_dir,
    #                                       "ref_swapped_contigs.fasta")))
    #     self.to_be_removed.append(outfile)

    # def test_fail_append_replacement_contigs(self):
    #     outfile = os.path.join(self.test_dir, "contigs_minus_NODE_3.fasta")
    #     remove_bad_contig(infile=self.bad_fasta,
    #                       outfile=outfile,
    #                       bad_name="NODE_3_",
    #                       logger=logger)
    #     with self.assertRaises(ValueError):
    #         append_replacement_contigs(
    #             infile=self.good_fasta, outfile=outfile,
    #             name_list="NODE_4_:NODE_5_:NODE_45_".split(":"),
    #             logger=logger)
    #     self.to_be_removed.append(outfile)


    def tearDown(self):
        """ delete temp files if no errors
        """
        if len(self.to_be_removed) != 0:
            for filename in self.to_be_removed:
                if os.path.exists(filename):
                    os.unlink(filename)

if __name__ == '__main__':
    unittest.main()
