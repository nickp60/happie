# -*- coding: utf-8 -*-
"""
Created on Tue Aug 30 08:57:31 2016
@author: nicholas

"""
import sys
import logging
import shutil
# import subprocess
import os
import unittest
# import multiprocessing

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from argparse import Namespace
from unittest.mock import MagicMock, patch

# I hate this line but it works :(
sys.path.append(os.path.join(
    os.path.dirname(os.path.dirname(__file__)), "happie"))


from happie import happie as hh
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
            os.path.dirname(__file__), "references")
        self.ref_contigs = os.path.join(
            self.ref_dir,
            'GCA_000225105.2_ASM22510v2_genomic_fragments.fna')
        self.ref_prophet = os.path.join(
            self.ref_dir,
            'phages_coords')
        self.ref_mlplasmids = os.path.join(
            self.ref_dir,
            'results.txt')
        self.maxDiff = 2000
        self.to_be_removed = []
        if not os.path.exists(self.test_dir):
            os.makedirs(self.test_dir, exist_ok=True)

    def test_parse_prophet_results(self):
        ref = [
            ['prophet', 'prophage', 'AFDU01000019.1', 177, 24461],
            ['prophet', 'prophage', 'AFDU01000019.1', 808102, 848966],
            ['prophet', 'prophage', 'AFDU01000012.1', 216544, 249812],
            ['prophet', 'prophage', 'AFDU01000012.1', 352274, 368386],
            ['prophet', 'prophage', 'AFDU01000014.1', 159, 8432],
            ['prophet', 'prophage', 'AFDU01000015.1', 385, 39315],
            ['prophet', 'prophage', 'AFDU01000031.1', 57319, 73419],
            ['prophet', 'prophage', 'AFDU01000005.1', 317383, 360232],
            ['prophet', 'prophage', 'AFDU01000011.1', 2908, 15720]
        ]
        self.assertEqual(
            hh.parse_prophet_results(self.ref_prophet),
            ref
        )

    def test_condensce_regions(self):
        all_results = [
            ['prophet', 'prophage', 'A', '5', '10'],
            ['prophet', 'prophage', 'A', '1', '15'],
            ['prophet', 'prophage', 'B', '55', '66'],
            ['prophet', 'prophage', 'B', '34', '100'],
            ['mlplasmids', 'plasmid', '"A"', '0', '1214'],
        ]
        non_overlapping_results = hh.condensce_regions(all_results)
        print(non_overlapping_results)

    def test_parse_mlplasmids_results(self):
        ref = [
            ['mlplasmids', 'plasmid', '"AFDU01000034.1"', '0', '1214'],
        ]

        self.assertEqual(
            hh.parse_mlplasmids_results(self.ref_mlplasmids),
            ref
        )

    def test_write_annotated_mobile_genome(self):
        outf = os.path.join(self.test_dir, "mobile_genome.fasta")
        all_results = [
            ['prophet', 'prophage', 'AFDU01000019.1', '5', '10'],
            ['prophet', 'prophage', 'AFDU01000005.1', '1', '15'],
            ['prophet', 'prophage', 'AFDU01000012.1', '55', '66'],
            ['prophet', 'prophage', 'AFDU01000031.1', '34', '100'],
        ]

        seq_length, cgview_entries = hh.write_annotated_mobile_genome(
            contigs=prokka_fna,
            output_path=reference_mobile_genome_path,
            non_overlapping_results=non_overlapping_results,
            all_results=all_results)

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
