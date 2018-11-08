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
        description="fetch or update the docker images used as part of the happie pipeline",
        add_help=False)
    parser.add_argument("--virtualization", action="store",
                        help="Whether this will be run with Docker or Singularity",
                        choices=["docker", "singularity"],
                        default="docker")
    # parser.add_argument("-i", "--image_location", action="store",
    #                     help="where to store your container images")
    optional = parser.add_argument_group('optional arguments')
    optional.add_argument("-h", "--help",
                          action="help", default=argparse.SUPPRESS,
                          help="Displays this help message")
    optional.add_argument('--version',
                          action='version',
                          version='%(prog)s ' + __version__)
    args = parser.parse_args()
    return args


def install_image(args, image_name):
    if args.virtualization == "docker":
        cmd = "docker pull {image_name}".format(**locals())
    else:
        cmd = "singularity pull docker://{image_name}".format(**locals())
    print(cmd)
    subprocess.run([cmd],
                   shell=sys.platform != "win32",
                   stdout=subprocess.PIPE,
                   stderr=subprocess.PIPE,
                   check=True)

def install_prokka(args):
    # install_image(args, "ummidock/prokka")
    install_image(args, "blaxterlab/prokka")

def install_prophet(args):
    install_image(args, "nickp60/prophet")

def install_mlplasmids(args):
    install_image(args, "nickp60/mlplasmids")

def install_dimob(args):
    install_image(args, "brinkmanlab/islandpath")



def main(args=None):
    if args is None:
        args = get_args()
    # output_root = os.path.abspath(os.path.expanduser(args.output))
    install_prokka(args)
    install_prophet(args)
    install_mlplasmids(args)
    install_dimob(args)


if __name__ == "__main__":
    args = get_args()
    main(args)
