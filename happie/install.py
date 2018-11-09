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
import yaml
import pkg_resources

from Bio.Seq import Seq
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

from . import __version__
from . import shared_methods as sm

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


def test_exe_exists(args):
    if args.virtualization == "docker":
        cmd = "docker --version"
    else:
        cmd = "singularity --version"
    subprocess.run([cmd],
                   shell=sys.platform != "win32",
                   stdout=subprocess.PIPE,
                   stderr=subprocess.PIPE,
                   check=True)

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


def install_programs(args, config):
    images_dict = sm.parse_docker_images(config)
    for k, v in images_dict.items():
        print("installing %s", k)
        install_image(args, v['image'])


def main(args=None):
    if args is None:
        args = get_args()
    # output_root = os.path.abspath(os.path.expanduser(args.output))
    try:
        test_exe_exists(args)
    except Exception as e:
        print(e)
        print("Error: %s executable is not install or not in PATH",
              args.virtualization)
        sys.exit(1)
    config = sm.get_containers_manifest()
    install_programs(args, config)


if __name__ == "__main__":
    args = get_args()
    main(args)
