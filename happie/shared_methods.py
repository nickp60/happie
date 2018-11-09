#!/usr/bin/env python3
#-*- coding: utf-8 -*-
import sys
import os
import yaml
import pkg_resources


def get_containers_manifest():
    resource_package = pkg_resources.Requirement.parse("happie")
    print(resource_package)
    print(pkg_resources.resource_listdir("happie", "data"))
    datapath = pkg_resources.resource_filename(resource_package, 'happie/data/containers.yaml')
    print(datapath)
    with open(datapath, "r") as inf:
        return yaml.load(inf)


def parse_docker_images(config):
    images_dict = {}
    try:
        programs = config['programs']
    except KeyError:
        print("Error with config file!  no 'program' header found")
        sys.exit(1)
    for k, v in programs.items():
        image_name = "{0}/{1}:{2}".format(
            v['maintainer'], v['name'], v['version'])
        images_dict[k] = {"image": image_name,
                          "exe": v["exe"]}
    return images_dict
