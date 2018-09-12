"""
Setup for mobilephone
"""

# Always prefer setuptools over distutils
from setuptools import setup, find_packages
# To use a consistent encoding
import re
from codecs import open
from os import path
import sys

try: # for pip >= 10
    from pip._internal.req import parse_requirements
except ImportError: # for pip <= 9.0.3
    from pip.req import parse_requirements

here = path.abspath(path.dirname(__file__))

VERSIONFILE = "mp/_version.py"
verstrline = open(VERSIONFILE, "rt").read()
VSRE = r"^__version__ = ['\"]([^'\"]*)['\"]"
mo = re.search(VSRE, verstrline, re.M)
if mo:
    verstr = mo.group(1)
else:
    raise RuntimeError("Unable to find version string in %s." % (VERSIONFILE,))


if sys.version_info <= (3, 0):
    sys.stderr.write("ERROR: mobilephone requires Python 3.5 " +
                     "or above...exiting.\n")
    sys.exit(1)

## parse requirements file
install_reqs = parse_requirements("requirements.txt",
                                  session=False)
requirements = [str(ir.req) for ir in install_reqs]

setup(
    name='mobilephone',
    version=verstr,

    description='mobilephone: just the mobile elements',
    # long_description=long_description,
    long_description="""
    check out the GitHub
    repo for the real README.md file
    """,

    url='https://github.com/nickp60/mobilephone',

    # Author details
    author='Nick Waters',
    author_email='nickp60@gmail.com',
    license='MIT',
    # handle requirments
    install_requires=requirements,
    # See https://pypi.python.org/pypi?%3Aaction=list_classifiers
    classifiers=[
        'Development Status :: 3 - Alpha',
        'Intended Audience :: Developers',
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python :: 3.5',
        'Programming Language :: Python :: 3.6',
    ],
    keywords='bioinformatics evolution genomics development',
    packages=find_packages(),
    # include_package_data=True,
    # package_data={
    #    '': [path.join(__name__, "tests", "refs/*")],
    # },
    entry_points={
       'console_scripts': [
           'mobilephone=mp.mobilephone:main',
       ],
    },
)
