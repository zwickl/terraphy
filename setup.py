#! /usr/bin/env python

##############################################################################
##  Phylogenomic tools package.
##
##  Copyright 2013 Derrick J. Zwickl
##  All rights reserved.
##
##  See "LICENSE.txt" for terms and conditions of usage.
##
##  If you use this work or any portion thereof in published work,
##  please cite it as:
##
##     REF
##
##############################################################################

import sys
import os

#HACK OF DENDROPY setup.py BELOW

###############################################################################
# setuptools/distutils/etc. import and configuration

try:
    import ez_setup
    try:
        ez_setup_path = " ('" + os.path.abspath(ez_setup.__file__) + "')"
    except OSError:
        ez_setup_path = ""
    sys.stderr.write("using ez_setup%s\n" %  ez_setup_path)
    ez_setup.use_setuptools()
    import setuptools
    try:
        setuptools_path = " ('" +  os.path.abspath(setuptools.__file__) + "')"
    except OSError:
        setuptools_path = ""
    sys.stderr.write("using setuptools%s\n" % setuptools_path)
    from setuptools import setup, find_packages
except ImportError, e:
    sys.stderr.write("using distutils\n")
    from distutils.core import setup
    sys.stderr.write("using canned package list\n")
    PACKAGES = ['terraphy']
    EXTRA_KWARGS = {}
else:
    sys.stderr.write("searching for packages\n")
    PACKAGES = find_packages()
    EXTRA_KWARGS = dict(
        install_requires = ['setuptools', 'dendropy', 'python-graph-core'], 
        include_package_data=True
    )

#this would be necessary to read gffs and rewrite alignments based on them (biopython too)
#install_requires = ['setuptools', 'dendropy', 'distribute>=0.6.49', 'bcbio-gff'],

PACKAGE_DIRS = [p.replace(".", os.path.sep) for p in PACKAGES]
PACKAGE_INFO = [("% 40s : %s" % p) for p in zip(PACKAGES, PACKAGE_DIRS)]
sys.stderr.write("packages identified:\n%s\n" % ("\n".join(PACKAGE_INFO)))
ENTRY_POINTS = {}

###############################################################################
# Script paths

SCRIPT_SUBPATHS = [
    ['scripts', 'terraphy.main.py'],
    ['scripts', 'maketriplets.py'],
]
SCRIPTS = [os.path.join(*i) for i in SCRIPT_SUBPATHS]
sys.stderr.write("\nscripts identified: %s\n" % "\n\t".join(SCRIPTS))

'''
try:
    mpl_needed = '1.2'
    import matplotlib
    mpl_ver = matplotlib.__version__
    if not float(mpl_ver[:3]) >= float(mpl_needed):
        sys.stderr.write('*******************\nNOTE: matplotlib VERSION %s OR LATER IS REQUIRED (%s found).\nIT IS NECESSARY FOR ALL PLOTTING FUNCTIONALITY.\nYOU MAY NEED TO INSTALL IT MANUALLY (see http://matplotlib.org/downloads.html)\n*******************\n' % (mpl_needed, mpl_ver))

except ImportError:
    sys.stderr.write('*******************\nNOTE: matplotlib PACKAGE COULD NOT BE IMPORTED.\nIT IS NECESSARY FOR ALL PLOTTING FUNCTIONALITY.\nYOU MAY NEED TO INSTALL IT MANUALLY (see http://matplotlib.org/downloads.html)\n*******************\n')
'''
'''don't need biopython yet
try:
    bp_needed = '1.5'
    import Bio
    bp_ver = Bio.__version__
    if not float(bp_ver[:3]) >= float(bp_needed):
        sys.stderr.write('*******************\nNOTE: biopython VERSION %s OR LATER IS REQUIRED (%s found).\nIT IS NECESSARY FOR ALL PLOTTING FUNCTIONALITY.\nYOU MAY NEED TO INSTALL IT MANUALLY (see http://biopython.org/wiki/Download)\n*******************\n' % (bp_needed, bp_ver))

except ImportError:
    sys.stderr.write('*******************\nNOTE: matplotlib PACKAGE COULD NOT BE IMPORTED.\nIT IS NECESSARY FOR ALL PLOTTING FUNCTIONALITY.\nYOU MAY NEED TO INSTALL IT MANUALLY (see http://matplotlib.org/downloads.html)\n*******************\n')
'''

###############################################################################
# setuptools/distuils command extensions

try:
    from setuptools import Command
except ImportError:
    sys.stderr.write("setuptools.Command could not be imported: setuptools extensions not available\n")
else:
    sys.stderr.write("setuptools command extensions are available\n")
    command_hook = "distutils.commands"
    ENTRY_POINTS[command_hook] = []

    ###########################################################################
    # coverage
    '''
    from dendropy.test.support import coverage_analysis
    if coverage_analysis.DENDROPY_COVERAGE_ANALYSIS_AVAILABLE:
        sys.stderr.write("coverage analysis available ('python setup.py coverage')\n")
        ENTRY_POINTS[command_hook].append("coverage = dendropy.test.support.coverage_analysis:CoverageAnalysis")
    else:
        sys.stderr.write("coverage analysis not available\n")
    '''


###############################################################################
# Main setup

__version__ = 1.0
EXTRA_KWARGS["zip_safe"] = True

### compose long description ###
long_description = open('README.txt').read()

setup(name='terraphy',
      version=__version__,
      author='Derrick Zwickl',
      author_email='zwickl@email.arizona.edu',
      url='https://github.com/zwickl/terraphy',
      description='Scripts for analyzing phylogenetic terraces.',
      license='GPL3',
      packages=PACKAGES,
      package_dir=dict(zip(PACKAGES, PACKAGE_DIRS)),
      scripts = SCRIPTS,
      long_description=long_description,
      entry_points = ENTRY_POINTS,
      classifiers = [
            "Environment :: Console",
            "Intended Audience :: Developers",
            "Intended Audience :: Science/Research",
            "License :: OSI Approved :: GPL3",
            "Natural Language :: English",
            "Operating System :: OS Independent",
            "Programming Language :: Python",
            "Topic :: Scientific/Engineering :: Bio-Informatics",
            ],
      keywords='phylogenetics phylogeny phylogenies phylogenomics evolution evolutionary biology systematics phyloinformatics bioinformatics',
      **EXTRA_KWARGS
      )
