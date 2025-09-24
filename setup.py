#!/usr/bin/env python
# Licensed under a 3-clause BSD style license - see LICENSE.rst
# Started from https://github.com/desihub/prospect/blob/master/setup.py
#
# Standard imports
#
import glob
import os
import sys
#
# setuptools' sdist command ignores MANIFEST.in
#
from distutils.command.sdist import sdist as DistutilsSdist
from setuptools import setup, find_packages
#
# DESI support code.
#
have_desiutil = False
try:
    import desiutil.setup as ds
except ImportError:
    have_desiutil = False
#
# Begin setup
#
setup_keywords = dict()
#
# THESE SETTINGS NEED TO BE CHANGED FOR EVERY PRODUCT.
#
setup_keywords['name'] = 'LSS'
setup_keywords['description'] = 'LSS catalog package'
setup_keywords['author'] = 'DESI Collaboration'
setup_keywords['author_email'] = 'desi-data@desi.lbl.gov'
setup_keywords['license'] = 'BSD'
setup_keywords['url'] = 'https://github.com/desihub/LSS'
#
# END OF SETTINGS THAT NEED TO BE CHANGED.
#
if have_desiutil:
    setup_keywords['version'] = ds.get_version(setup_keywords['name'])
else:
    try:
        with open(os.path.join('py', setup_keywords['name'], '_version.py')) as v:
            setup_keywords['version'] = v.read().split('=')[1].strip().strip("'").strip('"')
    except FileNotFoundError:
        setup_keywords['version'] = '0.0.1'
#
# Use README.rst as long_description.
#
setup_keywords['long_description'] = ''
if os.path.exists('README.rst'):
    with open('README.rst') as readme:
        setup_keywords['long_description'] = readme.read()
#
# Set other keywords for the setup function.  These are automated, & should
# be left alone unless you are an expert.
#
# Treat everything in bin/ except *.rst as a script to be installed.
#
if os.path.isdir('bin'):
    setup_keywords['scripts'] = [fname for fname in glob.glob(os.path.join('bin', '*'))
        if os.path.basename(fname).startswith('mkCat_') and not os.path.basename(fname).endswith('.rst')]

setup_keywords['provides'] = [setup_keywords['name']]
setup_keywords['python_requires'] = '>=3.5'
setup_keywords['zip_safe'] = False
setup_keywords['use_2to3'] = False
setup_keywords['packages'] = find_packages('py')
setup_keywords['package_dir'] = {'': 'py'}

setup_keywords['cmdclass'] = {'sdist': DistutilsSdist}
if have_desiutil:
    setup_keywords['cmdclass']['module_file'] = ds.DesiModule
    setup_keywords['cmdclass']['version'] = ds.DesiVersion
    setup_keywords['cmdclass']['test'] = ds.DesiTest
    setup_keywords['cmdclass']['api'] = ds.DesiAPI
#setup_keywords['test_suite'] = '{name}.test.{name}_test_suite'.format(**setup_keywords)
#
# Add internal data directories.
#
setup_keywords['package_data'] = {'LSS': ['data/*.dat']}
#
# Run setup command.
#
setup(**setup_keywords)
