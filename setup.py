#!/usr/bin/python
"""
Duplex sequencing method processing tools

@author: Jose Fernandez Navarro <jc.fernandez.navarro@gmail.com>
"""

import os
import io
import glob
from setuptools import setup, find_packages

# Get the long description from the relevant file
here = os.path.abspath(os.path.dirname(__file__))
with io.open(os.path.join(here, 'README.md'), encoding='utf-8') as f:
    long_description = f.read()

setup(
  name='duplexseq',
  version="0.1.0",
  description=__doc__.split("\n", 1)[0],
  long_description=long_description,
  keywords='duplex sequencing processing pipeline',
  author='Jose Fernandez Navarro',
  author_email='jc.fernandez.navarro@gmail.com',
  license='MIT',
  packages=find_packages(),
  include_package_data=False,
  package_data={'': ['RELEASE-VERSION']},
  zip_safe=False,
  install_requires=[
    'setuptools',
    'pysam',
    'biopython'
  ],
  #test_suite = 'tests',
  scripts=glob.glob('*.py'),
  classifiers=[
    'Development Status :: 4 - Beta',
    'Intended Audience :: Science/Research',
    'Topic :: Software Development',
    'Topic :: Scientific/Engineering :: Bio-Informatics',
    'License :: MIT:: Copyright Jose Fernandez Navarro',
    'Programming Language :: Python :: 3.6',
    'Programming Language :: Python :: 3.7',
    'Environment :: Console',
  ],
)
