# Copyright (c) 2014-2016 Genome Research Ltd.
#
# This file is part of IVA.
#
# IVA is free software: you can redistribute it and/or modify it under
# the terms of the GNU General Public License as published by the Free Software
# Foundation; either version 3 of the License, or (at your option) any later
# version.
# This program is distributed in the hope that it will be useful, but WITHOUT
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
# FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
# details.
# You should have received a copy of the GNU General Public License along with
# this program. If not, see <http://www.gnu.org/licenses/>.

import os
import glob
import sys
import shutil
from setuptools import setup, find_packages


required_progs = [
    'kmc',
    'kmc_dump',
    'nucmer',
    'delta-filter',
    'show-coords',
    'samtools',
    'smalt',
]

found_all_progs = True
print('Checking dependencies found in path:')
for program in required_progs:
    if shutil.which(program) is None:
        found_all_progs = False
        found = ' NOT FOUND'
    else:
        found = ' OK'
    print(found, program, sep='\t')

if not found_all_progs:
    print('Cannot install because some dependencies not found.', file=sys.stderr)
    sys.exit(1)


setup(
    name='iva',
    version='1.0.6',
    description='Iterative Virus Assembler',
    packages = find_packages(),
    package_data={'iva': ['gage/*', 'ratt/*', 'read_trim/*', 'test_run_data/*']},
    author='Martin Hunt',
    author_email='path-help@sanger.ac.uk',
    url='https://github.com/sanger-pathogens/iva',
    scripts=glob.glob('scripts/*'),
    test_suite='nose.collector',
    tests_require=['nose >= 1.3'],
    install_requires=[
        'pyfastaq >= 3.10.0',
        'networkx >= 1.7',
        'pysam >= 0.8.1, <= 0.8.3',
    ],
    license='GPLv3',
    classifiers=[
        'Development Status :: 4 - Beta',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'Programming Language :: Python :: 3 :: Only',
        'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
    ],
)
