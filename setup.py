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
    version='0.11.1',
    description='Iterative Virus Assembler',
    packages = find_packages(),
    package_data={'iva': ['gage/*', 'ratt/*', 'read_trim/*']},
    author='Martin Hunt',
    author_email='mh12@sanger.ac.uk',
    url='https://github.com/sanger-pathogens/iva',
    scripts=glob.glob('scripts/*'),
    test_suite='nose.collector',
    install_requires=['nose >= 1.3', 'pyfastaq >= 3.0.1', 'networkx', 'pysam'],
    license='GPLv3',
)
