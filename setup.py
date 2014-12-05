import os
import glob
from setuptools import setup, find_packages

def read(fname):
    return open(os.path.join(os.path.dirname(__file__), fname)).read()

setup(
    name='iva',
    version='0.10.1',
    description='Iterative Virus Assembler',
    long_description=read('README.md'),
    packages = find_packages(),
    package_data={'iva': ['gage/*', 'ratt/*', 'read_trim/*']},
    author='Martin Hunt',
    author_email='mh12@sanger.ac.uk',
    url='https://github.com/sanger-pathogens/iva',
    scripts=glob.glob('scripts/*'),
    test_suite='nose.collector',
    install_requires=['nose >= 1.3', 'fastaq >= 1.6.0', 'networkx'],
    dependency_links=['http://github.com/sanger-pathogens/fastaq/tarball/master#egg=fastaq-1.7.0'],
    license='GPLv3',
)
