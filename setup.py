import os
import glob
from setuptools import setup, find_packages

def read(fname):
    return open(os.path.join(os.path.dirname(__file__), fname)).read()

setup(
    name='pymummer',
    version='0.0.1',
    description='Wrapper for MUMmer',
    long_description=read('README.md'),
    packages = find_packages(),
    author='Martin Hunt, Nishadi De Silva',
    author_email='path-help@sanger.ac.uk',
    url='https://github.com/sanger-pathogens/pymummer',
    test_suite='nose.collector',
    install_requires=['nose >= 1.3', 'fastaq >= 2.0.0'],
    dependency_links=['http://github.com/sanger-pathogens/fastaq/tarball/master#egg=fastaq-2.0.0'],
    license='GPLv3',
)
