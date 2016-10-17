import os
import glob
import shutil
import sys
from setuptools import setup, find_packages


required_progs = ['nucmer', 'show-coords', 'show-snps', 'delta-filter']
found_all_progs = True
print('Checking MUMmer programs found in path:')

for program in required_progs:
    if shutil.which(program) is None:
        found_all_progs = False
        found = '   NOT FOUND'
    else:
        found = '          OK'

    print(found, program, sep='\t')


if not found_all_progs:
    print('Cannot install because some programs from the MUMer package not found.', file=sys.stderr)
    sys.exit(1)


setup(
    name='pymummer',
    version='0.9.0',
    description='Wrapper for MUMmer',
    packages = find_packages(),
    author='Martin Hunt, Nishadi De Silva',
    author_email='path-help@sanger.ac.uk',
    url='https://github.com/sanger-pathogens/pymummer',
    test_suite='nose.collector',
    install_requires=['pyfastaq >= 3.10.0'],
    tests_require=['nose >= 1.3'],
    license='GPLv3',
    classifiers=[
        'Development Status :: 4 - Beta',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'Programming Language :: Python :: 3 :: Only',
        'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
    ],
)
