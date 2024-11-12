# pymummer

Python3 wrapper for running MUMmer and parsing the output. 

[![Build Status](https://travis-ci.org/sanger-pathogens/pymummer.svg?branch=master)](https://travis-ci.org/sanger-pathogens/pymummer)   
[![License: GPL v3](https://img.shields.io/badge/License-GPL%20v3-brightgreen.svg)](https://github.com/sanger-pathogens/pymummer/blob/master/LICENSE)    
[![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat-square)](http://bioconda.github.io/recipes/pymummer/README.html)   
[![Container ready](https://img.shields.io/badge/container-ready-brightgreen.svg)](https://quay.io/repository/biocontainers/pymummer)  

## Contents
  * [Introduction](#introduction)
  * [Installation](#installation)
    * [Required dependencies](#required-dependencies)
    * [Homebrew/LinuxBrew](#homebrewlinuxbrew)
    * [Using pip](#using-pip)
    * [Running the tests](#running-the-tests)
  * [Usage (for developers)](#usage-for-developers)
    * [pymummer nucmer class](#pymummer-nucmer-class)
    * [pymummer coords\_file class](#pymummer-coords_file-class)
    * [pymummer alignment class](#pymummer-alignment-class)
  * [License](#license)
  * [Feedback/Issues](#feedbackissues)

## Introduction
Runs MUMmer and parses the output.

## Installation
pymummer has the following dependencies:

### Required dependencies
* [MUMmer](http://mummer.sourceforge.net/manual/#installation)

There are a number of ways to install pymummer and details are provided below. If you encounter an issue when installing pymummer please contact your local system administrator. If you encounter a bug please log it [here](https://github.com/sanger-pathogens/pymummer/issues) or email us at path-help@sanger.ac.uk.

### Homebrew/LinuxBrew
```
brew tap homebrew/python
brew install pymummer
```

### Conda

We have provided conda environment recipes in this repo that can be used to create a fresh environment with the required dependencies. After creating a new env you can `pip install` pymummer from pypi or this repo using the commands in the next section.

```bash
# Create pymummer env 
conda env create -f environment.yml
# Activate env
conda activate pymummer
# Install pymummer
pip install pymummer
```

If you are using an M-Series Mac (ARM64 processor) you will need to create a mock Intel environment to install Mummer4 from Bioconda.

```bash
# Apple ARM64 Macs only
# Create mock Intel env
conda env create -f env_osx64.yml
# Activate env
conda activate pymummer-osx64
# Install pymummer
pip install pymummer
```

### Pip install

Install from PyPi

```bash
pip3 install pymummer
```

Or pip install the latest development version directly from this repo.

```bash
pip3 install git+https://github.com/sanger-pathogens/pymummer.git
```

### Running the tests
The test can be run from the top level directory: 

```
pytest tests
```

## Usage (for developers)

Example showing how pymummer can be used to run nucmer on a fasta file and parse the output file to produce a set of alignment objects:
```
from pymummer import coords_file, alignment, nucmer
...
runner = nucmer.Runner(reference_file, query_file, results_file) 
runner.run()
file_reader = coords_file.reader(results_file)
alignments = [coord for coord in file_reader if not coord.is_self_hit()] #Remove self hits
...
```
### pymummer nucmer class

Wraps the `nucmer`, `delta-filter`, `show-coords` and `show-snps` commands. 

Arguments:
```
__ref__			reference file  
__query__			query file  
__outfile__		output file  
__min\_id__		min\_id for delta-filter command (default None)  
__min\_length__	min\_length for delta-filter command (default None)  
__breaklen__		breaklen for nucmer command (nucmer's default is 200)   
__coords\_header__	print header in show-coords output (default True)  
__maxmatch__		maxmatch for nucmer (default False)  
__show\_snps__		run show-snps (default False)  
__snps\_header__ 	print header in show-snps output (default True)  
```

### pymummer coords_file class

Parses the nucmer output and populate an alignment object for each hit in the output
  
### pymummer alignment class

Check attributes of a hit, swap the reference and query, check if it's a self hit and so on

## License
pymummer is free software, licensed under [GPLv3](https://github.com/sanger-pathogens/pymummer/blob/master/LICENSE).

## Feedback/Issues
Please report any issues to the [issues page](https://github.com/sanger-pathogens/pymummer/issues) or email path-help@sanger.ac.uk.