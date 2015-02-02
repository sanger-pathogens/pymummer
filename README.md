pymummer
========

Python3 wrapper for running MUMmer and parsing the output. 

Installation
------------

###Pre-requisites###

The MUMmer package must be installed.
Instructions to install MUMmer can be found [here](http://mummer.sourceforge.net/manual/#installation)
    
		
###Installation###

Install with

    pip3 install pymummer


Usage (for developers)
----------------------

Example showing how pymummer can be used to run nucmer on a fasta file and parse the output file to produce a set of alignment objects:

	from pymummer import coords_file, alignment, nucmer
	...
	runner = nucmer.Runner(reference_file, query_file, results_file) 
	runner.run()
	file_reader = coords_file.reader(results_file)
	alignments = [coord for coord in file_reader if not coord.is_self_hit()] #Remove self hits
	...

###pymummer nucmer class###

Wraps the `nucmer`, `delta-filter`, `show-coords` and `show-snps` commands. 

Arguments:

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

###pymummer promer class###

[TODO]

###pymummer coords_file class###

Parses the nucmer output and populate an alignment object for each hit in the output

  
###pymummer alignment class###

Check attributes of a hit, swap the reference and query, check if it's a self hit and so on

Contact
-------

Authors: Martin Hunt, Nishadi De Silva

Affiliation: Wellcome Trust Sanger Institute, Hinxton, UK

Email: path-help@sanger.ac.uk
      
