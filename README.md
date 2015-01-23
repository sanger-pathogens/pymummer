pymummer
========

Python3 module for running MUMmer and parsing the output. 

Installation
------------

#Pre-requisites#

__1.	MUMmer__

Instructions to install MUMmer can be found [here](http://mummer.sourceforge.net/manual/#installation)
    
__2.	nose__
	
Install using your chosen method of installing python modules. Instructions [here](https://nose.readthedocs.org/en/latest/). We used pip
	
	pip3 install nose

__3.	fastaq__ 
	
Get a clone of the code from [Github](https://github.com/sanger-pathogens/Fastaq):
	
	git clone git@github.com:sanger-pathogens/Fastaq.git - replace with where tar is
		
Run the tests (needs nose): 
	
	python3 setup.py test
		
Install: 
	
	python3 setup.py install
		
#Installation#

Download code:
	
	replace with where tar is
		
Run the tests (needs nose): 
	
	python3 setup.py test
		
Install: 
	
	python3 setup.py install


Usage (for developers)
----------------------

Example showing how you can run nucmer on a fasta file and parse the output file to produce a set of alignment objects:

	from pymummer import coords_file, alignment, nucmer
	...
	runner = nucmer.Runner(reference_file, query_file, results_file) 
	runner.run()
	file_reader = coords_file.reader(results_file)
	alignments = [coord for coord in file_reader if not coord.is_self_hit()] #Remove self hits
	...

__pymummer nucmer__

Wraps the `nucmer`, `delta-filter`, `show-coords` and `show-snps` commands. 

Arguments:

__ref__			reference file  
__query__			query file  
__outfile__		output file  
__min_id__		min\_id for delta-filter command (default None)  
__min_length__	min\_length for delta-filter command (default None)  
__breaklen__		breaklen for nucmer command (nucmer's default is 200)   
__coords_header__	print header in show-coords output (default True)  
__maxmatch__		maxmatch for nucmer (default False)  
__show_snps__		run show-snps (default False)  
__snps_header__ 	print header in show-snps output (default True)  

__pymummer promer__

[TODO]

__pymummer coords_file__

Parses the nucmer output and populate an alignment object for each hit in the output

  
__pymummer alignment__

Do useful things with each hit like swap the reference and query, 
check if it's a self hit and if the reference & query are on the same strand 
      
