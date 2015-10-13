import pyfastaq
from pymummer import alignment

class Error (Exception): pass

def reader(fname):
    '''Helper function to open the results file (coords file) and create alignment objects with the values in it'''
    f = pyfastaq.utils.open_file_read(fname)

    for line in f:
        if line.startswith('[') or (not '\t' in line):
            continue

        yield alignment.Alignment(line)

    pyfastaq.utils.close(f)


def convert_to_msp_crunch(infile, outfile, ref_fai=None, qry_fai=None):
    '''Converts a coords file to a file in MSPcrunch format (for use with ACT, most likely).
       ACT ignores sequence names in the crunch file, and just looks at the numbers.
       To make a compatible file, the coords all must be shifted appropriately, which
       can be done by providing both the ref_fai and qry_fai options.
       Both or neither of these must be used, otherwise an error will be thrown.'''
    fai_files = {ref_fai, qry_fai}
    if None in fai_files and len(fai_files) != 1:
       print(fai_files)
       raise Error('Error in convert_to_msp_crunch. Must use both of ref_fai and qry_fai, or neither of them')

    if ref_fai is not None:
        assert qry_fai is not None
        ref_offsets = pyfastaq.tasks.length_offsets_from_fai(ref_fai)
        qry_offsets = pyfastaq.tasks.length_offsets_from_fai(qry_fai)
    
    file_reader = reader(infile)
    f_out = pyfastaq.utils.open_file_write(outfile)

    for aln in file_reader:
        if ref_fai is not None:
           aln.ref_start += ref_offsets[aln.ref_name] 
           aln.ref_end += ref_offsets[aln.ref_name] 
           aln.qry_start += qry_offsets[aln.qry_name] 
           aln.qry_end += qry_offsets[aln.qry_name] 

        print(aln.to_msp_crunch(), file=f_out)

    pyfastaq.utils.close(f_out)
