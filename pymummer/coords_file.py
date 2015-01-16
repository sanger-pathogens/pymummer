import fastaq
from pymummer import alignment

def reader(fname):
    '''Helper function to open the results file (coords file) and create alignment objects with the values in it'''
    f = fastaq.utils.open_file_read(fname)

    for line in f:
        if line.startswith('[') or (not '\t' in line):
            continue

        yield alignment.Alignment(line)

    fastaq.utils.close(f)

