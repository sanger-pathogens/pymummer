import fastaq
from pymummer import snp

def reader(fname):
    f = fastaq.utils.open_file_read(fname)

    for line in f:
        if line.startswith('[') or (not '\t' in line):
            continue

        yield snp.Snp(line)

    fastaq.utils.close(f)

