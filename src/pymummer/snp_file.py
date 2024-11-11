import pyfastaq
from pymummer import snp, variant

def reader(fname):
    f = pyfastaq.utils.open_file_read(fname)

    for line in f:
        if line.startswith('[') or (not '\t' in line):
            continue

        yield snp.Snp(line)

    pyfastaq.utils.close(f)


def get_all_variants(fname):
    variants = []
    fr = reader(fname)
    for nucmer_snp in fr:
        if len(variants) == 0 or not variants[-1].update_indel(nucmer_snp):
            variants.append(variant.Variant(nucmer_snp))

    return variants

