class Error (Exception): pass

SNP = 1
DEL = 2
INS = 3

var_types = {
    1: 'SNP',
    2: 'DEL',
    3: 'INS',
}


class Variant:
    def __init__(self, snp):
        if snp.ref_base == '.':
            self.var_type = INS
            self.qry_base = snp.qry_base
            self.ref_base = '.'
        elif snp.qry_base == '.':
            self.var_type = DEL
            self.qry_base = '.'
            self.ref_base = snp.ref_base
        elif '.' not in [snp.ref_base, snp.qry_base]:
            self.var_type = SNP
            self.ref_base = snp.ref_base
            self.qry_base = snp.qry_base
        else:
            raise Error('Error constructing Variant from pymummer.snp.Snp:' + str(snp))

        self.ref_start = snp.ref_pos
        self.ref_end = snp.ref_pos
        self.ref_length = snp.ref_length
        self.ref_name = snp.ref_name
        self.qry_start = snp.qry_pos
        self.qry_end = snp.qry_pos
        self.qry_length = snp.qry_length
        self.qry_name = snp.qry_name


