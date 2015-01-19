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
        '''Create a Variant object from a pymummer.snp.Snp object'''
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


    def __eq__(self, other):
        return type(other) is type(self) and self.__dict__ == other.__dict__


    def __str__(self):
        return '\t'.join([
            str(self.ref_start + 1),
            str(self.ref_end + 1),
            str(self.ref_length),
            str(self.ref_name),
            self.ref_base,
            str(self.qry_start + 1),
            str(self.qry_end + 1),
            str(self.qry_length),
            str(self.qry_name),
            self.qry_base
        ])

    def update_indel(self, nucmer_snp):
        '''Indels are reported over multiple lines, 1 base insertion or deletion per line. This method extends the current variant by 1 base if it's an indel and adjacent to the new SNP and returns True. If the current variant is a SNP, does nothing and returns False'''
        new_variant = Variant(nucmer_snp)
        if self.var_type not in [INS, DEL] \
          or self.var_type != new_variant.var_type \
          or self.qry_name != new_variant.qry_name \
          or self.ref_name != new_variant.ref_name:
            return False
        if self.var_type == INS \
          and self.ref_start == new_variant.ref_start \
          and self.qry_end + 1 == new_variant.qry_start:
            self.qry_base += new_variant.qry_base
            self.qry_end += 1
            return True
        if self.var_type == DEL \
          and self.qry_start == new_variant.qry_start \
          and self.ref_end + 1 == new_variant.ref_start:
            self.ref_base += new_variant.ref_base
            self.ref_end += 1
            return True

        return False

