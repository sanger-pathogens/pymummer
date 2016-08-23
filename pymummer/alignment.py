import pyfastaq
from pymummer import variant

class Error (Exception): pass

class Alignment:
    def __init__(self, line):
        '''Constructs Alignment object from a line of show-coords -dTlro'''
        # nucmer:
        # [S1]  [E1]    [S2]    [E2]    [LEN 1] [LEN 2] [% IDY] [LEN R] [LEN Q] [FRM]   [TAGS]
        #1162    25768   24536   4   24607   24533   99.32   640851  24536   1   -1  ref qry   [CONTAINS]

        # promer:
        #[S1]    [E1]    [S2]    [E2]    [LEN 1] [LEN 2] [% IDY] [% SIM] [% STP] [LEN R] [LEN Q] [FRM]   [TAGS]
        # 1   1398    4891054 4892445 1398    1392    89.55   93.18   0.21    1398    5349013 1   1   ref qry    [CONTAINED]

        fields = line.rstrip().split('\t')

        try:
            self.ref_start = int(fields[0]) - 1
            self.ref_end = int(fields[1]) - 1
            self.qry_start = int(fields[2]) - 1
            self.qry_end = int(fields[3]) - 1
            self.hit_length_ref = int(fields[4])
            self.hit_length_qry = int(fields[5])
            self.percent_identity = float(fields[6])

            if len(fields) >= 15:  # promer has more fields
                self.ref_length = int(fields[9])
                self.qry_length = int(fields[10])
                self.frame = int(fields[11])
                self.ref_name = fields[13]
                self.qry_name = fields[14]
            else:
                self.ref_length = int(fields[7])
                self.qry_length = int(fields[8])
                self.frame = int(fields[9])
                self.ref_name = fields[11]
                self.qry_name = fields[12]
        except:
            raise Error('Error reading this nucmer line:\n' + line)


    def __eq__(self, other):
        return type(other) is type(self) and self.__dict__ == other.__dict__


    def __hash__(self):
        return hash((self.ref_start, self.ref_end, self.qry_start, self.qry_end, self.hit_length_ref, self.hit_length_qry, self.percent_identity, self.ref_length, self.qry_length, self.frame, self.ref_name, self.qry_name))


    def _swap(self):
        '''Swaps the alignment so that the reference becomes the query and vice-versa. Swaps their names, coordinates etc. The frame is not changed'''
        self.ref_start, self.qry_start = self.qry_start, self.ref_start
        self.ref_end, self.qry_end = self.qry_end, self.ref_end
        self.hit_length_ref, self.hit_length_qry = self.hit_length_qry, self.hit_length_ref
        self.ref_length, self.qry_length = self.qry_length, self.ref_length
        self.ref_name, self.qry_name = self.qry_name, self.ref_name


    def qry_coords(self):
        '''Returns a pyfastaq.intervals.Interval object of the start and end coordinates in the query sequence'''
        return pyfastaq.intervals.Interval(min(self.qry_start, self.qry_end), max(self.qry_start, self.qry_end))


    def ref_coords(self):
        '''Returns a pyfastaq.intervals.Interval object of the start and end coordinates in the reference sequence'''
        return pyfastaq.intervals.Interval(min(self.ref_start, self.ref_end), max(self.ref_start, self.ref_end))


    def on_same_strand(self):
        '''Returns true iff the direction of the alignment is the same in the reference and the query'''
        return (self.ref_start < self.ref_end) == (self.qry_start < self.qry_end)


    def is_self_hit(self):
        '''Returns true iff the alignment is of a sequence to itself: names and all coordinates are the same and 100 percent identity'''
        return self.ref_name == self.qry_name \
                and self.ref_start == self.qry_start \
                and self.ref_end == self.qry_end \
                and self.percent_identity == 100


    def reverse_query(self):
        '''Changes the coordinates as if the query sequence has been reverse complemented'''
        self.qry_start = self.qry_length - self.qry_start - 1
        self.qry_end = self.qry_length - self.qry_end - 1


    def reverse_reference(self):
        '''Changes the coordinates as if the reference sequence has been reverse complemented'''
        self.ref_start = self.ref_length - self.ref_start - 1
        self.ref_end = self.ref_length - self.ref_end - 1


    def __str__(self):
        '''Returns a tab delimited string containing the values of this alignment object'''
        return '\t'.join(str(x) for x in
            [self.ref_start + 1,
            self.ref_end + 1,
            self.qry_start + 1,
            self.qry_end + 1,
            self.hit_length_ref,
            self.hit_length_qry,
            '{0:.2f}'.format(self.percent_identity),
            self.ref_length,
            self.qry_length,
            self.frame,
            self.ref_name,
            self.qry_name])


    def to_msp_crunch(self):
        '''Returns the alignment as a line in MSPcrunch format. The columns are space-separated and are:
           1. score
           2. percent identity
           3. match start in the query sequence
           4. match end in the query sequence
           5. query sequence name
           6. subject sequence start
           7. subject sequence end
           8. subject sequence name'''

        # we don't know the alignment score. Estimate it. This approximates 1 for a match.
        aln_score = int(self.percent_identity * 0.005 * (self.hit_length_ref + self.hit_length_qry))

        return ' '.join(str(x) for x in [
                aln_score,
                '{0:.2f}'.format(self.percent_identity),
                self.qry_start + 1,
                self.qry_end + 1,
                self.qry_name,
                self.ref_start + 1,
                self.ref_end + 1,
                self.ref_name
        ])


    def qry_coords_from_ref_coord(self, ref_coord, variant_list):
        '''Given a reference position and a list of variants ([variant.Variant]),
           works out the position in the query sequence, accounting for indels.
           Returns a tuple: (position, True|False), where second element is whether
           or not the ref_coord lies in an indel. If it is, then
           returns the corresponding start position
           of the indel in the query'''
        if self.ref_coords().distance_to_point(ref_coord) > 0:
            raise Error('Cannot get query coord in qry_coords_from_ref_coord because given ref_coord ' + str(ref_coord) + ' does not lie in nucmer alignment:\n' + str(self))

        indel_variant_indexes = []

        for i in range(len(variant_list)):
            if variant_list[i].var_type not in {variant.INS, variant.DEL}:
                continue
            if variant_list[i].ref_start <= ref_coord <= variant_list[i].ref_end:
                return variant_list[i].qry_start, True
            elif variant_list[i].ref_start < ref_coord:
                indel_variant_indexes.append(i)

        distance = ref_coord - min(self.ref_start, self.ref_end)

        for i in indel_variant_indexes:
            if variant_list[i].var_type == variant.INS:
                distance += len(variant_list[i].qry_base)
            else:
                assert variant_list[i].var_type == variant.DEL
                distance -= len(variant_list[i].ref_base)

        if self.on_same_strand():
            return min(self.qry_start, self.qry_end) + distance, False
        else:
            return max(self.qry_start, self.qry_end) - distance, False


    def ref_coords_from_qry_coord(self, qry_coord, variant_list):
        '''Given a qryerence position and a list of variants ([variant.Variant]),
           works out the position in the ref sequence, accounting for indels.
           Returns a tuple: (position, True|False), where second element is whether
           or not the qry_coord lies in an indel. If it is, then
           returns the corresponding start position
           of the indel in the ref'''
        if self.qry_coords().distance_to_point(qry_coord) > 0:
            raise Error('Cannot get ref coord in ref_coords_from_qry_coord because given qry_coord ' + str(qry_coord) + ' does not lie in nucmer alignment:\n' + str(self))

        indel_variant_indexes = []

        for i in range(len(variant_list)):
            if variant_list[i].var_type not in {variant.INS, variant.DEL}:
                continue
            if variant_list[i].qry_start <= qry_coord <= variant_list[i].qry_end:
                return variant_list[i].ref_start, True
            elif variant_list[i].qry_start < qry_coord:
                indel_variant_indexes.append(i)

        distance = qry_coord - min(self.qry_start, self.qry_end)

        for i in indel_variant_indexes:
            if variant_list[i].var_type == variant.DEL:
                distance += len(variant_list[i].ref_base)
            else:
                assert variant_list[i].var_type == variant.INS
                distance -= len(variant_list[i].qry_base)

        if self.on_same_strand():
            return min(self.ref_start, self.ref_end) + distance, False
        else:
            return max(self.ref_start, self.ref_end) - distance, False

