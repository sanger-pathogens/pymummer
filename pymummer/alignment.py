import pyfastaq

class Error (Exception): pass

class Alignment:
    def __init__(self, line):
        '''Constructs Alignment object from a line of show-coords -dTlro'''
        # [S1]  [E1]    [S2]    [E2]    [LEN 1] [LEN 2] [% IDY] [LEN R] [LEN Q] [FRM]   [TAGS]
        #1162    25768   24536   4   24607   24533   99.32   640851  24536   1   -1  MAL1    NODE_25757_length_24482_cov_18.920391   [CONTAINS]

        try:
            l = line.rstrip().split('\t')
            self.ref_start = int(l[0]) - 1
            self.ref_end = int(l[1]) - 1
            self.qry_start = int(l[2]) - 1
            self.qry_end = int(l[3]) - 1
            self.hit_length_ref = int(l[4])
            self.hit_length_qry = int(l[5])
            self.percent_identity = float(l[6])
            self.ref_length = int(l[7])
            self.qry_length = int(l[8])
            self.frame = int(l[9])
            self.ref_name = l[11]
            self.qry_name = l[12]

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

