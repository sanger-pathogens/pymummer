import os
import tempfile
import shutil
import fastaq

class Error (Exception): pass

def run_nucmer(
  query,
  ref,
  outfile,
  min_id=95,
  min_length=100,
  breaklen=200,
  ):
    query = os.path.abspath(query)
    ref = os.path.abspath(ref)
    outfile = os.path.abspath(outfile)
    tmpdir = tempfile.mkdtemp(prefix='tmp.run_nucmer.', dir=os.getcwd())
    original_dir = os.getcwd()
    os.chdir(tmpdir)
    script = 'run_nucmer.sh'
    f = fastaq.utils.open_file_write(script)
    print('nucmer --maxmatch -p p -b', breaklen, ref, query, file=f)
    print('delta-filter -i', min_id, '-l', min_length, 'p.delta > p.delta.filter', file=f)
    print('show-coords -dTlro p.delta.filter >', outfile, file=f)
    fastaq.utils.close(f)
    common.syscall('bash ' + script)
    os.chdir(original_dir)
    shutil.rmtree(tmpdir)


def file_reader(fname):
    f = fastaq.utils.open_file_read(fname)
    in_header = True

    for line in f:
        if in_header:
            if line.startswith('['):
                in_header = False
            continue
        yield NucmerHit(line)

    fastaq.utils.close(f)


class NucmerHit:
    def __init__(self, line):
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
        self.ref_start, self.qry_start = self.qry_start, self.ref_start
        self.ref_end, self.qry_end = self.qry_end, self.ref_end
        self.hit_length_ref, self.hit_length_qry = self.hit_length_qry, self.hit_length_ref
        self.ref_length, self.qry_length = self.qry_length, self.ref_length
        self.ref_name, self.qry_name = self.qry_name, self.ref_name


    def qry_coords(self):
        return fastaq.intervals.Interval(min(self.qry_start, self.qry_end), max(self.qry_start, self.qry_end))


    def ref_coords(self):
        return fastaq.intervals.Interval(min(self.ref_start, self.ref_end), max(self.ref_start, self.ref_end))


    def on_same_strand(self):
        return (self.ref_start < self.ref_end) == (self.qry_start < self.qry_end)


    def sort(self):
        if self.ref_name > self.qry_name:
            self._swap()

        if self.ref_start > self.ref_end:
            self.ref_start, self.ref_end = self.ref_end, self.ref_start
            self.qry_start, self.qry_end = self.qry_end, self.qry_start


    def is_self_hit(self):
        return self.ref_name == self.qry_name \
                and self.ref_start == self.qry_start \
                and self.ref_end == self.qry_end \
                and self.percent_identity == 100


    def __str__(self):
        return '\t'.join(str(x) for x in
            [self.ref_start,
            self.ref_end,
            self.qry_start,
            self.qry_end,
            self.hit_length_ref,
            self.hit_length_qry,
            '{0:.2f}'.format(self.percent_identity),
            self.ref_length,
            self.qry_length,
            self.frame,
            self.ref_name,
            self.qry_name])

