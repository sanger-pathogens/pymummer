import os
import tempfile
import shutil
import pyfastaq
from pymummer import syscall

class Error (Exception): pass


class Runner:
    '''Handy reference for all the arguments needed for nucmer, delta-filter, show-coords, show-snps'''
    def __init__(
      self,
      ref,
      query,
      outfile,
      min_id=None,
      min_length=None,
      breaklen=None,
      coords_header=True,
      diagdiff=None,
      maxmatch=False,
      mincluster=None,
      simplify=True,
      show_snps=False,
      snps_header=True,
      verbose=False,
      promer=False,
   ):
        self.qry = query
        self.ref = ref
        self.outfile = outfile
        self.min_id = min_id
        self.min_length = min_length
        self.breaklen = breaklen
        self.diagdiff = diagdiff
        self.coords_header = coords_header
        self.maxmatch = maxmatch
        self.mincluster = mincluster
        self.simplify = simplify
        self.show_snps = show_snps
        self.snps_header = snps_header
        self.verbose = verbose
        self.use_promer = promer



    def _nucmer_command(self, ref, qry, outprefix):
        '''Construct the nucmer command'''
        if self.use_promer:
            command = 'promer'
        else:
            command = 'nucmer'

        command += ' -p ' + outprefix

        if self.breaklen is not None:
            command += ' -b ' + str(self.breaklen)

        if self.diagdiff is not None and not self.use_promer:
            command += ' -D ' + str(self.diagdiff)

        if self.maxmatch:
            command += ' --maxmatch'

        if self.mincluster is not None:
            command += ' -c ' + str(self.mincluster)

        if not self.simplify and not self.use_promer:
        	command += ' --nosimplify'

        return command + ' ' + ref + ' ' + qry



    def _delta_filter_command(self, infile, outfile):
        '''Construct delta-filter command'''
        command = 'delta-filter'

        if self.min_id is not None:
            command += ' -i ' + str(self.min_id)

        if self.min_length is not None:
            command += ' -l ' + str(self.min_length)

        return command + ' ' + infile + ' > ' + outfile


    def _show_coords_command(self, infile, outfile):
        '''Construct show-coords command'''
        command = 'show-coords -dTlro'

        if not self.coords_header:
            command += ' -H'

        return command + ' ' + infile + ' > ' + outfile


    def _show_snps_command(self, infile, outfile):
        command = 'show-snps -TClr'

        if not self.snps_header:
            command += ' -H'

        return command + ' ' + infile + ' > ' + outfile


    def _write_script(self, script_name, ref, qry, outfile):
        '''Write commands into a bash script'''
        f = pyfastaq.utils.open_file_write(script_name)
        print(self._nucmer_command(ref, qry, 'p'), file=f)
        print(self._delta_filter_command('p.delta', 'p.delta.filter'), file=f)
        print(self._show_coords_command('p.delta.filter', outfile), file=f)
        if self.show_snps:
            print(self._show_snps_command('p.delta.filter', outfile + '.snps'), file=f)
        pyfastaq.utils.close(f)


    def run(self):
        '''
        Change to a temp directory
        Run bash script containing commands
        Place results in specified output file
        Clean up temp directory
        '''
        qry = os.path.abspath(self.qry)
        ref = os.path.abspath(self.ref)
        outfile = os.path.abspath(self.outfile)
        tmpdir = tempfile.mkdtemp(prefix='tmp.run_nucmer.', dir=os.getcwd())
        original_dir = os.getcwd()
        os.chdir(tmpdir)
        script = 'run_nucmer.sh'
        self._write_script(script, ref, qry, outfile)
        syscall.run('bash ' + script, verbose=self.verbose)
        os.chdir(original_dir)
        shutil.rmtree(tmpdir)

