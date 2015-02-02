import os
import sys
import subprocess

class Error (Exception): pass


def decode(x):
    try:
        s = x.decode()
    except:
        return x
    return s


def run(cmd, verbose=False):
    if verbose:
        print('Running command:', cmd, flush=True)
    try:
        output = subprocess.check_output(cmd, shell=True, stderr=subprocess.STDOUT)
    except subprocess.CalledProcessError as error:
        print('The following command failed with exit code', error.returncode, file=sys.stderr)
        print(cmd, file=sys.stderr)
        print('\nThe output was:\n', file=sys.stderr)
        print(error.output.decode(), file=sys.stderr)
        raise Error('Error running command:', cmd)
            
    if verbose:
        print(decode(output))

