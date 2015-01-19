import unittest
import os
import fastaq
from pymummer import snp_file, snp

modules_dir = os.path.dirname(os.path.abspath(snp_file.__file__))
data_dir = os.path.join(modules_dir, 'tests', 'data')

class TestUtils(unittest.TestCase):
    def test_snp_file(self):
        '''test coords_file'''
        expected = [
            '\t'.join(['133', 'G', '.', '122', '1', '122', '500', '489', '1', '1', 'ref', 'qry']),
            '\t'.join(['143', '.', 'C', '131', '1', '132', '500', '489', '1', '1', 'ref', 'qry']),
            '\t'.join(['253', 'T', 'A', '242', '120', '242', '500', '489', '1', '1', 'ref', 'qry'])
        ]

        expected = [snp.Snp(x) for x in expected]

        infiles = [os.path.join(data_dir, 'snp_file_test_with_header.snps'), os.path.join(data_dir, 'snp_file_test_no_header.snps')]

        for fname in infiles:
            fr = snp_file.reader(fname)
            snps = [x for x in fr]
            self.assertEqual(snps, expected)

