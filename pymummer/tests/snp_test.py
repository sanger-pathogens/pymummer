import unittest
import os
from pymummer import snp


modules_dir = os.path.dirname(os.path.abspath(snp.__file__))
data_dir = os.path.join(modules_dir, 'tests', 'data')


class TestSnp(unittest.TestCase):
    def test_str_no_c_option(self):
        '''Test __str__ with format with no -C option'''
        l_in = ['187', 'A', 'C', '269', '187', '187', '654', '853', '1', '1', 'ref_name', 'qry_name']
        s = snp.Snp('\t'.join(l_in))
        expected = '\t'.join(['187', 'A', 'C', '269', '654', '853', 'ref_name', 'qry_name'])
        self.assertEqual(str(s), expected)


    def test_str_with_c_option(self):
        '''Test __str__ with format with -C option'''
        l_in = ['187', 'A', 'C', '269', '187', '187', '0', '0', '654', '853', '1', '1', 'ref_name', 'qry_name']
        s = snp.Snp('\t'.join(l_in))
        expected = '\t'.join(['187', 'A', 'C', '269', '654', '853', 'ref_name', 'qry_name'])
        self.assertEqual(str(s), expected)
