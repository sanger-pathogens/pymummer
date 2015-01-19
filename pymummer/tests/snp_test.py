import unittest
import os
from pymummer import snp


modules_dir = os.path.dirname(os.path.abspath(snp.__file__))
data_dir = os.path.join(modules_dir, 'tests', 'data')


class TestSnp(unittest.TestCase):
    def test_str(self):
        '''Test __str__'''
        l_in = ['187', 'A', 'C', '269', '187', '187', '654', '853', '1', '1', 'ref_name', 'qry_name']
        # only use columns 0-3, 6-7, 10-11
        l_out = l_in[:4] + l_in[6:8] + l_in[10:]
        s = snp.Snp('\t'.join(l_in))
        self.assertEqual(str(s), '\t'.join(l_out))
