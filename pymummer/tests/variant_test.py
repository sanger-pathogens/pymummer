import unittest
import os

from pymummer import variant, snp

modules_dir = os.path.dirname(os.path.abspath(variant.__file__))
data_dir = os.path.join(modules_dir, 'tests', 'data')


class TestVariant(unittest.TestCase):
    def test_init(self):
        '''Test init gets correct vairnat type'''
        lines = [
            ['42', 'T', 'A', '42', '42', '42', '1000', '1000', '1', '1', 'ref', 'ref'],
            ['242', 'G', '.', '241', '1', '241', '1000', '1000', '1', '1', 'ref', 'ref'],
            ['300', '.', 'G', '298', '0', '298', '1000', '1000', '1', '1', 'ref', 'ref']
        ]

        variants = [variant.Variant(snp.Snp('\t'.join(x))) for x in lines]
        expected = [variant.SNP, variant.DEL, variant.INS]
        for i in range(len(lines)):
            self.assertEqual(variants[i].var_type, expected[i])

