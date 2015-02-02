import unittest
import os
from pymummer import coords_file, alignment

modules_dir = os.path.dirname(os.path.abspath(coords_file.__file__))
data_dir = os.path.join(modules_dir, 'tests', 'data')

class TestUtils(unittest.TestCase):
    def test_coords_file(self):
        '''test coords_file'''
        expected = [
            '\t'.join(['61', '900', '1', '840', '840', '840', '99.76', '1000', '840', '1', '1', 'test_ref1', 'test_qry1', '[CONTAINS]']),
            '\t'.join(['62', '901', '2', '841', '841', '850', '99.66', '999', '839', '1', '1', 'test_ref2', 'test_qry2', '[CONTAINS]']),
            '\t'.join(['63', '902', '3', '842', '842', '860', '99.56', '998', '838', '1', '1', 'test_ref3', 'test_qry3', '[CONTAINS]'])
        ]
        expected = [alignment.Alignment(x) for x in expected]

        infiles = [os.path.join(data_dir, 'coords_file_test_with_header.coords'), os.path.join(data_dir, 'coords_file_test_no_header.coords')]

        for fname in infiles:
            fr = coords_file.reader(fname)
            alignments = [x for x in fr]
            self.assertEqual(alignments, expected)

