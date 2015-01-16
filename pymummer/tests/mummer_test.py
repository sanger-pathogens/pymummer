import unittest
import os
import fastaq
from pymummer import mummer

modules_dir = os.path.dirname(os.path.abspath(mummer.__file__))
data_dir = os.path.join(modules_dir, 'tests', 'data')

class TestMummer(unittest.TestCase):
    def test_run_nucmer(self):
        '''test run_nucmer'''
        # TODO
        pass


    def test_swap(self):
        ''' test swap'''
        # TODO
        pass


    def test_sort(self):
        '''test sort'''
        # TODO
        pass


    def test_file_read(self):
        '''test file_read'''
        # TODO
        pass


    def test_qry_coords(self):
        '''Test qry_coords'''
        hits = ['\t'.join(['1', '100', '1', '100', '100', '100', '100.00', '1000', '1000', '1', '1', 'ref', 'qry']),
                '\t'.join(['1', '101', '100', '1', '100', '100', '100.00', '1000', '1000', '1', '1', 'ref', 'qry'])
        ]
        for h in hits:
            m = mummer.NucmerHit(h)
            self.assertEqual(fastaq.intervals.Interval(0,99), m.qry_coords())


    def test_ref_coords(self):
        '''Test ref_coords'''
        hits = ['\t'.join(['1', '100', '1', '100', '100', '100', '100.00', '1000', '1000', '1', '1', 'ref', 'ref']),
                '\t'.join(['100', '1', '100', '1', '100', '100', '100.00', '1000', '1000', '1', '1', 'ref', 'ref'])
        ]
        for h in hits:
            m = mummer.NucmerHit(h)
            self.assertEqual(fastaq.intervals.Interval(0,99), m.ref_coords())


    def test_on_same_strand(self):
        '''test on_same_strand'''
        self.assertTrue(mummer.NucmerHit('\t'.join(['1', '100', '1', '100', '100', '100', '100.00', '1000', '1000', '1', '1', 'ref', 'ref'])).on_same_strand())
        self.assertTrue(mummer.NucmerHit('\t'.join(['100', '1', '100', '1', '100', '100', '100.00', '1000', '1000', '1', '1', 'ref', 'ref'])).on_same_strand())
        self.assertFalse(mummer.NucmerHit('\t'.join(['1', '100', '100', '1', '100', '100', '100.00', '1000', '1000', '1', '1', 'ref', 'ref'])).on_same_strand())
        self.assertFalse(mummer.NucmerHit('\t'.join(['100', '1', '1', '100', '100', '100', '100.00', '1000', '1000', '1', '1', 'ref', 'ref'])).on_same_strand())


    def test_is_self_hit(self):
        '''Test is_self_hit'''
        tests = [('\t'.join(['1', '100', '1', '100', '100', '100', '100.00', '1000', '1000', '1', '1', 'ref', 'ref']), True),
            ('\t'.join(['1', '101', '1', '100', '100', '100', '100.00', '1000', '1000', '1', '1', 'ref', 'ref']), False),
            ('\t'.join(['2', '100', '1', '100', '100', '100', '100.00', '1000', '1000', '1', '1', 'ref', 'ref']), False),
            ('\t'.join(['1', '100', '1', '100', '100', '100', '100.00', '1000', '1000', '1', '1', 'ref', 'ref2']), False),
            ('\t'.join(['1', '100', '1', '100', '100', '100', '99.9', '1000', '1000', '1', '1', 'ref', 'ref']), False),
        ]

        for t in tests:
            m = mummer.NucmerHit(t[0])
            self.assertEqual(m.is_self_hit(), t[1])
        pass



