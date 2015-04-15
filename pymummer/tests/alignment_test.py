import unittest
import os
import pyfastaq
from pymummer import alignment

modules_dir = os.path.dirname(os.path.abspath(alignment.__file__))
data_dir = os.path.join(modules_dir, 'tests', 'data')

class TestNucmer(unittest.TestCase):
    def test_init_nucmer(self):
        '''test __init__ nucmer'''
        line = '\t'.join(['1', '100', '2', '200', '101', '202', '42.42', '123', '456', '-1', '0', 'ref', 'qry', '[FOO]'])
        a = alignment.Alignment(line)
        self.assertEqual(a.ref_start, 0)
        self.assertEqual(a.ref_end, 99)
        self.assertEqual(a.qry_start, 1)
        self.assertEqual(a.qry_end, 199)
        self.assertEqual(a.hit_length_ref, 101)
        self.assertEqual(a.hit_length_qry, 202)
        self.assertEqual(a.percent_identity, 42.42)
        self.assertEqual(a.ref_length, 123)
        self.assertEqual(a.qry_length, 456)
        self.assertEqual(a.frame, -1)
        self.assertEqual(a.ref_name, 'ref')
        self.assertEqual(a.qry_name, 'qry')


    def test_init_promer(self):
        '''test __init__ promer'''
        line = '\t'.join(['1', '1398', '4891054', '4892445', '1398', '1392', '89.55', '93.18', '0.21', '1398', '5349013', '1', '1', 'ref', 'qry', '[CONTAINED]'])

        a = alignment.Alignment(line)
        self.assertEqual(a.ref_start, 0)
        self.assertEqual(a.ref_end, 1397)
        self.assertEqual(a.qry_start, 4891053)
        self.assertEqual(a.qry_end, 4892444)
        self.assertEqual(a.hit_length_ref, 1398)
        self.assertEqual(a.hit_length_qry, 1392)
        self.assertEqual(a.percent_identity, 89.55)
        self.assertEqual(a.ref_length, 1398)
        self.assertEqual(a.qry_length, 5349013)
        self.assertEqual(a.frame, 1)
        self.assertEqual(a.ref_name, 'ref')
        self.assertEqual(a.qry_name, 'qry')


    def test_swap(self):
        ''' test swap'''
        l_in = ['1', '100', '2', '200', '101', '202', '42.42', '123', '456', '-1', '0', 'ref', 'qry']
        l_out = ['2', '200', '1', '100', '202', '101', '42.42', '456', '123', '-1', '0', 'qry', 'ref']
        a_in = alignment.Alignment('\t'.join(l_in))
        a_in._swap()
        self.assertEqual(a_in, alignment.Alignment('\t'.join(l_out)))
        a_in._swap()
        self.assertEqual(a_in, alignment.Alignment('\t'.join(l_in)))


    def test_qry_coords(self):
        '''Test qry_coords'''
        hits = ['\t'.join(['1', '100', '1', '100', '100', '100', '100.00', '1000', '1000', '1', '1', 'ref', 'qry']),
                '\t'.join(['1', '101', '100', '1', '100', '100', '100.00', '1000', '1000', '1', '1', 'ref', 'qry'])
        ]
        for h in hits:
            a = alignment.Alignment(h)
            self.assertEqual(pyfastaq.intervals.Interval(0,99), a.qry_coords())


    def test_ref_coords(self):
        '''Test ref_coords'''
        hits = ['\t'.join(['1', '100', '1', '100', '100', '100', '100.00', '1000', '1000', '1', '1', 'ref', 'ref']),
                '\t'.join(['100', '1', '100', '1', '100', '100', '100.00', '1000', '1000', '1', '1', 'ref', 'ref'])
        ]
        for h in hits:
            a = alignment.Alignment(h)
            self.assertEqual(pyfastaq.intervals.Interval(0,99), a.ref_coords())


    def test_on_same_strand(self):
        '''test on_same_strand'''
        self.assertTrue(alignment.Alignment('\t'.join(['1', '100', '1', '100', '100', '100', '100.00', '1000', '1000', '1', '1', 'ref', 'ref'])).on_same_strand())
        self.assertTrue(alignment.Alignment('\t'.join(['100', '1', '100', '1', '100', '100', '100.00', '1000', '1000', '1', '1', 'ref', 'ref'])).on_same_strand())
        self.assertFalse(alignment.Alignment('\t'.join(['1', '100', '100', '1', '100', '100', '100.00', '1000', '1000', '1', '1', 'ref', 'ref'])).on_same_strand())
        self.assertFalse(alignment.Alignment('\t'.join(['100', '1', '1', '100', '100', '100', '100.00', '1000', '1000', '1', '1', 'ref', 'ref'])).on_same_strand())


    def test_is_self_hit(self):
        '''Test is_self_hit'''
        tests = [('\t'.join(['1', '100', '1', '100', '100', '100', '100.00', '1000', '1000', '1', '1', 'ref', 'ref']), True),
            ('\t'.join(['1', '101', '1', '100', '100', '100', '100.00', '1000', '1000', '1', '1', 'ref', 'ref']), False),
            ('\t'.join(['2', '100', '1', '100', '100', '100', '100.00', '1000', '1000', '1', '1', 'ref', 'ref']), False),
            ('\t'.join(['1', '100', '1', '100', '100', '100', '100.00', '1000', '1000', '1', '1', 'ref', 'ref2']), False),
            ('\t'.join(['1', '100', '1', '100', '100', '100', '99.9', '1000', '1000', '1', '1', 'ref', 'ref']), False),
        ]

        for t in tests:
            a = alignment.Alignment(t[0])
            self.assertEqual(a.is_self_hit(), t[1])


    def test_str(self):
        '''Test __str__'''
        l_in = ['1', '100', '2', '200', '101', '202', '42.42', '123', '456', '-1', '0', 'ref', 'qry']
        # the 10th column (counting from zero) is ignored and so not output by __str__
        l_out = l_in[:10] + l_in[11:]
        a = alignment.Alignment('\t'.join(l_in))
        self.assertEqual(str(a), '\t'.join(l_out))

