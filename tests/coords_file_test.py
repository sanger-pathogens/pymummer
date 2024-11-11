import unittest
import os
import filecmp
from pymummer import coords_file, alignment

data_dir = "tests/data"


class TestCoordsFile(unittest.TestCase):
    def test_coords_file(self):
        """test coords_file"""
        expected = [
            "\t".join(
                [
                    "61",
                    "900",
                    "1",
                    "840",
                    "840",
                    "840",
                    "99.76",
                    "1000",
                    "840",
                    "1",
                    "1",
                    "test_ref1",
                    "test_qry1",
                    "[CONTAINS]",
                ]
            ),
            "\t".join(
                [
                    "62",
                    "901",
                    "2",
                    "841",
                    "841",
                    "850",
                    "99.66",
                    "999",
                    "839",
                    "1",
                    "1",
                    "test_ref2",
                    "test_qry2",
                    "[CONTAINS]",
                ]
            ),
            "\t".join(
                [
                    "63",
                    "902",
                    "3",
                    "842",
                    "842",
                    "860",
                    "99.56",
                    "998",
                    "838",
                    "1",
                    "1",
                    "test_ref3",
                    "test_qry3",
                    "[CONTAINS]",
                ]
            ),
        ]
        expected = [alignment.Alignment(x) for x in expected]

        infiles = [
            os.path.join(data_dir, "coords_file_test_with_header.coords"),
            os.path.join(data_dir, "coords_file_test_no_header.coords"),
        ]

        for fname in infiles:
            fr = coords_file.reader(fname)
            alignments = [x for x in fr]
            self.assertEqual(alignments, expected)

    def test_convert_to_msp_crunch_no_offset(self):
        """Test convert_to_msp_crunch with no offsets"""
        infile = os.path.join(data_dir, "coords_file_test_convert_to_msp_crunch.coords")
        expected = os.path.join(
            data_dir, "coords_file_test_convert_to_msp_crunch.no_offset.crunch"
        )
        tmpfile = "tmp.test_convert_to_msp_crunch_no_offset.crunch"
        coords_file.convert_to_msp_crunch(infile, tmpfile)
        self.assertTrue(filecmp.cmp(expected, tmpfile, shallow=False))
        os.unlink(tmpfile)

    def test_convert_to_msp_crunch_with_offset(self):
        """Test convert_to_msp_crunch with offsets"""
        infile = os.path.join(data_dir, "coords_file_test_convert_to_msp_crunch.coords")
        ref_fai = os.path.join(
            data_dir, "coords_file_test_convert_to_msp_crunch.ref.fa.fai"
        )
        qry_fai = os.path.join(
            data_dir, "coords_file_test_convert_to_msp_crunch.qry.fa.fai"
        )
        expected = os.path.join(
            data_dir, "coords_file_test_convert_to_msp_crunch.with_offset.crunch"
        )
        tmpfile = "tmp.test_convert_to_msp_crunch_with_offset.crunch"

        with self.assertRaises(coords_file.Error):
            coords_file.convert_to_msp_crunch(infile, tmpfile, ref_fai=ref_fai)
            coords_file.convert_to_msp_crunch(infile, tmpfile, qry_fai=qry_fai)

        coords_file.convert_to_msp_crunch(
            infile, tmpfile, ref_fai=ref_fai, qry_fai=qry_fai
        )
        self.assertTrue(filecmp.cmp(expected, tmpfile, shallow=False))
        os.unlink(tmpfile)
