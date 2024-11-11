import unittest
import copy
import os
from pymummer import variant, snp

data_dir = "tests/data"


class TestVariant(unittest.TestCase):
    def test_init(self):
        """Test init gets correct variant type"""
        lines = [
            ["42", "T", "A", "42", "42", "42", "1000", "1000", "1", "1", "ref", "ref"],
            [
                "242",
                "G",
                ".",
                "241",
                "1",
                "241",
                "1000",
                "1000",
                "1",
                "1",
                "ref",
                "ref",
            ],
            [
                "300",
                ".",
                "G",
                "298",
                "0",
                "298",
                "1000",
                "1000",
                "1",
                "1",
                "ref",
                "ref",
            ],
        ]

        variants = [variant.Variant(snp.Snp("\t".join(x))) for x in lines]
        expected = [variant.SNP, variant.DEL, variant.INS]
        for i in range(len(lines)):
            self.assertEqual(variants[i].var_type, expected[i])

    def test_update_indel_no_change(self):
        """Test update_indel does nothing in the right cases"""
        initial_vars = [
            snp.Snp(
                "\t".join(
                    [
                        "42",
                        "A",
                        "C",
                        "100",
                        "x",
                        "x",
                        "300",
                        "400",
                        "x",
                        "1",
                        "ref",
                        "qry",
                    ]
                )
            ),
            snp.Snp(
                "\t".join(
                    [
                        "42",
                        "A",
                        "C",
                        "100",
                        "x",
                        "x",
                        "300",
                        "400",
                        "x",
                        "1",
                        "ref",
                        "qry",
                    ]
                )
            ),
            snp.Snp(
                "\t".join(
                    [
                        "42",
                        "A",
                        ".",
                        "100",
                        "x",
                        "x",
                        "300",
                        "400",
                        "x",
                        "1",
                        "ref",
                        "qry",
                    ]
                )
            ),
            snp.Snp(
                "\t".join(
                    [
                        "42",
                        "A",
                        ".",
                        "100",
                        "x",
                        "x",
                        "300",
                        "400",
                        "x",
                        "1",
                        "ref",
                        "qry",
                    ]
                )
            ),
            snp.Snp(
                "\t".join(
                    [
                        "42",
                        "A",
                        ".",
                        "100",
                        "x",
                        "x",
                        "300",
                        "400",
                        "x",
                        "1",
                        "ref",
                        "qry",
                    ]
                )
            ),
            snp.Snp(
                "\t".join(
                    [
                        "42",
                        "A",
                        ".",
                        "100",
                        "x",
                        "x",
                        "300",
                        "400",
                        "x",
                        "1",
                        "ref",
                        "qry",
                    ]
                )
            ),
            snp.Snp(
                "\t".join(
                    [
                        "42",
                        "A",
                        ".",
                        "100",
                        "x",
                        "x",
                        "300",
                        "400",
                        "x",
                        "1",
                        "ref",
                        "qry",
                    ]
                )
            ),
            snp.Snp(
                "\t".join(
                    [
                        "42",
                        ".",
                        "A",
                        "100",
                        "x",
                        "x",
                        "300",
                        "400",
                        "x",
                        "1",
                        "ref",
                        "qry",
                    ]
                )
            ),
            snp.Snp(
                "\t".join(
                    [
                        "42",
                        ".",
                        "A",
                        "100",
                        "x",
                        "x",
                        "300",
                        "400",
                        "x",
                        "1",
                        "ref",
                        "qry",
                    ]
                )
            ),
            snp.Snp(
                "\t".join(
                    [
                        "42",
                        ".",
                        "A",
                        "100",
                        "x",
                        "x",
                        "300",
                        "400",
                        "x",
                        "1",
                        "ref",
                        "qry",
                    ]
                )
            ),
            snp.Snp(
                "\t".join(
                    [
                        "42",
                        ".",
                        "A",
                        "100",
                        "x",
                        "x",
                        "300",
                        "400",
                        "x",
                        "1",
                        "ref",
                        "qry",
                    ]
                )
            ),
            snp.Snp(
                "\t".join(
                    [
                        "42",
                        ".",
                        "A",
                        "100",
                        "x",
                        "x",
                        "300",
                        "400",
                        "x",
                        "1",
                        "ref",
                        "qry",
                    ]
                )
            ),
        ]

        to_add = [
            snp.Snp(
                "\t".join(
                    [
                        "142",
                        "A",
                        ".",
                        "1000",
                        "x",
                        "x",
                        "2000",
                        "3000",
                        "x",
                        "1",
                        "ref",
                        "qry",
                    ]
                )
            ),
            snp.Snp(
                "\t".join(
                    [
                        "142",
                        ".",
                        "A",
                        "1000",
                        "x",
                        "x",
                        "2000",
                        "3000",
                        "x",
                        "1",
                        "ref",
                        "qry",
                    ]
                )
            ),
            snp.Snp(
                "\t".join(
                    [
                        "43",
                        "A",
                        ".",
                        "100",
                        "x",
                        "x",
                        "300",
                        "400",
                        "x",
                        "1",
                        "ref2",
                        "qry",
                    ]
                )
            ),
            snp.Snp(
                "\t".join(
                    [
                        "43",
                        "A",
                        ".",
                        "100",
                        "x",
                        "x",
                        "300",
                        "400",
                        "x",
                        "1",
                        "ref",
                        "qry2",
                    ]
                )
            ),
            snp.Snp(
                "\t".join(
                    [
                        "44",
                        "A",
                        ".",
                        "100",
                        "x",
                        "x",
                        "300",
                        "400",
                        "x",
                        "1",
                        "ref",
                        "qry",
                    ]
                )
            ),
            snp.Snp(
                "\t".join(
                    [
                        "42",
                        "A",
                        ".",
                        "100",
                        "x",
                        "x",
                        "300",
                        "400",
                        "x",
                        "1",
                        "ref",
                        "qry",
                    ]
                )
            ),
            snp.Snp(
                "\t".join(
                    [
                        "43",
                        ".",
                        "A",
                        "100",
                        "x",
                        "x",
                        "300",
                        "400",
                        "x",
                        "1",
                        "ref",
                        "qry",
                    ]
                )
            ),
            snp.Snp(
                "\t".join(
                    [
                        "43",
                        ".",
                        "A",
                        "100",
                        "x",
                        "x",
                        "300",
                        "400",
                        "x",
                        "1",
                        "ref2",
                        "qry",
                    ]
                )
            ),
            snp.Snp(
                "\t".join(
                    [
                        "43",
                        ".",
                        "A",
                        "100",
                        "x",
                        "x",
                        "300",
                        "400",
                        "x",
                        "1",
                        "ref",
                        "qry2",
                    ]
                )
            ),
            snp.Snp(
                "\t".join(
                    [
                        "44",
                        ".",
                        "A",
                        "100",
                        "x",
                        "x",
                        "300",
                        "400",
                        "x",
                        "1",
                        "ref",
                        "qry",
                    ]
                )
            ),
            snp.Snp(
                "\t".join(
                    [
                        "42",
                        ".",
                        "A",
                        "100",
                        "x",
                        "x",
                        "300",
                        "400",
                        "x",
                        "1",
                        "ref",
                        "qry",
                    ]
                )
            ),
            snp.Snp(
                "\t".join(
                    [
                        "42",
                        "A",
                        ".",
                        "100",
                        "x",
                        "x",
                        "300",
                        "400",
                        "x",
                        "1",
                        "ref",
                        "qry",
                    ]
                )
            ),
        ]

        assert len(initial_vars) == len(to_add)

        for i in range(len(initial_vars)):
            var = variant.Variant(initial_vars[i])
            var_original = copy.copy(var)
            self.assertFalse(var.update_indel(to_add[i]))
            self.assertEqual(var, var_original)

    def test_update_indel_insertion(self):
        """Test update_indel extends insertions correctly"""
        insertion = variant.Variant(
            snp.Snp(
                "\t".join(
                    [
                        "42",
                        ".",
                        "A",
                        "100",
                        "x",
                        "x",
                        "300",
                        "400",
                        "x",
                        "-1",
                        "ref",
                        "qry",
                    ]
                )
            )
        )
        to_add = snp.Snp(
            "\t".join(
                ["42", ".", "C", "101", "x", "x", "300", "400", "x", "-1", "ref", "qry"]
            )
        )
        expected = copy.copy(insertion)
        # coords stored zero-based, so subtract 1 from the real expected coords
        expected.ref_start = 41
        expected.ref_end = 41
        expected.ref_length = 300
        expected.ref_name = "ref"
        expected.ref_base = "."
        expected.qry_start = 99
        expected.qry_end = 100
        expected.qry_length = 400
        expected.qry_name = "qry"
        expected.qry_base = "AC"
        self.assertTrue(insertion.update_indel(to_add))
        self.assertEqual(expected, insertion)

    def test_update_indel_deletion(self):
        """Test update_indel extends deletions correctly"""
        deletion = variant.Variant(
            snp.Snp(
                "\t".join(
                    [
                        "42",
                        "A",
                        ".",
                        "100",
                        "x",
                        "x",
                        "300",
                        "400",
                        "x",
                        "1",
                        "ref",
                        "qry",
                    ]
                )
            )
        )
        to_add = snp.Snp(
            "\t".join(
                ["43", "C", ".", "100", "x", "x", "300", "400", "x", "1", "ref", "qry"]
            )
        )
        expected = copy.copy(deletion)
        # coords stored zero-based, so subtract 1 from the real expected coords
        expected.ref_start = 41
        expected.ref_end = 42
        expected.ref_length = 300
        expected.ref_name = "ref"
        expected.ref_base = "AC"
        expected.qry_start = 99
        expected.qry_end = 99
        expected.qry_length = 400
        expected.qry_name = "qry"
        expected.qry_base = "."
        self.assertTrue(deletion.update_indel(to_add))
        self.assertEqual(expected, deletion)
