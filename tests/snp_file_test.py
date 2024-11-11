import unittest
import os
from pymummer import snp_file, snp, variant

data_dir = "tests/data"


class TestUtils(unittest.TestCase):
    def test_snp_file(self):
        """test coords_file"""
        expected = [
            "\t".join(
                [
                    "133",
                    "G",
                    ".",
                    "122",
                    "1",
                    "122",
                    "500",
                    "489",
                    "1",
                    "1",
                    "ref",
                    "qry",
                ]
            ),
            "\t".join(
                [
                    "143",
                    ".",
                    "C",
                    "131",
                    "1",
                    "132",
                    "500",
                    "489",
                    "1",
                    "1",
                    "ref",
                    "qry",
                ]
            ),
            "\t".join(
                [
                    "253",
                    "T",
                    "A",
                    "242",
                    "120",
                    "242",
                    "500",
                    "489",
                    "1",
                    "1",
                    "ref",
                    "qry",
                ]
            ),
        ]

        expected = [snp.Snp(x) for x in expected]

        infiles = [
            os.path.join(data_dir, "snp_file_test_with_header.snps"),
            os.path.join(data_dir, "snp_file_test_no_header.snps"),
        ]

        for fname in infiles:
            fr = snp_file.reader(fname)
            snps = [x for x in fr]
            self.assertEqual(snps, expected)

    def test_get_all_variants(self):
        """Test load all variants from file"""
        deletion_snps = [
            "\t".join(
                [
                    "125",
                    "T",
                    ".",
                    "124",
                    "1",
                    "124",
                    "500",
                    "497",
                    "1",
                    "1",
                    "ref1",
                    "qry1",
                ]
            ),
            "\t".join(
                [
                    "126",
                    "A",
                    ".",
                    "124",
                    "1",
                    "124",
                    "500",
                    "497",
                    "1",
                    "1",
                    "ref1",
                    "qry1",
                ]
            ),
            "\t".join(
                [
                    "127",
                    "C",
                    ".",
                    "124",
                    "1",
                    "124",
                    "500",
                    "497",
                    "1",
                    "1",
                    "ref1",
                    "qry1",
                ]
            ),
        ]
        deletion_snps = [snp.Snp(x) for x in deletion_snps]
        deletion_variant = variant.Variant(deletion_snps[0])
        deletion_variant.update_indel(deletion_snps[1])
        deletion_variant.update_indel(deletion_snps[2])

        just_a_snp = "\t".join(
            [
                "386",
                "C",
                "T",
                "383",
                "115",
                "115",
                "500",
                "497",
                "1",
                "1",
                "ref1",
                "qry1",
            ]
        )
        snp_variant = variant.Variant(snp.Snp(just_a_snp))

        insertion_snps = [
            "\t".join(
                [
                    "479",
                    ".",
                    "G",
                    "480",
                    "0",
                    "22",
                    "500",
                    "504",
                    "1",
                    "1",
                    "ref2",
                    "qry2",
                ]
            ),
            "\t".join(
                [
                    "479",
                    ".",
                    "A",
                    "481",
                    "0",
                    "22",
                    "500",
                    "504",
                    "1",
                    "1",
                    "ref2",
                    "qry2",
                ]
            ),
            "\t".join(
                [
                    "479",
                    ".",
                    "T",
                    "482",
                    "0",
                    "22",
                    "500",
                    "504",
                    "1",
                    "1",
                    "ref2",
                    "qry2",
                ]
            ),
            "\t".join(
                [
                    "479",
                    ".",
                    "A",
                    "483",
                    "0",
                    "22",
                    "500",
                    "504",
                    "1",
                    "1",
                    "ref2",
                    "qry2",
                ]
            ),
        ]
        insertion_snps = [snp.Snp(x) for x in insertion_snps]
        insertion_variant = variant.Variant(insertion_snps[0])
        for i in range(1, len(insertion_snps)):
            insertion_variant.update_indel(insertion_snps[i])

        variants_from_file = snp_file.get_all_variants(
            os.path.join(data_dir, "snp_file_test_get_all_variants.snps")
        )
        self.assertEqual(len(variants_from_file), 3)
        self.assertEqual(variants_from_file[0], deletion_variant)
        self.assertEqual(variants_from_file[1], snp_variant)
        self.assertEqual(variants_from_file[2], insertion_variant)
