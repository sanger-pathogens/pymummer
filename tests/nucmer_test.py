import unittest
import os
import filecmp
from pymummer import nucmer

data_dir = "tests/data"


class TestRunner(unittest.TestCase):
    def test_nucmer_command(self):
        """test _nucmer_command"""
        tests = [
            [nucmer.Runner("ref", "qry", "outfile"), "nucmer -p pre ref qry"],
            [
                nucmer.Runner("ref", "qry", "outfile", breaklen=42),
                "nucmer -p pre -b 42 ref qry",
            ],
            [
                nucmer.Runner("ref", "qry", "outfile", diagdiff=11),
                "nucmer -p pre -D 11 ref qry",
            ],
            [
                nucmer.Runner("ref", "qry", "outfile", diagdiff=11, promer=True),
                "promer -p pre ref qry",
            ],
            [
                nucmer.Runner("ref", "qry", "outfile", maxmatch=True),
                "nucmer -p pre --maxmatch ref qry",
            ],
            [
                nucmer.Runner("ref", "qry", "outfile", mincluster=42),
                "nucmer -p pre -c 42 ref qry",
            ],
            [
                nucmer.Runner("ref", "qry", "outfile", simplify=False),
                "nucmer -p pre --nosimplify ref qry",
            ],
            [
                nucmer.Runner("ref", "qry", "outfile", promer=True),
                "promer -p pre ref qry",
            ],
            [
                nucmer.Runner("ref", "qry", "outfile", promer=True, breaklen=42),
                "promer -p pre -b 42 ref qry",
            ],
            [
                nucmer.Runner("ref", "qry", "outfile", promer=True, maxmatch=True),
                "promer -p pre --maxmatch ref qry",
            ],
            [
                nucmer.Runner("ref", "qry", "outfile", promer=True, simplify=False),
                "promer -p pre ref qry",
            ],
        ]

        for l in tests:
            self.assertEqual(l[0]._nucmer_command("ref", "qry", "pre"), l[1])

    def test_delta_filter_command(self):
        """test _delta_filter_command"""
        tests = [
            [nucmer.Runner("ref", "qry", "outfile"), "delta-filter infile > outfile"],
            [
                nucmer.Runner("ref", "qry", "outfile", min_id=42),
                "delta-filter -i 42 infile > outfile",
            ],
            [
                nucmer.Runner("ref", "qry", "outfile", min_length=43),
                "delta-filter -l 43 infile > outfile",
            ],
        ]

        for l in tests:
            self.assertEqual(l[0]._delta_filter_command("infile", "outfile"), l[1])

    def test_show_coords_command(self):
        """test _show_coords_command"""
        tests = [
            [
                nucmer.Runner("ref", "qry", "outfile", coords_header=False),
                "show-coords -dTlro -H infile > outfile",
            ],
            [
                nucmer.Runner("ref", "qry", "outfile"),
                "show-coords -dTlro infile > outfile",
            ],
        ]

        for l in tests:
            self.assertEqual(l[0]._show_coords_command("infile", "outfile"), l[1])

    def test_show_snps_command(self):
        """test _show_snps_command"""
        tests = [
            [
                nucmer.Runner("ref", "qry", "outfile", snps_header=False),
                "show-snps -TClr -H infile > outfile",
            ],
            [
                nucmer.Runner("ref", "qry", "outfile"),
                "show-snps -TClr infile > outfile",
            ],
            [
                nucmer.Runner("ref", "qry", "outfile", show_snps_C=False),
                "show-snps -Tlr infile > outfile",
            ],
        ]

        for nuc_obj, expected in tests:
            self.assertEqual(nuc_obj._show_snps_command("infile", "outfile"), expected)

    def test_write_script_no_snps(self):
        """test _write_script no snps"""
        tmp_script = "tmp.script.sh"
        r = nucmer.Runner("ref", "qry", "outfile")
        r._write_script(tmp_script, "ref", "qry", "outfile")
        expected = os.path.join(data_dir, "nucmer_test_write_script_no_snps.sh")
        self.assertTrue(filecmp.cmp(expected, tmp_script, shallow=False))
        os.unlink(tmp_script)

    def test_write_script_with_snps(self):
        """test _write_script with snps"""
        tmp_script = "tmp.script.sh"
        r = nucmer.Runner("ref", "qry", "outfile", show_snps="outfile.snps")
        r._write_script(tmp_script, "ref", "qry", "outfile")
        expected = os.path.join(data_dir, "nucmer_test_write_script_with_snps.sh")
        self.assertTrue(filecmp.cmp(expected, tmp_script, shallow=False))
        os.unlink(tmp_script)

    def test_run_nucmer(self):
        """test run_nucmer"""
        qry = os.path.join(data_dir, "nucmer_test_qry.fa")
        ref = os.path.join(data_dir, "nucmer_test_ref.fa")
        tmp_out = "tmp.nucmer.out"
        runner = nucmer.Runner(
            ref, qry, tmp_out, coords_header=False, show_snps=True, snps_header=False
        )
        runner.run()
        expected = os.path.join(data_dir, "nucmer_test_out.coords")
        self.assertTrue(filecmp.cmp(tmp_out, expected, shallow=False))
        self.assertTrue(
            filecmp.cmp(tmp_out + ".snps", expected + ".snps", shallow=False)
        )
        os.unlink(tmp_out)
        os.unlink(tmp_out + ".snps")
