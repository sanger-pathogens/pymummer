import unittest
import os
from pymummer import syscall


class TestSyscall(unittest.TestCase):
    def test_run_fail(self):
        """Test that run raises error when command fails"""
        with self.assertRaises(syscall.Error):
            syscall.run("notacommandandthrowerror")

    def test_run_ok(self):
        """Test run is ok on command that works"""
        syscall.run("ls")
