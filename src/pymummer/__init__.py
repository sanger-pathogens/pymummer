import shutil
import sys

from pymummer._version import __version__

required_progs = ["nucmer", "show-coords", "show-snps", "delta-filter"]
missing_progs = [prog for prog in required_progs if shutil.which(prog) is None]

if missing_progs:
    print(
        "The following required programs from the MUMmer package are not found in your PATH:",
        file=sys.stderr,
    )
    for prog in missing_progs:
        print(f"   NOT FOUND: {prog}", file=sys.stderr)
    raise ImportError(
        "Some required MUMmer programs are missing. Please install them and ensure they are in your PATH."
    )

__all__ = [
    "alignment",
    "coords_file",
    "nucmer",
    "snp",
    "snp_file",
    "syscall",
    "variant",
]

from pymummer import *
