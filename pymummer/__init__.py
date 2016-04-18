from pkg_resources import get_distribution

try:
    __version__ = get_distribution('pymummer').version
except:
    __version__ = 'local'


__all__ = [
    'alignment',
    'coords_file',
    'nucmer',
    'snp',
    'snp_file',
    'syscall',
    'variant',
]

from pymummer import *
