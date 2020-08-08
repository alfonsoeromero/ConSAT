import bz2
import gzip
import sys
from io import IOBase
from urllib.request import urlopen


def open_anything(fname, *args, **kwds):
    """Opens the given file. The file may be given as a file object
    or a filename. If the filename ends in ``.bz2`` or ``.gz``, it will
    automatically be decompressed on the fly. If the filename starts
    with ``http://``, ``https://`` or ``ftp://`` and there is no
    other argument given, the remote URL will be opened for reading.
    A single dash in place of the filename means the standard input.
    """
    if isinstance(fname, IOBase):
        infile = fname
    elif fname == "-":
        infile = sys.stdin
    elif (fname.startswith("http://") or fname.startswith("ftp://") or
          fname.startswith("https://")) and not kwds and not args:
        infile = urlopen(fname)
    elif fname.endswith(".bz2"):
        infile = bz2.BZ2File(fname, *args, **kwds)
    elif fname.endswith(".gz"):
        infile = gzip.GzipFile(fname, *args, **kwds)
    else:
        infile = open(fname, *args, **kwds)
    return infile
