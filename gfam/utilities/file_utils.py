import os
import platform
import sys
from contextlib import contextmanager
from shutil import rmtree
from tempfile import mkdtemp

__author__ = "Tamas Nepusz, Alfonso E. Romero"
__email__ = "tamas@cs.rhul.ac.uk, aeromero@cs.rhul.ac.uk"
__copyright__ = "Copyright (c) 2014, Tamas Nepusz"
__license__ = "GPL"

__all__ = ["redirected",
           "search_file", "temporary_dir"]


# pylint:disable-msg=W0613
# W0613: unused argument.
@contextmanager
def redirected(stdin=None, stdout=None, stderr=None):
    """Context manager that temporarily redirects some of the input/output
    streams.

    `stdin` is the new standard input, `stdout` is the new standard output,
    `stderr` is the new standard error. ``None`` means to leave the
    corresponding stream unchanged.

    Example::

        with redirected(stdout=open("test.txt", "w")):
            print "Yay, redirected to a file!"
    """
    loc = locals()
    stream_names = ["stdin", "stdout", "stderr"]
    old_streams = {}
    try:
        for sname in stream_names:
            stream = loc.get(sname, None)
            if stream is not None and stream != getattr(sys, sname):
                old_streams[sname] = getattr(sys, sname)
                setattr(sys, sname, loc.get(sname, None))
        yield
    finally:
        for sname, stream in old_streams.items():
            setattr(sys, sname, stream)


def search_file(filename, search_path=None, executable=True):
    """Finds the given `filename` in the given search path.
    If `executable` and we are on Windows, ``.exe`` will be
    appended to the filename. Returns the full path of the
    file or ``None`` if it is not found on the path."""
    if not search_path:
        search_path = os.environ["PATH"]

    if executable and platform.system() == "Windows":
        filename = filename+".exe"

    exists = os.path.exists
    join = os.path.join

    for path in search_path.split(os.pathsep):
        fullpath = join(path, filename)
        if exists(fullpath):
            return os.path.abspath(fullpath)

    return None


@contextmanager
def temporary_dir(*args, **kwds):
    """Context manager that creates a temporary directory when entering the
    context and removes it when exiting.

    Every argument is passed on to `tempfile.mkdtemp` except a keyword argument
    named `change` which tells whether we should change to the newly created
    temporary directory or not. The current directory will be restored when
    exiting the context manager."""
    change = "change" in kwds
    if change:
        del kwds["change"]

    name = mkdtemp(*args, **kwds)
    try:
        if change:
            old_dir = os.getcwd()
            os.chdir(name)
        yield name
    finally:
        if change:
            os.chdir(old_dir)
        rmtree(name)
