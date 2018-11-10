#!/usr/bin/env python
"""Binary for the 'blast_filter' app"""

import sys
from gfam.scripts.blast_filter import BlastFilterApp

if __name__ == "__main__":
    sys.exit(BlastFilterApp().run())
