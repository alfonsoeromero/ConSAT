#!/usr/bin/env python
"""Binary for the 'blast_all' app"""

import sys
from gfam.scripts.blast_all import AllAgainstAllBLASTApp

if __name__ == "__main__":
    sys.exit(AllAgainstAllBLASTApp().run())
