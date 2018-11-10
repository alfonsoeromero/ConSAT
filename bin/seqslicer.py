#!/usr/bin/env python
"""Binary for the 'seqslicer' app"""

import sys
from gfam.scripts.seqslicer import SeqSlicerApp

if __name__ == "__main__":
    sys.exit(SeqSlicerApp().run())
