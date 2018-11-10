#!/usr/bin/env python
"""Binary for the 'get_text' app"""

import sys
from gfam.scripts.get_text import GetText

if __name__ == "__main__":
    sys.exit(GetText().run())
