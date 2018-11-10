#!/usr/bin/env python
"""Binary for the 'find_unassigned' app"""

import sys
from gfam.scripts.find_unassigned import FindUnassignedApp

if __name__ == "__main__":
    sys.exit(FindUnassignedApp().run())
