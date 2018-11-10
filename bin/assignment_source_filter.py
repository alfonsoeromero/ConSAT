#!/usr/bin/env python
"""Binary for the 'assignment_source_filter' app"""

import sys
from gfam.scripts.assignment_source_filter import AssignmentSourceFilterApp


if __name__ == "__main__":
    sys.exit(AssignmentSourceFilterApp().run())
