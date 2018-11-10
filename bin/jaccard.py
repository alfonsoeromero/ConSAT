#!/usr/bin/env python
"""Binary for the 'jaccard' app"""

import sys
from gfam.scripts.jaccard import JaccardSimilarityApp

if __name__ == "__main__":
    sys.exit(JaccardSimilarityApp().run())
