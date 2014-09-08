import os
import sys

# Ensure our scripts are imported before any others since we can't depend upon
# the packaging middleware to handle relative imports, or should we vendorize
# pip or easy_install?
sys.path = [(os.path.dirname(os.path.abspath(__file__)))] + sys.path

import entry_point

if __name__ == '__main__':
    entry_point.main()
