import os
import sys

# Ensure our scripts are imported before any others since we can't depend upon
# the packaging middleware to handle relative imports, or should we vendorize
# pip or easy_install?
sys.path = [(os.path.dirname(os.path.abspath(__file__)))] + sys.path


entry_point = __import__(
    "glycresoft_ms2_classification.entry_point").entry_point


def main():
    entry_point.main()

if __name__ == '__main__':
    main()
