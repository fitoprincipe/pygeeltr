#!/usr/bin/env python
# -*- coding: utf-8 -*-
from pygeeltr.tests import landtrendr
import unittest
import sys

if __name__ == "__main__":
    if len(sys.argv) > 3:
        landtrendr.LandTrendr.FOLDER = sys.argv.pop()
        landtrendr.LandTrendr.USER = sys.argv.pop()

    unittest.main(landtrendr)