# -*- coding: utf-8 -*-

""" Package to apply LandTrendr in Google Earth Engine """

from __future__ import absolute_import, division, print_function
from ._version import __version__
from . import landtrendr, classification, ipytools
from .landtrendr import LandTrendr

__all__ = (
    "__title__", "__summary__", "__uri__", "__version__", "__author__",
    "__email__", "__license__", "__copyright__",
)

__title__ = "gee-landtrendr"
__summary__ = "Apply LandTrendr Algorithm in Google Earth Engine time series"
__uri__ = "https://github.com/fitoprincipe/pygeeltr"

__author__ = "Rodrigo E. Principe"
__email__ = "rprincipe@ciefap.org.ar"

__license__ = "GNU GENERAL PUBLIC LICENSE, Version 3"
__copyright__ = "Rodrigo E. Principe"