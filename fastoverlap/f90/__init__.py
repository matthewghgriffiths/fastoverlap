# -*- coding: utf-8 -*-

try:
    import fastbulk
    import fastclusters
    import libbnb
    have_fortran = True
except ImportError:
    have_fortran = False