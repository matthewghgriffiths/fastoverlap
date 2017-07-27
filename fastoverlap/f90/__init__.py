# -*- coding: utf-8 -*-

try:
    import fastbulk
    have_fastbulk = True
except ImportError:
    have_fastbulk = False

try:
    import fastclusters
    have_fastclusters = True
except ImportError:
    have_fastclusters = False

try:
    import libbnb
    have_libbnb = True
except ImportError:
    have_libbnb = False

have_fortran = all((have_fastbulk,have_fastclusters,have_libbnb))
