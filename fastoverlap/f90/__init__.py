# -*- coding: utf-8 -*-

try:
    import fastoverlap.f90.fastbulk
    have_fastbulk = True
except ImportError:
    have_fastbulk = False

try:
    import fastoverlap.f90.fastclusters
    have_fastclusters = True
except ImportError:
    have_fastclusters = False

try:
    import fastoverlap.f90.libbnb
    have_libbnb = True
except ImportError:
    have_libbnb = False

have_fortran = all((have_fastbulk,have_fastclusters,have_libbnb))
