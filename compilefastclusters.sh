cd ./fastoverlap/f90
f2py -c fastclusters.f90 -m fastclusters --link-lapack --link-fftw
