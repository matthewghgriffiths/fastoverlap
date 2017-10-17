
mkdir -p build
cd build
#gfortran -c -fPIC ../fastoverlap/f90/kdtree2.f90
#ar crs libkdtree.a kdtree2.o

gfortran -c -fPIC ../fastoverlap/f90/priorityqueue.f90
ar crs libqueue.a priorityqueue.o

cd ../fastoverlap/f90
f2py -c bnbalign.f90 --fcompiler=gfortran -L../../build -I../../build -lkdtree -lqueue -m libbnb --link-lapack
f2py -c fastbulk.f90 -m fastbulk --link-lapack --link-fftw
f2py -c fastclusters.f90 -m fastclusters --link-lapack --link-fftw
