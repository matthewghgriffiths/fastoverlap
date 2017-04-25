
cd build
gfortran -c -fPIC ../branchnbound/kdtree2.f90
ar crs libkdtree.a kdtree2.o

gfortran -c -fPIC ../branchnbound/priorityqueue.f90
ar crs libqueue.a priorityqueue.o

cd ../branchnbound
f2py -c bnbalign.f90 --fcompiler=gfortran -L../build -I../build -lkdtree -lqueue -m libbnb --link-lapack

cd ../fastoverlap
f2py -c fastbulk.f90 -m fastbulk --link-lapack --link-fftw
f2py -c fastcluster.f90 -m fastcluster --link-lapack --link-fftw
