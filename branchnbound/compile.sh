
gfortran -c -fPIC kdtree2.f90
ar crs libkdtree.a kdtree2.o

gfortran -c -fPIC priorityqueue.f90
ar crs libqueue.a priorityqueue.o

f2py -c bnbalign.f90 --fcompiler=gfortran -L. -I. -lkdtree -lqueue -m libbnb --link-lapack
