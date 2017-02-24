# FASTOVERLAP
Algorithms for fast alignment of structures in finite and periodic systems.

## INSTALLATION

## Required packages

### for compilation:

1. fortran compiler
2. fftw

### python packages:

1. numpy:
  We use numpy everywhere for doing numerical work. It also installs f2py which is used to compile fortran code into modules callable by python.

2. scipy:
  For some of the optimizers and various scientific tools

3. hungarian:
  For permutational alignment

## Compilation

To use fortran extenstions using numpy.distutils

```
$ python setup.py build --fcompiler=gfortran
$ python setup.py install [--user]
```

Building in place (note PYTHONPATH must then be modified)

```
$ python setup.py build_ext -i --fcompiler=gfortran
```
