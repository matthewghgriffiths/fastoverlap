# FASTOVERLAP
Algorithms for fast alignment of structures in finite and periodic systems.

## INSTALLATION

## Required packages

### for compilation:

1. fortran compiler
2. fftw
3. lapack

### python packages:

1. numpy:
  We use numpy everywhere for doing numerical work. It also installs f2py which is used to compile fortran code into modules callable by python.

2. scipy:
  For some of the optimizers and various scientific tools

3. hungarian:
  For permutational alignment

## Compilation

To compile f2py modules run

```
$ ./compile.sh
```
