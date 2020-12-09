# FASTOVERLAP
Algorithms for fast alignment of structures in finite and periodic systems. These methods are also implemented in the [Cambridge Energy Landscape Software](http://www-wales.ch.cam.ac.uk/software.html) in `GMIN` and `OPTIM`.

FASTOVERLAP can be run 'as is' with no compilation required using `PeriodicAlign` and `SphericalAlign`. The Fortran modules need to be compiled to peform branch and bound alignments with `BranchnBoundAlignment`. 

These algorithms have been detailed in the paper,

**Optimal Alignment of Structures for Finite and Periodic Systems** Matthew Griffiths, Samuel P. Niblett, and David J. Wales _Journal of Chemical
Theory and Computation_ **2017** 13(10), 4914-4931, doi:[10.1021/acs.jctc.7b00543](http://dx.doi.org/10.1021/acs.jctc.7b00543)

If you use this module, please cite the above paper.

## INSTALLATION

## Required packages

### python packages:

1. `numpy`:
  We use numpy everywhere for doing numerical work. It also installs `f2py` which is used to compile fortran code into modules callable by python.

2. `scipy`:
  For some of the optimizers and various scientific tools

3. [`munkres`](https://pypi.org/project/munkres/) or [`pele`](http://pele-python.github.io/pele/):
  To solve the linear assignment problem for permutational alignment.

### for compilation:

1. Fortran compiler
2. `fftw`
3. `lapack`

## Compilation

The modules can be compiled using `numpy.distutils` to build the Fortran Modules in place

```
$ python setup.py build_ext -i
```

alternative the `f2py` modules can be directly compiled by running

```
$ ./compile.sh
```


