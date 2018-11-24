# -*- coding: utf-8 -*-
import setuptools
from numpy.distutils.core import setup
from numpy.distutils.core import Extension
from numpy.distutils.system_info import get_info
from numpy.distutils.ccompiler import new_compiler
from numpy.distutils.fcompiler import new_fcompiler


lapack = get_info('lapack')
fftw = get_info('fftw')
queue = dict(
    libraries=['queue'], library_dirs=['build'], include_dirs=['build'])


def merge_libs(*libs):
    lib = {}
    for k in set(k for l in libs for k in l):
        for l in libs:
            v = l.get(k)
            if isinstance(v, list):
                lib.setdefault(k, []).extend(v)
    return lib


extra_compile_args=[
    '-Wall', '-Wextra', '-pedantic', '-funroll-loops', '-O2']
fastbulk_ext = Extension(name='fastoverlap.f90.fastbulk',
                         sources=['fastoverlap/f90/fastbulk.f90'],
                         extra_compile_args=extra_compile_args,
                         **merge_libs(fftw, lapack))
fastcluster_ext = Extension(name='fastoverlap.f90.fastclusters',
                            sources=['fastoverlap/f90/fastclusters.f90'],
                            extra_compile_args=extra_compile_args,
                            **merge_libs(fftw, lapack))
bnb_ext = Extension(name='fastoverlap.f90.libbnb',
                    sources=['fastoverlap/f90/bnbalign.f90'],
                    extra_compile_args=extra_compile_args,
                    **merge_libs(lapack, queue))


if __name__ == "__main__":
    # compiling static fortran priority queue library
    fcompiler = new_fcompiler(compiler='gfortran')
    ccompiler = new_compiler()
    fcompiler.customize()
    queue_objs = fcompiler.compile(
        sources=["fastoverlap/f90/priorityqueue.f90"],
        output_dir='build',
        extra_preargs=['-c', '-fPIC'])
    ccompiler.create_static_lib(
        queue_objs, "queue",
        output_dir='build', debug=1)

    setup(name                 = 'fastoverlap',
          version              = '0.1',
          description          = (
              'Algorithms for fast alignment of atomic'
              'structures in finite and periodic systems'),
          url                  = 'https://github.com/matthewghgriffiths/fastoverlap',
          author               = 'Matthew Griffiths',
          author_email         ='matthewghgriffiths@gmail.com',
          license              ='GNU General Public License',
          packages             = setuptools.find_packages(),
          package_data         = {'': ['*.f90']},
          include_package_data = True,
          ext_modules=[fastbulk_ext, fastcluster_ext, bnb_ext])
