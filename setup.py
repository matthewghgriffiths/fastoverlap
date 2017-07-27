# -*- coding: utf-8 -*-



from numpy.distutils.core import Extension
from numpy.distutils.system_info import get_info

fast_ext = Extension(name = 'fastbulk', 
                     sources=['fastoverlap/fastbulk.f90', 'fastoverlap/minperm.f90'],
                     **get_info('fftw'))

if __name__ == "__main__":
    from numpy.distutils.core import setup
    setup(name='fastoverlap',
          version='0.1',
          description='Algorithms for fast alignment of atomic structures in finite and periodic systems',
          url='https://github.com/matthewghgriffiths/fastoverlap',
          author='Matthew Griffiths',
          author_email='matthewghgriffiths@gmail.com',
          license='GNU General Public License',
          packages=['fastoverlap'],
          #install_requires=['numpy', 'scipy', 'hungarian'],
          #zip_safe=False,
          ext_modules = [fast_ext])