# -*- coding: utf-8 -*-

from setuptools import setup

setup(name='fastoverlap',
      version='0.1',
      description='Algorithms for fast alignment of atomic structures in finite and periodic systems',
      url='https://github.com/matthewghgriffiths/fastoverlap',
      author='Matthew Griffiths',
      author_email='matthewghgriffiths@gmail.com',
      license='GNU General Public License',
      packages=['fastoverlap'],
      install_requires=['numpy', 'scipy', 'hungarian'],
      zip_safe=False)