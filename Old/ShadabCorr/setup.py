# python f-setup.py build_ext --inplace
#   cython f.pyx -> f.cpp
#   g++ -c f.cpp -> f.o
#   g++ -c fc.cpp -> fc.o
#   link f.o fc.o -> f.so

# distutils uses the Makefile distutils.sysconfig.get_makefile_filename()
# for compiling and linking: a sea of options.

# http://docs.python.org/distutils/introduction.html
# http://docs.python.org/distutils/apiref.html  20 pages ...
# http://stackoverflow.com/questions/tagged/distutils+python

import numpy
from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext
# from Cython.Build import cythonize

ext_modules = [Extension(
    name="Scorr",
    sources=["calc.pyx", "corr2dc.c"],
    # extra_objects=["fc.o"],  # if you compile fc.cpp separately
    include_dirs = [numpy.get_include()],  # .../site-packages/numpy/core/include
    language="c",
    # libraries=
    extra_compile_args = ['-O3'],
    # extra_link_args = "...".)
    )]

setup(
    name = 'Scorr',
    cmdclass = {'build_ext': build_ext},
    ext_modules = ext_modules,)
