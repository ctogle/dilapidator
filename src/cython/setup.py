from distutils.core import setup
from Cython.Build import cythonize

import numpy

from distutils.extension import Extension
from Cython.Distutils import build_ext

vectorext = ('dp_vector',[Extension('dp_vector',['dp_vector.pyx','dp_vector.pxd'])])
quatrnext = ('dp_quaternion',[Extension('dp_quaternion',['dp_quaternion.pyx'])])
raycsnext = ('dp_ray',[Extension('dp_ray',['dp_ray.pyx'])])

def ext_setup(ext_name, ext_modules):
    setup(name = ext_name, 
      	ext_modules = ext_modules, 
      	cmdclass = {'build_ext': build_ext}, 
      	include_dirs = [numpy.get_include(),'.'])

ext_setup(*vectorext)
ext_setup(*quatrnext)
ext_setup(*raycsnext)

#python setup.py build_ext --inplace


