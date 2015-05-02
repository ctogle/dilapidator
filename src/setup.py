from distutils.core import setup
from Cython.Build import cythonize

import os,numpy,appdirs

from distutils.extension import Extension
from Cython.Distutils import build_ext

#from setuptools import setup,Extension

core_modules = []
ext_modules = [
    #Extension('dp_utils', ['support/mp_utils.c']), 
    Extension('dp_vector',['cython/dp_vector.c']), 
    #Extension('dp_bboxes', ['support/mp_bboxes.c']), 
    #Extension('dp_terrain', ['support/mp_terrain.c']), 
            ]

resourcesdir = os.path.join(appdirs.user_data_dir(),'dilap_resources')
resourcesrcd = os.path.join(os.getcwd(),'resources')
resourcefils = []

for rpath in os.walk(resourcesrcd):
    rsrcp = rpath[0][len(os.getcwd())+1:]
    for rfile in rpath[2]:
        rfi = '/'.join([rsrcp,rfile])
        resourcefils.append(rfi)

pkgs = [
    'dilap',
    'dilap.io',
    'dilap.core',
    'dilap.primitive',
    'dilap.generate',
]

setup(
    name="dilapidator",
    version = '1.0',
    description = "dilapidator python pkg",
    author = "ctogle",
    author_email = "cogle@vt.edu",
    license = "MIT License",
    long_description = 'procedural model construction/dilapidation', 
    packages = pkgs, 
    py_modules = core_modules, 
    ext_modules = ext_modules, 
    include_dirs = [numpy.get_include()], 
    data_files=[(resourcesdir,resourcefils)], 
    )




