#!/usr/bin/env python3
from distutils.core import setup
from distutils.extension import Extension
from Cython.Build import cythonize
from Cython.Distutils import build_ext
import appdirs
import shutil
import numpy
import os


userdata = os.path.join(appdirs.user_data_dir(),'dilap_resources')
rsrcdata = os.path.join(os.getcwd(),'resources')
shutil.rmtree(userdata)
shutil.copytree(rsrcdata,userdata)


scripts = [
    'bin/dilap-gen',
    'bin/dilap-serve', ]


pkgs = [
    'dilap.geometry',
    'dilap.topology',
    'dilap.modeling',
    'dilap.worldly',
    'dilap',
    'dilap.io',
    'dilap.core', ]


extnames = [
    'dilap.geometry.vec3',
    'dilap.geometry.quat',
    'dilap.geometry.ray3',
    'dilap.geometry.curve',
    'dilap.geometry.tform',
    'dilap.geometry.pointset',
    'dilap.geometry.tools',
    'dilap.geometry.triangulate', ]
ext = lambda name : Extension(name, [name.replace('.','/')+'.pyx', name.replace('.','/')+'.pxd'])
exts = [ext(name) for name in extnames]


def install(*ags):
    setup(
        script_args=ags,
        name='dilapidator',
        version='1.0',
        description='procedural mesh generation',
        author='ctogle',
        author_email='curtis.t.ogle@gmail.com',
        url='http://github.com/ctogle/dilapidator.git',
        license='MIT License',
        long_description='procedural mesh generation for export to other applications (.obj, .fbx)',
        scripts=scripts,
        packages=pkgs, 
        ext_modules=exts,
        cmdclass={'build_ext': build_ext},
        include_dirs=[numpy.get_include()],
        data_files=[], )


if __name__ == '__main__':
    install('build_ext','install','--user')
