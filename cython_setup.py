from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext

ext_modules =[
    Extension('csie',
              ['csie.pyx'])
]

setup(
    name= 'GravLens Model Modules',
    cmdclass= {'build_ext':build_ext},
    ext_modules= ext_modules)
