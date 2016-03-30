from __future__ import division, absolute_import, print_function
from setuptools import setup
#from numpy.distutils.core import setup


def build_ext(config):
    sief_sources = ['gravpy/sie/fortran/sief.f90']
    config.add_extension(name='sief',sources=sief_sources)
    
def setup_gravpy():

    from numpy.distutils.misc_util import Configuration
    
    config_dict = build_ext(Configuration('gravpy',parent_package=None,top_path=None))
    
    setup(name='gravpy',
          version='0.1',
          description='A general gravitational lens solver writtin in python',
          author='Sourabh Cheedella',
          author_email='cheedella.sourabh@gmail',
          install_requires=['numpy','scipy','numexpr','matplotlib'],
          packages =['gravpy.sie','gravpy.sie.fortran','gravpy.sie.python',
                     'gravpy.alpha','gravpy.alpha.fortran','gravpy.alpha.python',
                     'gravpy.nfw','gravpy.nfw.fortran'],
      #ext_modules= [sie_f]
    )

if __name__ == '__main__':
    setup_gravpy()



