from os.path import abspath
from setuptools import setup, Extension

setup(name='libspud', version='1.1.3',
      description='Python bindings for libspud',
      ext_modules=[Extension('libspud', sources=['libspud.c'],
                             libraries=["spud"],
                             library_dirs=[abspath("..")],
                             include_dirs=[abspath("../include")])])
