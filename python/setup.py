from distutils.core import setup
from distutils.extension import Extension
import os
import os.path

try:
  destdir = os.environ["DESTDIR"]
except KeyError:
  destdir = ""

setup(
      name='fluidity',
      version='0.1',
      description="Fluidity python files",
      author = "The ICOM team",
      author_email = "patrick.farrell06@imperial.ac.uk",
      url = "http://amcg.ese.ic.ac.uk",
      packages = ['fluidity', 'fluidity.diagnostics'],
      package_dir = {'fluidity': 'fluidity'},
      py_modules = ['fluidity_tools', 'vtktools']
     )


