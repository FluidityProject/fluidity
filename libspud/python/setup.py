import os.path

from distutils.core import Extension
from distutils.core import setup

module1 = Extension(
    "libspud",
    sources=["libspud.c"],
    libraries=["spud"],
    library_dirs=[os.path.abspath("..")],
    include_dirs=[os.path.abspath("../include")],
)

setup(
    name="libspud",
    version="1.1.3",
    description="Python bindings for libspud",
    ext_modules=[module1],
)
