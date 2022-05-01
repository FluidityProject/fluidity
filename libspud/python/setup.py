from os.path import abspath

from setuptools import Extension, setup

setup(
    ext_modules=[
        Extension(
            name="libspud",
            sources=["libspud.c"],
            include_dirs=[abspath("../include")],
            library_dirs=[abspath("../../build/lib")],
            libraries=["spud"],
        )
    ],
)
