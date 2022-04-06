from setuptools import Extension, setup

setup(
    ext_modules=[
        Extension(
            name='libspud', sources=['libspud.c'], include_dirs=['../include'],
            library_dirs=['../../build/lib'], libraries=['spud'],
        )
    ]
)
