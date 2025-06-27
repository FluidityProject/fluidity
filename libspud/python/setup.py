from setuptools import Extension, setup

setup(
    ext_modules=[
        Extension(
            name="libspud",
            sources=["libspud.c"],
            include_dirs=["/home/fluidity/skramer/git/fluidity/libspud/include"],
            library_dirs=["/home/fluidity/skramer/git/fluidity/build/lib"],
            libraries=["spud"],
            extra_link_args=["-Wl,--enable-new-dtags,-R/home/fluidity/skramer/git/fluidity/build/lib"],
        )
    ],
)
