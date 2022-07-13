import glob
import os.path
from distutils.core import setup

try:
    destdir = os.environ["DESTDIR"]
except KeyError:
    destdir = ""

setup(
    name="dxdiff",
    version="1.0",
    description="An XML aware diff tool.",
    author="The ICOM team",
    author_email="fraser.waters08@imperial.ac.uk",
    url="http://amcg.ese.ic.ac.uk",
    packages=["dxdiff"],
    scripts=["dxdiff/dxdiff"],
)
