from setuptools import setup, find_packages
import os

versionfile=os.popen("svn info |grep Revision|cut -d: -f2", "r")

version="3.4.svn"+versionfile.read().strip()

versionfile.close()

setup(
    name = "fluidity",
    version = version,
    py_modules = ["fluidity_tools","vtktools", "tvtktools"],
)
