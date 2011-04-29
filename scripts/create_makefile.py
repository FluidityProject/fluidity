#!/usr/bin/python
import re
import os
import sys

def trysystem(command):
    if not(os.system(command)==0):
        raise OSError

def create_refcounts():

    refcounts_raw=os.popen("grep include.*Ref *.F90").read()

    refcounts=re.findall(r'"Reference_count.*?"',refcounts_raw)

    trysystem("make "+" ".join(refcounts))

def strip_makefile(filename):
    
    input=file(filename,'r').readlines()

    for i,l in enumerate(input):
        if l.startswith("#####Automatically generated dependencies. Do not edit below this line#####"):
            break

    file(filename,'w').writelines(input[:i+1])

def generate_dependencies():
    import glob

    fortran=set(glob.glob("*.F90"))
    
    setsize=len(s)+1 # Make sure loop executes once.
    
    # Loop as long as we are making progress.
    while len(s)<setsize:
        discards=set()
        for f in fortran:
            os.popen

sys.stderr.write("Making clean\n")
trysystem("make clean")
sys.stderr.write("Creating reference counts\n")
create_refcounts()
sys.stderr.write("Stripping Makefile\n")
strip_makefile("Makefile")
