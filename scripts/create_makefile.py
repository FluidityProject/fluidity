#!/usr/bin/python
import re
import os
import sys
import glob

exceptions=set(["Refcount_interface_templates.F90",
                "Refcount_templates.F90",
                "testshapefunctions.F90",
                "test_element_numbering.F90",
                "Residual_estimation.F90",
                "mmpde.F90"])

def trysystem(command):
    if not(os.system(command)==0):
        raise OSError

def create_refcounts():

    refcounts_raw=os.popen("grep include.*Ref *.F90").read()

    refcounts=re.findall(r'"Reference_count.*?"',refcounts_raw)

    if len(refcounts)>0:
        trysystem("make "+" ".join(refcounts))

def strip_makefile(filename):
    
    input=file(filename,'r').readlines()

    for i,l in enumerate(input):
        if l.startswith("#####Automatically generated dependencies. Do not edit below this line#####"):
            break

    file(filename,'w').writelines(input[:i+1])

def rebuild_makefile(filename, dependencies):
    
    outfile=file(filename,'a')

    outfile.writelines(dependencies)

def strip_absolute_paths(dependencies):
    # We are only interested in dependencies from within the Fluidity tree
    # so we dump the ones from outside. 
    new_dep=[]
    firstrow=True
    
    for dep in dependencies:
        tmp=[]
        for d in dep.split():
            if not(d.startswith("/")):
                tmp.append(d)

        if (len(tmp)>0 and tmp[0]!="\\"):
            if firstrow:
                firstrow=False
                tmp=["\n"]+tmp
            else:
                tmp=["   "]+tmp
                
            # [1:] strips the spurious leading space.
            new_dep.append((tmp[0]+" ".join(tmp[1:]))+"\n")

    
    return new_dep

def split_module_dependency(dependencies):
    # Take a rule foo.o bar.mod: blah blah
    # and produce an extra rule:
    # bar.mod: foo.o
    #       @true
    targets=[]

    try:
        for dep in dependencies:
            for d in dep.split():
                if d!="\\":
                    targets.append(d)
                if len(targets)==2:
                    raise StopIteration
    except StopIteration:
        if targets[1].endswith(".mod:"):
            return [targets[1]+" "+targets[0]+"\n",
                    "	@true\n"]+dependencies
    
    return dependencies
        

def generate_dependencies(fortran):
    import os.path
    
    setsize=len(fortran)+1 # Make sure loop executes once.
    
    dependencies=[]
    # Loop as long as we are making progress.
    while len(fortran)<setsize and len(fortran)>0:
        print("Generating dependencies, "+`len(fortran)`+" to go.")
        setsize=len(fortran)        

        discards=set()
        for f in fortran:            
            print "  "+f
            object=os.path.splitext(f)[0]+".o"
            os.system("rm "+object+" 2>\\dev\\null || true")
            pipe=os.popen("make GENFLAGS=-M "+object)
            stdout=pipe.readlines()
            if pipe.close()==None:
                # The first two lines are the command line.
                dependencies+=split_module_dependency(
                    strip_absolute_paths(stdout[2:]))+["\n"]
                discards.add(f)

        fortran.difference_update(discards)
            

    if len(fortran)>0:
        print "Failed to generate all dependencies. Failures:"
        print str(fortran)
        raise OSError

    return dependencies

if __name__=='__main__':
    sys.stderr.write("Making clean\n")
    trysystem("make clean")

    sys.stderr.write("Listing F90 files\n")
    fortran=set(glob.glob("*.F90")).difference(exceptions)

    sys.stderr.write("Creating reference counts\n")
    create_refcounts()

    sys.stderr.write("Stripping Makefile\n")
    strip_makefile("Makefile")

    dependencies=generate_dependencies(fortran)

    sys.stderr.write("Inserting dependencies in Makefile\n")
    rebuild_makefile("Makefile", dependencies)

    sys.stderr.write("Stripping Makefile.in\n")
    strip_makefile("Makefile.in")

    sys.stderr.write("Inserting dependencies in Makefile.in\n")
    rebuild_makefile("Makefile.in", dependencies)
