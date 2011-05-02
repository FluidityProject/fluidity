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

class dependency_list(object):
    def __init__(self, obj, source, dep_strings):
        self.obj=obj
        self.source=source
        
        self.targets=set()
        self.deps=set()

        intargets=True
        for dep in dep_strings:
            for d in dep.split():
                # Drop continuaton characters.
                if d=="\\":
                    continue
                
                if intargets:
                    if d[-1]==":":
                        self.targets.add(d[:-1])
                        intargets=False
                    else:
                        self.targets.add(d)

                else:
                    self.deps.add(d)
        # Treat the .o and the .F90 specially.
        self.targets.remove(obj)
        self.deps.remove(source)
        # Gfortran produces spurious circular dependencies if the .F90
        # contains both a module and routines which use that module.
        self.deps.difference_update(self.targets)
        
    def remove_dep_by_rule(self, f):
        
        discards=set()
        for dep in self.deps:
            if f(dep):
                discards.add(dep)
        self.deps.difference_update(discards)
    
    def as_strings(self):
        
        out=[]

        # Special rule to fake .mod dependency on .o.
        for t in self.targets:
            if t.endswith(".mod"):
                out+=wrap(t+": "+self.obj)+["\n"]
                out+=["\t@true\n","\n"]

        # Main rule.
        out+=wrap(self.obj+" "
                  +" ".join(self.targets)+": "
                  +self.source+" "
                  +" ".join(self.deps))+["\n","\n"]
        
        return out

                           
def wrap(string):
    """Wrap dependencies string according to makefile conventions"""
    linelen=78
    
    lines=[]
    line=""

    for substring in string.split():
        # Always put one thing in a line: 
        if len(line) == 0:
            line = substring
        elif len(line) + len(substring)+1 <= linelen:
            line += " " + substring
        else:
            lines.append(line+" \\\n")
            line="   " + substring
    # Put the last line in.
    lines.append(line)

    return lines

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
            obj=os.path.splitext(f)[0]+".o"
            os.system("rm "+obj+" 2>\\dev\\null || true")

            pipe=os.popen('make GENFLAGS="-M -MF '+obj+'_dependencies" '+obj)
            stdout=pipe.readlines()
            if pipe.close()==None:
                
                this_deps=dependency_list(
                    obj,
                    f,
                    file(obj+"_dependencies","r").readlines())
                # Get rid of absolute paths.
                # We are only interested in dependencies from within the
                # Fluidity tree  so we dump the ones from outside.  
                this_deps.remove_dep_by_rule(lambda x: x.startswith("/"))
                # Remove dependencies on ../confdefs.h because it causes
                # lots of spurious rebuilds.
                this_deps.remove_dep_by_rule(lambda x:
                    x.startswith("../confdefs.h"))

                dependencies+=this_deps.as_strings()
                #split_module_dependency(
                #    strip_absolute_paths(this_deps))+["\n"]
                discards.add(f)
            os.system("rm "+obj+"_dependencies 2>\\dev\\null || true")
        fortran.difference_update(discards)
            

    if len(fortran)>0:
        print "Failed to generate all dependencies. Failures:"
        print str(fortran)
        raise OSError

    return dependencies

def handle_options():
    from optparse import OptionParser
    optparser=OptionParser(usage='usage: %prog [options] <filename>',
                           add_help_option=True,
                           description="""Use gfortran (>=4.5) to automatically """ + 
                           """generate makefile module dependencies.""")

    optparser.add_option("--exclude",
                  help="list of .F90 files to exclude from consideration.",
                  action="store", type="string", dest="exclude", default="")

    (options, argv) = optparser.parse_args()

    return options

if __name__=='__main__':
    options=handle_options()

    sys.stderr.write("Making clean\n")
    trysystem("make clean")

    sys.stderr.write("Listing F90 files\n")
    fortran=set(glob.glob("*.F90")).difference(exceptions)\
        .difference(set(options.exclude.split()))

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
