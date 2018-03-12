#!/usr/bin/env python

import re
import glob
import sys
import os
import io
import subprocess

# These files are ignored in the scanning, due to being too weird.
SKIP_LIST=['Diagnostic_Fields_Interfaces.F90','Diagnostic_Fields_New.F90']

def read_fortran_file(filename):
    """ Parse a Fortran 90 file and extract the use statements."""

    data=""

    file_handle=open(filename)
    flag = True
    for line in file_handle.readlines():
        clean_line=" ".join(line.strip().rstrip().split(" "))+'\n'
        if clean_line.lower().startswith("interface"):
            flag = False
        if flag:
            data += clean_line
        if clean_line.lower().startswith("end interface"):
            flag = True
    file_handle.close()

    module_name=re.search('(?<=module )\w+',data,re.I)
    if module_name:
        module_name = module_name.group(0).lower()
        modflag = True
    else:
        modflag = False

    modules_used = re.findall('(?<=^use )\w+',data,re.MULTILINE+re.I)
    modules_used = [module.lower() for module in modules_used]
    module_text = {}

    for module in modules_used:
        text='&'
        test =',.+'
        while text[-1]=='&':
            text = re.search("(?<=^use %s)%s"%(module,test),data,re.MULTILINE)
            if text:
                text = text.group(0)
            else:
                text=""
                break
            test+=r"\n.+"
        module_text[module] = text
    

    return module_name, modules_used, module_text, modflag

def build_trees(filenames):
    """Build the hierarchy of use statements for an input list of filenames"""

    data = dict()
    namelist = dict()

    for filename in filenames:
        filename = os.path.relpath(filename)
    
        name, children, module_text, test = read_fortran_file(filename)
    
        if test:
            namelist[filename]=children, name
        if name:
            data[name]=dict(filename=filename,
                            directly_used_modules=set(children),
                            indirectly_used_modules=set(children),
                            module_text=module_text,
                            scan=True)

    names = tuple(data.keys())
    for name in names:
        add_children(data,name)

    return data, namelist

def print_misplaced(data, namelist, silent, args):
    """Find and list any modules which have use statements out of order."""

    def child_cmp(nm1,nm2):
        """ module comparision for precidence."""
        if nm1 in data.setdefault(nm2,{'indirectly_used_modules':[]})['indirectly_used_modules']:
            return -1
        elif nm2 in data.setdefault(nm1,{'indirectly_used_modules':[]})['indirectly_used_modules']:
            return 1
        else:
            return 0


    def insert_sort(modules, comp):
        """Perform an insertion sort to check that modules are always inserted before any others which use them."""
        out = []

        if modules:
            out.insert(0,modules[0])
        for module in modules[1:]:
            k=len(out)
            for j, val in enumerate(out):
                if comp(module, val) == -1:
                    k=j
                    break
            out.insert(k,module)
            if k!=len(out)-1:
                insert_sort(out,comp)

        return out

  
    order = data.keys()

    out = ""

    high_filenames=[]
    if args:
        high_filenames.extend(args)
    else:
        for nm in order:
            if data[nm]['filename']:
                high_filenames.append(data[nm]['filename'])

    clean_high_filenames = set([os.path.relpath(filename)
                      for filename in high_filenames])



    for filename in clean_high_filenames:
        if (os.path.basename(filename) in SKIP_LIST
            or filename not in namelist):
            continue
        new_order = insert_sort(namelist[filename][0],child_cmp)

        if namelist[filename][1] in data:
            module_text = data[namelist[filename][1]]['module_text']
        else:
            module_text={}

        if new_order != namelist[filename][0]:
            out+=filename+':\n'

            for module in new_order:
                out+= '  use '+module+module_text.setdefault(module,"")+'\n'

    if not silent:
        print(out)
    return out

def add_children(data,name):
    """Add the use statements inherited from the previous generation."""
    for child in data[name]['directly_used_modules']:
        if child in data:
            if data[child]['scan']:
                add_children(data,child)
            data[name]['indirectly_used_modules']=data[name]['indirectly_used_modules'].union(data[child]['indirectly_used_modules'])
        else:
            data[child]=dict(filename=None,
                            directly_used_modules=set(),
                            indirectly_used_modules=set(),
                            module_text="",
                            scan=False)
    data[name]['scan']=False

def make_graphs(tree, namelist):
    """ Output the inheritance tree for named files."""

    from PIL import Image

    high_filenames=[]

    if len(sys.argv)>2:
        high_filenames.extend(sys.argv[2:])

    for filename in high_filenames:
        filename = os.path.relpath(filename)

        outstring = "digraph %s{\n"%namelist[filename][1]
        visited = set()

        def add_used(name, out, visited):

            if name in visited:
                return out

            for module in tree[name]['directly_used_modules']:
                if module not in visited or (name == namelist[filename][1]
                                             and module in namelist[filename][0]):
                    out += "%s->%s\n"%(name,module)
                    out = add_used(module, out, visited)

            visited.add(name)
            return out

        outstring = add_used(namelist[filename][1], outstring, visited)
        outstring += "}\n"

    dot = subprocess.Popen(['dot','-Ttiff'],
                           stdin = subprocess.PIPE,
                           stdout = subprocess.PIPE)
    raw = dot.communicate(input=outstring)[0]
    Image.open(io.BytesIO(raw)).show()

    return outstring



if __name__ == "__main__":
    import optparse

    parser = optparse.OptionParser()
    parser.add_option("-s", "--silent", action="store_true", dest="silent", help="silence output", default=False)
    parser.add_option("-g", "--graph", action="store_true", dest="graph", help="generate graph of use statement inheritance.", default=False)
    parser.add_option("-p", "--path", dest="path", help="set path", default=os.path.dirname(os.path.realpath(__file__))+"/../*/*.F90")

    (options, args) = parser.parse_args()

    tree, namelist = build_trees(glob.glob(options.path))
    out = print_misplaced(tree, namelist, options.silent or options.graph, args)

    if options.graph:

        make_graphs(tree, namelist)

    if out:
        sys.exit(10)
    else:
        sys.exit(0)
