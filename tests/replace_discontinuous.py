#!/usr/bin/env python

from optparse import OptionParser
import re
import sys
import os.path
import libspud as s

#####################################################################
# Script starts here.
optparser=OptionParser(usage='usage: %prog  <flml filename>',
                       add_help_option=True,
                       description="""This takes an flml file """ + 
                       """and updates all discontinuous meshes.""")
(options, argv) = optparser.parse_args()

if len(argv)<1:
    optparser.print_help()
    sys.exit(1)

flml=argv[0]

s.load_options(flml)

for i in range(s.option_count("/geometry/mesh")):
    path="/geometry/mesh["+`i`+"']/from_mesh"
    
    try:
        continuity=s.get_option(path+"/mesh_continuity")
    except s.SpudKeyError:
        continue

    if (continuity=="discontinuous"):
        try:
            s.add_option(path+"/mesh_shape/element_type")
        except s.SpudNewKeyWarning:
            pass
        try:
            s.set_option(path+"/mesh_shape/element_type",
                         "discontinuous lagrangian")
        except s.SpudNewKeyWarning:
            pass
        
        s.delete_option(path+"/mesh_continuity")
    

s.write_options(flml)
