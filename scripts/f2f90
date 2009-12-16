#!/usr/bin/env python
# This utility reads .F files and spits out (hopefully valid) .F90 files.
from optparse import OptionParser
import re
import sys

#####################################################################
# Script starts here.
optparser=OptionParser(usage='usage: %prog [options] <filename>',
                       add_help_option=True,
                       description="""This converts fixed form fortran to free
form. Once you have used it you probably want to re-indent the resulting F90
file (use emacs, for example). This program will reformat comments and
continuation lines appropriately and rewrite all numbered do loops as
do...end do""")

(options, argv) = optparser.parse_args()

if len(argv)<1:
    optparser.print_help()
    sys.exit(1)

ffilename=argv[0]

ffile=file(ffilename,'r')
fortran=ffile.readlines()

do_re=re.compile("^ \s*do\s+(?P<num>\d*)(?P<loop>.*)",re.IGNORECASE)
continue_re=re.compile("^\s*(?P<num>\d*)\s*continue",re.IGNORECASE)
procedure_re=re.compile("^\s*\w*\s*(function|subroutine)\s*\w*",re.IGNORECASE)

# Need this to check if continue lines end multiple numbered do loops
numbered_loop_count = []

for i in range(len(fortran)):

    try:
        # Deal with leading comments.
        if (fortran[i][0] in "cC!*"):
            fortran[i]="!"+fortran[i][1:]
        # Check if line is a comment
        elif (fortran[i].lstrip()[0] == "!"):
            pass
        # Deal with #include
        elif (fortran[i][0] == "#"):
            pass
        # Deal with continuation lines.   
        elif (not fortran[i][5] in " \n"):
            fortran[i]="     &"+fortran[i][6:]
            # Find previous non-comment line.
            j=i-1            
            while (len(fortran[j].lstrip()) == 0 or fortran[j].lstrip()[0] in "!#"):
                j=j-1
            # Deal with trailing comments before continuation lines
            split = fortran[j].split("!")
            if(len(split) > 1):
              fortran[j] = split[0] + "& "
              for section in split[1:]:
                fortran[j] += "!"
                fortran[j] += section
            else:
              fortran[j]=fortran[j][:-1]+"&\n"

    except IndexError:
        # Handle short lines.
        pass
        
    # Reset the number counting for each subroutine/function
    if(procedure_re.match(fortran[i])):
      numbered_loop_count = []
   
    # Deal with numbered do loops.
    domatch = do_re.match(fortran[i])
    if domatch:
        match = domatch.group('num')
        fortran[i] = "      do " + domatch.group('loop')
                    
        if not match == "":
          fortran[i] += "! Was loop " + match
          numbered_loop_count.append(match)
          
        fortran[i] += "\n"
            
    continuematch = continue_re.match(fortran[i])
    if continuematch:
      match = continuematch.group('num')
      if(numbered_loop_count.count(match) > 0):
        fortran[i] = ""
        for j in range(numbered_loop_count.count(match)):
          fortran[i] += "      end do ! Was loop " + match + "\n"
          numbered_loop_count.remove(match)

f90file=file(ffilename+"90",'w')

f90file.writelines(fortran)
