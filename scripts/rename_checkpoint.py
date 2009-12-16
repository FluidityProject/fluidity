#!/usr/bin/env python

from optparse import OptionParser
import glob
import shutil
import os
import sys

#####################################################################
# Script starts here.
optparser=OptionParser(usage='usage: %prog [options] <base_filename> <index>',
                           add_help_option=True,
                           description="""This takes a list of vtu files in the working directory produced """ +
                           """from a serial checkpointed flml file with names base_filename_checkpoint_i.vtu """ +
                           """for all i and renames them as base_filename_index+i.vtu.  """ +
                           """ """ + 
                           """Can additionally take a list of vtu and pvtu files in the current directory produced """ +
                           """from a checkpointed parallel flml file with names base_filename_checkpoint_i_j.vtu """ +
                           """and base_filename_checkpoint_i.pvtu for all i (index) and j (processor number) """ +
                           """and renames them as base_filename_index+i_j.vtu and base_filename_index+i.pvtu.  """ +
                           """ """ + 
                           """WARNING: This may overwrite files if the backup filenames being written to exist already!""")

optparser.add_option("-v", "--verbose", help="print filenames being moved", action = "store_true", dest = "verbose", default = False)

(options, argv) = optparser.parse_args()

if len(argv)<2:
   optparser.print_help()
   sys.exit(1)

if argv[0][-4:]==".vtu":
   base_filename = os.path.basename(argv[0][:-4])
elif argv[0][-5:]==".pvtu":
   base_filename = os.path.basename(argv[0][:-5])
else:
   base_filename = os.path.basename(argv[0])

verbose = options.verbose

index = int(argv[1])

filelist = glob.glob(base_filename+"_checkpoint_*[0-9].vtu")+glob.glob(base_filename+"_checkpoint_*[0-9].pvtu")
for i in range(len(filelist)):
  if filelist[i][-4:]==".vtu":
    filesplit  = filelist[i].split(".vtu")[0].split(base_filename+"_checkpoint_")[-1].split("_")
    # serial vtus
    if(len(filesplit)==1):
      newindex = index + int(filesplit[0])
      newfilename = base_filename+"_"+str(newindex)+".vtu"
      if(os.path.exists(newfilename)):
        if(verbose): print "backing up", newfilename, "to", newfilename+".bak"
        shutil.move(newfilename, newfilename+".bak")
      if(verbose): print "moving", filelist[i], "to", newfilename
      shutil.move(filelist[i], newfilename)
    # parallel vtus
    elif(len(filesplit)==2):
      newindex = index + int(filesplit[0])
      newfilename = base_filename+"_"+str(newindex)+"_"+filesplit[1]+".vtu"
      if(os.path.exists(newfilename)):
        if(verbose): print "backing up", newfilename, "to", newfilename+".bak"
        shutil.move(newfilename, newfilename+".bak")
      if(verbose): print "moving", filelist[i], "to", newfilename
      shutil.move(filelist[i], newfilename)
  elif filelist[i][-5:]==".pvtu":
    filesplit  = filelist[i].split(".pvtu")[0].split(base_filename+"_checkpoint_")[-1].split("_")
    # parallel pvtus
    if(len(filesplit)==1):
      newindex = index + int(filesplit[0])
      newfilename = base_filename+"_"+str(newindex)+".pvtu"
      if(os.path.exists(newfilename)):
        if(verbose): print "backing up", newfilename, "to", newfilename+".bak"
        shutil.move(newfilename, newfilename+".bak")
      if(verbose): print "moving", filelist[i], "to", newfilename
      shutil.move(filelist[i], newfilename)
      # must also adjust content of pvtu so it points at moved parallel vtus
      pvtufile = file(newfilename,'r')
      lines = pvtufile.readlines()
      pvtufile.close()
      for line in range(len(lines)):
        lineindex=lines[line].find(base_filename+"_checkpoint_"+filesplit[0])
        if(lineindex!=-1):
          processorsplit = lines[line][lineindex:].split(".vtu")[0].split(base_filename+"_checkpoint_")[-1].split("_")
          if(len(processorsplit)==2):
            newline = lines[line][:lineindex]+base_filename+"_"+str(newindex)+"_"+str(processorsplit[-1])+".vtu"+lines[line][lineindex:].split(".vtu")[-1]
            lines[line] = newline
      pvtufile = file(newfilename,'w')
      pvtufile.writelines(lines)
      pvtufile.close()




