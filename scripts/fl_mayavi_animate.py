#!/usr/bin/env python
#    Copyright (C) 2006 Imperial College London and others.
#
#    Please see the AUTHORS file in the main source directory for a full list
#    of copyright holders.
#
#    Prof. C Pain
#    Applied Modelling and Computation Group
#    Department of Earth Science and Engineering
#    Imperial College London
#
#    C.Pain@Imperial.ac.uk
#
#    This library is free software; you can redistribute it and/or
#    modify it under the terms of the GNU Lesser General Public
#    License as published by the Free Software Foundation; either
#    version 2.1 of the License, or (at your option) any later version.
#
#    This library is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#    Lesser General Public License for more details.
#
#    You should have received a copy of the GNU Lesser General Public
#    License along with this library; if not, write to the Free Software
#    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307
#    USA

import mayavi
import getopt, sys, os, time
import mayavi.Common, mayavi.Base.Objects, mayavi.Base.ModuleManager, mayavi.Base.DataVizManager
import mayavi.Sources.VtkDataReader, mayavi.Sources.PLOT3DReader, mayavi.Sources.VRMLImporter
import mayavi.Sources.mv3DSImporter, mayavi.Sources.VtkData
import string

def usage ():
    msg="""Usage:\n\nfl_mayavi_animate -f <fl-file> [-z <mv-file>] [options] <start> <finish>

    Where <start> <finish> are the first and last dump-id's of the
    files to be rendered and <mv-file> is the MayaVi visualization
    file to be used.

    Valid options are one or more of the following:

    -s
    --fps

         Frames per second in output animation (default: 5)

    -f
    --flfile

         Fluidity file

    -z
    --vizfile

         The MayaVi visualization file that is to be used for all the
         data-files

    -r
    --recycle

         Recycles existing frames - doesn't overwrite

    --gif

         Generate animated gif (default)

    --png

         Generate pngs but no movie (gif or avi)

    --avi

         Generate AVI file
    """
    return msg

def get_file (root, id):
    # Check for vtk file
    if os.path.isfile( root+"_%d"%id+".vtk" ):
        return root+"_%d"%id+".vtk"

    # Check for vtu file
    if os.path.isfile( root+"_%d"%id+".vtu" ):
        return root+"_%d"%id+".vtu"

    # Check for pvtu file
    if os.path.isfile( root+"_%d"%id+".pvtu" ):
        return root+"_%d"%id+".pvtu"

    raise Exception("id " + str(id) + " not found")

def main ():
  try:
    opts, args = getopt.getopt(sys.argv[1:], "hf:z:rs:", ["help", "flfile", "vizfile=", "", "recycle", "gif", "avi", "fps", "png"])
  except getopt.GetoptError:
    # print help information and exit:
    print "ERROR: Bad arguments!"
    print usage()
    sys.exit(2)

  # Collect options.
  recycle = 0
  flfile = None
  mvfile = None
  movie = "gif"
  fps = 5
  maximize = False
  for o, a in opts:
    if o in ("-h", "--help"):
      print usage()
      sys.exit()
    if o in ("-z", "--vizfile"):
      mvfile = a
    if o in ("-r", "--recycle"):
      recycle = 1
    if o in ("-f", "--flfile"):
      flfile = a
    if o == "--avi":
      movie = "avi"
    if o == "--png":
      movie = "non"
    if o in ("-s", "--fps"):
      try:
        fps = float(a)
        assert(fps > 0.0)
      except:
        print "ERROR: Bad arguments!"
        print usage()
        sys.exit(2)
      
  # Check for fluidity file
  if not flfile:
    print "ERROR: No fluidity file given."
    print usage()
    sys.exit()

  # Check for MayaVi visualisation file
  if not mvfile:
    # Guess
    mvfile = flfile+".mv"
    if not os.path.isfile( mvfile ):
      print "ERROR: No MayaVi visualisation file found."
      print usage()
      sys.exit()


  try:
    start  = string.atoi( args[0] )
  except:
    print "ERROR: Bad starting id %s"%args[0]
    print usage()
    sys.exit()

  try:
    finish = string.atoi( args[1] )
  except:
    print "ERROR: Bad last id %s"%args[1]
    print usage()
    sys.exit()

  # Movie file that we'll be writting to
  moviefile = flfile+"."+movie

  # instantiate a MayaVi
  v = mayavi.mayavi()

  # load the visualisation
  v.load_visualization(mvfile)

  # grab the DataVizManager list
  dvms = v.get_dvm_names()

  # list frames to be animated
  frames = ""
  for id in range(start, finish+1):
    # Get frame name
    frame='./'+flfile+'%04d.png'%id
    if (movie=="avi"):
      frames = frames+","+frame
    else:
      frames = frames+" "+frame

    # Are we recycling frames
    if recycle:
      if os.path.isfile(frame):
        continue

    # Make a vtk file
    datafile = get_file (flfile, id)

    # go through all the DVM's
    for i in dvms:

      # grab a handle to the DVM
      dvm = v.mayavi.data_viz_mgr[i]

      # follow the pipeline, and load in the new data file
      ds  = dvm.get_data_source ()
      rdr = ds.get_reader ()
      rdr.SetFileName (datafile)

      ds.reread_file ()
      v.Render ()

      # write this image to the disk
      v.renwin.save_png (frame)

  # Bail out if no movie desired
  if (movie=="non"):
      return 0

  # Now generate the avi file
  pid          = os.getpid ()
  tmpfile_name = os.environ['PWD']+"/.tmpmvfile%d"%pid
  tmpfile      = open(tmpfile_name, 'w')

  tmpfile.write( "#!/bin/sh\n")
  if (movie=="avi"):
    tmpfile.write( "mencoder mf://"+flfile+"*.png -mf fps=" + str(fps) + ":type=png -ovc copy -o "+flfile+".avi\n")
  else:
    tmpfile.write( "convert -delay " + str(100.0 / fps) + " -loop 0 "+frames+" "+flfile+".gif\n" )
  tmpfile.write( "rm -f %s\n"%tmpfile_name )
  tmpfile.close()
  os.chmod(tmpfile_name, 0700)

  # Generate animation
  if os.fork() == 0:
    os.execl( tmpfile_name, "" )
  else:
    os.wait()

  return 0

if __name__ == "__main__":
    main()
