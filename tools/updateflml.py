#!/usr/bin/env python

# Copyright (C) 2013 Columbia University in the City of New York and others.
#
# Please see the AUTHORS file in the main source directory for a full list
# of contributors.
#
# This file is part of TerraFERMA.
#
# TerraFERMA is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# TerraFERMA is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with TerraFERMA. If not, see <http://www.gnu.org/licenses/>.

import glob
import os
import sys
import argparse
import string
import shutil
from lxml import etree

parser = argparse.ArgumentParser( \
                       description="""Updates flml files so they pass validation.""")
parser.add_argument('file', action='store', metavar='file', type=str, nargs='+',
                    help='specify a filename or expression')
parser.add_argument('-r', '--recursive', metavar='depth', action='store', type=int, dest='recurse', nargs='?', default=None,
                    required=False, const=-1,
                    help='recursively search the directory tree for files (if no depth is specified full recursion will be used)')
args = parser.parse_args()

filenames = []
for f in args.file:

  if args.recurse is None:
    for filename in glob.glob(f):
      ext = filename.split('.')[-1]
      if ext == "flml":
        filenames.append(os.path.join(os.curdir, filename))
      else:
        print "Don't know how to deal with extension "+ext+".  Only know about flmls."
        sys.exit(1)
  else:

    if os.path.isabs(f):
      dirname = os.path.dirname(f)
    else:
      dirname = os.curdir
    dirname = os.path.normpath(dirname)

    for root, dirnames, files in os.walk(dirname, topdown=True):
      depth = string.count(root, os.path.sep)
      for filename in glob.glob1(os.path.join(root, os.path.dirname(f)), os.path.split(f)[-1]):
        ext = filename.split('.')[-1]
        if ext == "flml":
          filenames.append(os.path.join(os.path.join(root, os.path.dirname(f)), filename))
        else:
          print "Don't know how to deal with extension "+ext+".  Only know about flmls."
          sys.exit(1)
      if depth == args.recurse:
        dirnames[:] = []

def update_extrusion_with_layers():
  global changed
  elements = tree.findall("//extrude/regions")
  if len(elements) != 0:
    for element in elements:
      if len(element.findall('layers'))>0:
          continue
      layer = etree.SubElement(element, "layers", attrib={"name": "WholeDepth"})
      for child in element.getchildren():
          print child.tag
          # region_ids remains directly under regions
          if child.tag == 'region_ids' or child.tag == 'layers':
              continue
          # all other items go under the new layer
          layer.append(child)
    changed = True

for filename in filenames:
  print 'Updating: '+filename
  tree = etree.parse(filename)
  global changed
  changed = False
  update_extrusion_with_layers()
  if changed:
    backupfilename = os.path.join(os.path.split(filename)[0], os.path.split(filename)[-1].split('.')[0]+'.flml.bak')
    shutil.copy2(filename, backupfilename)
    tree.write(filename, encoding='utf-8', xml_declaration=True)
