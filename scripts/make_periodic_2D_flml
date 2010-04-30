#!/usr/bin/env python
#encoding: utf-8

#This script takes in a 2+1 periodic flml and converts it into a 2D periodic flml (with no fields).
#This is to allow a parallel periodic 2D mesh to be made using flredecomp
#which can then be read into the original 2+1 periodic flml and then run in parallel 
#(see test case periodic_2plus1_parallel)
#a.mcvicar 30/04/10

from lxml import etree
import optparse

#input file is name of 2plus1 periodic flml
optionParser = optparse.OptionParser( \
  usage = "%prog [OPTIONS] ... PROJECT", \
  add_help_option = True, \
  description = "makes a periodic 2D flml from a periodic 2plus1 flml")

optionParser.add_option("-o", "--optimisation", action = "store_true", dest = "optimisation", help = "Enable optimisation", default = False)
optionParser.add_option("-v", "--verbose", action = "store_true", dest = "verbose", help = "Verbose mode", default = False)

opts, args = optionParser.parse_args()

filename = args[0]

print 'file loaded = ', filename[:-5]

#load filename
tree=etree.parse(filename)
root=tree.getroot()

#list of elements that are allowed
allowed=['simulation_name','geometry','timestepping','io','problem_type']

#find tags of elements that need to be changed
todeletefield=[]
todelete=[]
count=0
for child in root.getchildren():
  if (child.tag == 'material_phase'):
    mp_count = count
    for field in root[count].getchildren():
      todeletefield.append(field)
  elif(child.tag == 'geometry'):
    geo_count = count
  elif not(child.tag in allowed):
    todelete.append(child)
  count=count+1

#Note: have to loop separately as uncertain of the order added in the flml
#locate where 2D mesh is loaded
for i in range(len(root[geo_count].getchildren())):
  for j in range(len(root[geo_count][i].getchildren())):
    if (root[geo_count][i][j].tag == 'from_file'):
      mesh_file=root[geo_count][i].get('name')
      file_load = [i,j]

#locate where periodic boundaries are applied
#locate mesh that has been extruded
todeletefield2=[]
for i in range(len(root[geo_count].getchildren())):
  if (root[geo_count][i].tag != 'mesh') and (root[geo_count][i].tag != 'dimension') and (root[geo_count][i].tag != 'quadrature'):
    todeletefield2.append(root[geo_count][i])
  for j in range(len(root[geo_count][i].getchildren())):
    for k in range(len(root[geo_count][i][j].getchildren())):
      if (root[geo_count][i][j][k].tag == 'extrude'):
	#extrude_loc = [i,j,k]
	todeletefield2.append(root[geo_count][i])
      if (root[geo_count][i][j][k].tag == 'periodic_boundary_conditions'):
	if (root[geo_count][i][j][0].get('name')==mesh_file):
	  periodic_bd_loc = [i,j,k]
	  periodc_bd_mesh =root[geo_count][i].get('name')
          #get coordinate map of periodic boundary
	  for m in range(len(root[geo_count][i][j][k].getchildren())):
	    if (root[geo_count][i][j][k][m].tag=='coordinate_map'):
	      coordinate_map=root[geo_count][i][j][k][m][0].text

#field = root[geo_count].getchildren()[extrude_loc[0]]	  
#remove unwanted fields in geometry
for field in todeletefield2:
  root[geo_count].remove(field) 
 
#remove all material phase fields
for field in todeletefield:
  root[mp_count].remove(field) 

#remove all unwanted elements (not in allowed list)
for element in todelete:
  root.remove(element)

#set dimension to 2D
dim= root.xpath("geometry/dimension")
dimension=dim[0].getchildren()[0]
dimension.text='2'

#ensure meshes are pointing to 2D mesh
meshes=root.xpath("geometry/mesh/from_mesh/mesh")
for i in range(len(meshes)):
  if meshes[i].get('name') != mesh_file:
    meshes[i].set('name',periodc_bd_mesh)

#adjust periodic boundaries from 3D to 2D
corodinate_map=root.xpath("geometry/mesh/from_mesh/periodic_boundary_conditions/coordinate_map")
for i in range(len(corodinate_map)):
  corodinate_map[i][0].text=coordinate_map

#Note the inverse coordinate map is not changed as this is not required by periodise or flredecomp

tree.write(filename[:-5] + "_2D.flml")
