#!/usr/bin/env python3

# This script will plot the final numerical (solid line) and exact (dotted line)
# solutions to phi_x - diff*phi_xx = 1, phi = 0 x = 0,1 (Donea and Huerta pgs 38, 40, a = 1)
# for different Peclet numbers (runs of which should be in separate directories)
# and will give the grid Peclet number expected and as calculated in 
# fluidity (as an average over the values in the final output vtu)

# The directories where the vtus and flml are to be found
# need to be specified (below). N.B. it assumes they are in the same 
# directory (if they are not just alter the paths)

import vtktools
import numpy
import glob
from pylab import *
from lxml import etree
import re


################################################################################################
#----------------------------------------------------------------------------------------------#
################################################################################################

def GetFiles(directory):
# gets list of vtus and sorts them into ascending time order


  filelist_not_sorted = glob.glob('./'+directory+'/*.vtu')
  vtu_nos_not_sorted = []
  for file in filelist_not_sorted:
    vtu_nos_not_sorted.append(int(file.split('.vtu')[0].split('_')[-1]))
  
  vtu_nos_sorted = numpy.argsort(vtu_nos_not_sorted)
  filelist_sorted = []
  for i in vtu_nos_sorted:
    filelist_sorted.append(filelist_not_sorted[i])

  return filelist_sorted

################################################################################################
#----------------------------------------------------------------------------------------------#
################################################################################################

def not_comment(x):
# function to filter stream for use in GetmixingbinsBounds
  return not 'comment' in x.tag

################################################################################################
#----------------------------------------------------------------------------------------------#
################################################################################################

def GetDiffusivity(flml_name):
# returns mixing bin bounds from flml, courtesy of Mr Patrick Farrell

# We will be filtering the children of the elements later,
# to remove comments.

# The spud file to modify
  filename = flml_name

# The path to the node in the tree
  xpath = '/fluidity_options/material_phase[@name="fluid"]/scalar_field[@name="Tracer"]/prognostic/tensor_field[@name="Diffusivity"]/prescribed/value[@name="WholeMesh"]/isotropic/constant'

# Open it up
  tree = etree.parse(open(filename))

  node = tree.xpath(xpath)[0]

  child = filter(not_comment, node.getchildren())[0]

  return float(child.text)

################################################################################################
#----------------------------------------------------------------------------------------------#
################################################################################################


dx = 0.1
a = 1.0

directories  = []
directories.append('1')
directories.append('2')
directories.append('3')
#directories.append('path to directory 1')
#directories.append('path to directory 2')
#directories.append('path to directory 3')
#directories.append('path to directory 4')
# etc

line_colour = ['b','g','r','c','m','y','k']

x_val = [i*dx for i in range(0,11)]

figure(num=None, figsize = (16.5, 11.5))

for directory in directories:
    filelist = GetFiles(directory)

    data = vtktools.vtu(filelist[-1])

    Pec = average(data.GetScalarField('GridPecletNumber'))
    phi = data.GetScalarField('Tracer') 
    nu = GetDiffusivity(directory+'/grid_peclet_number.flml')
    plot(x_val, phi, line_colour[0]+'-', label = str(a*dx/nu)+', '+str(Pec))

    exact = [x-((1-math.exp(x/nu))/(1-math.exp(1.0/nu))) for x in x_val]
    plot(x_val, exact, line_colour[0]+'--') 
    line_colour.remove(line_colour[0])

xlabel('$x$')
ylabel('$\phi$')
grid('True')
legend(loc=0)
savefig('phivsx.png')

