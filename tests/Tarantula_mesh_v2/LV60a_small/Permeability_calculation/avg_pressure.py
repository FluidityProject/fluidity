#!/usr/bin/env python

# Calculate the flux through a plane.

import vtk
import numpy
from math import sqrt



def avg_p(filename):
	origin = (95.0, 0, 0.0)
	normal = (1.0, 0, 0)

        pre_fieldname = "Pressure"

	# Get the unstructured mesh.
	vtkfile = vtk.vtkXMLUnstructuredGridReader()
	vtkfile.SetFileName(filename)
	vtkfile.Update()
	ug = vtkfile.GetOutput()
	
	# Use an implicit plane to do the cutting.
	plane = vtk.vtkPlane()
	plane.SetOrigin(origin)
	plane.SetNormal(normal)
		
	# Set up the cutter.
	cutter = vtk.vtkCutter()
	cutter.SetInput(ug)
	cutter.SetCutFunction(plane)
	cutter.Update()
	
	surface = cutter.GetOutput()
	
	ncells = surface.GetNumberOfCells()
	
	tot_pre = 0.0
	
	for i in range(ncells):

	    n1 = surface.GetCell(i).GetPointId(0)

            pre = surface.GetPointData().GetArray(pre_fieldname).GetValue(n1)
            
            
            f1 = open("output.txt", 'a')
            
            print>>f1, pre

            tot_pre += pre
           
        avg_pre = tot_pre / ncells


	return avg_pre


if __name__=="__main__":
	from numpy import zeros
	import pylab as pl
	f1 = open("output.txt", 'a')


	
        fname = 'cube.vtu'

        average_pressure = avg_p(fname)

        print>>f1, 'AVERAGE PRESSURE', average_pressure

        print>>f1, '********************************************************'



