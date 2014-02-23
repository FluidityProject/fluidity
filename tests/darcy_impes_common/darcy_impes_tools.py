#!/usr/bin/env python

"""darcy_impes_tools.py

For interpolating between simulated fields and quasi-analytic solutions,
reading/writing reports of error norms and convergence rates, and batch
processing the darcy_impes tests.

About the interpolations: it is straightforward to compare the the 1D
simulations to the expected 1D quasi-analytical solution at the latter's
point locations; a piecewise-linear interpolant is constructed between
the simulation's points (see computeAbsErrors).  This is difficult to
extend to 2D and 3D, but not if the interpolation is done the other way
round - the analytical solution is used to interpolate at each
simulation point (interpUsingAnalytic).  The jump discontinuity of the
analytical solution will be ill-represented, but we can still expect
first order convergence of the l1-errors.  However, l2-errors are
deprecated.  Their reported convergence rates can be seen to be much
less than their l1 counterparts, presumably as a result of the squaring
of errors at the discontinuity.

"""

import sys
sys.path.append('../../python')

import fluidity_tools
import vtktools
import numpy
import re
import subprocess
import os


## FOR PUBLIC USE 

def interpUsingAnalytic(x, analyticFilename):
    """Interpolates at point x using the profile given by
    analyticFilename."""
    [xA, vA] = readAnalyticSolution(analyticFilename)
    return numpy.interp(x, xA, vA)

       
def findRate1D(key):
    """Opens and searches a 1D rates report and returns the value associated
    with key.
    """
    with open("error_rates_1d.txt", "r") as f:
        v = find(f, key)
        return v

        
def findRateND(key):
    """Opens and searches an ND rates report and returns the value associated
    with key.
    """
    with open("error_rates_nd.txt", "r") as f:
        v = find(f, key)
        return v

        
class TestHelper:
   def __init__(self, numericalFilenameStem, modelNameList,
                gridNameListPerDimension, analyticFilenameStem,
                fieldNameList):
       """Class to help run and postprocess tests.  numericalFilenameStem is a
       string to be appended with the elements in modelNameList and
       gridNameList, where the latter is one of three elements in the
       list gridNameListPerDimension.  analyticFilenameStem is a string to
       be appended with shorthand names for the fields in question,
       which in turn are translated from the elements in fieldNameList.
       If a fieldName is 'Phase2::Saturation', then the shorthand is
       saturation2.  """
       # constants
       self.normList = [1, 2]
       self.verbose = True
       # initialised from args
       self.numericalFilenameStem = numericalFilenameStem
       self.modelNameList = modelNameList
       self.gridNameListPerDimension = gridNameListPerDimension
       self.analyticFilenameStem = analyticFilenameStem
       self.fieldNameList = fieldNameList

       
   def computeAbsErrors(self, numericalSuffix, analyticSuffix, fieldName):
       """Computes errors between interpolated numerical values and analytic
       point values.  Also returns mean grid spacing."""
       sA = self.analyticFilenameStem+'_'+analyticSuffix+'.txt'
       sN = self.numericalFilenameStem+'_'+numericalSuffix+'_1.vtu'
       [xA, vA] = readAnalyticSolution(sA)
       [xN, vN] = readNumericalSolution(sN, fieldName)
       vNInterp = numpy.interp(xA, xN, vN)
       return numpy.abs(vNInterp - vA), (xA[-1] - xA[0])/(len(xA) - 1)

       
   def generateErrorReport1D(self):
       """Compares 1D numerical solutions with the appropriate analytical
       solutions, writing the error norms in a report.  The output
       report has two columns: the first has identifiers in the format
       modelName_1D_gridName_fieldID_normtype; the second has the
       values of the error norms.  """
       gridNameList = self.gridNameListPerDimension[0]
       with open("error_norms_1d.txt", "w") as f:
   
           # loop over numerical solutions
           for modelName in self.modelNameList:
               for gridName in gridNameList:
                   s1 = modelName+'_1d_'+gridName
   
                   # loop over fieldnames
                   for fieldName in self.fieldNameList:
                       # translate the fieldName to shorthand form
                       s2 = str.lower(re.sub('Phase(.)::(.*)', '\\2\\1',
                                             fieldName))
                       [errs, dx] = self.computeAbsErrors(s1, s2, fieldName)
                       
                       # loop over norm types
                       for i in self.normList:
                           # print the ID and value
                           key = '{0}_{1}_l{2:g}'.format(s1,s2,i)
                           v = computeNorm(errs, dx, i)
                           f.write('{0}{1:24.12e}\n'.format(key,v))
   
                           
   def generateConvergenceReport1D(self):
       """Can be run after generateErrorReport1D.  N.B. Has a
       different ordering and formatting of the IDs on the LHS."""
       gridNameList = self.gridNameListPerDimension[0]
       with open("error_rates_1d.txt", "w") as rates:
           with open("error_norms_1d.txt", "r") as norms:
               for modelName in self.modelNameList:
                   for fieldName in self.fieldNameList:
                           s2 = str.lower(re.sub('Phase(.)::(.*)', '\\2\\1',
                                                 fieldName))
                           for i in self.normList:
                               for j, gridName in enumerate(gridNameList):
                                   s1 = modelName+'_1d_'+gridName
                                   key = '{0}_{1}_l{2:g}'.format(s1,s2,i)
                                   v = find(norms, key)
                                   
                                   if gridName!='A':
                                       # print a new ID and the rate
                                       s3 = gridNameList[j-1]+gridName
                                       key = modelName+'_1d_'+s2+'_l'+str(i)+'_'+s3
                                       r = numpy.log2(v0/v)
                                       rates.write('{0}{1:12.6f}\n'.format(key,r))
                                   v0 = v

   
   def generateConvergenceReportND(self, dimensionList):
       """Produces a similar output to generateConvergenceReport1D, but uses
       information associated with multidimensional VTU files incorporating
       interpUsingAnalytic.  This is a more versatile approach as
       explained above.
       """
       with open("error_rates_nd.txt", "w") as rates:
           for modelName in self.modelNameList:
               # loop over dimensions
               for dim in dimensionList:
                   gridNameList = self.gridNameListPerDimension[dim-1]
                   
                   for fieldName in self.fieldNameList:
                       phaseIndex = re.sub('Phase(.)::(.*)', '\\1', fieldName)
                       subFieldName = re.sub('Phase(.)::(.*)', '\\2', fieldName)
                       s2 = str.lower(subFieldName)+phaseIndex
                       
                       for i in self.normList:
                           if i==1: calcType = "integral"
                           elif i==2: calcType = "l2norm"
   
                           for j, gridName in enumerate(gridNameList):
                               s1 = modelName+'_'+str(dim)+'d_'+gridName
                               v = fluidity_tools.stat_parser(self.numericalFilenameStem+'_'+s1+'.stat')['Phase'+phaseIndex]['Analytic'+subFieldName+'Error'][calcType][-1]
                               
                               if gridName!='A':
                                   # print a new ID and the rate
                                   s3 = gridNameList[j-1]+gridName
                                   key = modelName+'_'+str(dim)+'d_'+s2+'_l'+str(i)+'_'+s3
                                   r = numpy.log2(v0/v)
                                   rates.write('{0}{1:12.6f}\n'.format(key,r))
                               v0 = v
       
   
   def generateReports(self):
       """Calls all report generators."""
       self.generateErrorReport1D()
       self.generateConvergenceReport1D()
       self.generateConvergenceReportND([2, 3])


   def processFolder(self):
       """Preprocesses input files before running simulations for each
       variant."""
       if self.verbose: print "Processing "+self.numericalFilenameStem+":"
       for modelName in self.modelNameList:
           if self.verbose: print "  "+modelName+":"
   
           # loop over dimensions
           for i, gridNameList in enumerate(self.gridNameListPerDimension):
               dim = i + 1
               if self.verbose: print "    "+str(dim)+"d:"
               s1 = self.numericalFilenameStem+'_'+modelName+'_'+str(dim)+'d'
       
               # loop over grids
               for gridName in gridNameList:
                   if self.verbose: print "      grid "+gridName
                   s2 = s1+'_'+gridName+'.diml'
                   # for the first grid, an options file should exist.
                   # For the other grids, adapt it
                   if gridName!='A':
                       with open(s1+'_A.diml', 'r') as f_A:
                           f = open(s2, 'w')
                           for line in f_A:
                               f.write(line.replace('_A', '_'+gridName))
                           f.close()
                   # process and clean up 
                   subprocess.call(["../../bin/darcy_impes", s2])
                   if gridName!='A': os.remove(s2)
       

## HELPER FUNCTIONS 
    
def readAnalyticSolution(filename):
    """Converts a file containing a column of x-coordinates and a column of
    values into corresponding lists."""
    x = []; v = [];
    with open(filename, "r") as f:
       f.seek(0)
       for l in f:
          if not l.strip(): continue
          x.append(float(l.split()[0]))
          v.append(float(l.split()[1]))
    if not numpy.all(numpy.diff(x) > 0): 
        print ("Issue with analytic mesh: "
               "x-coordinates not monotonically increasing.")
        sys.exit()
    return x, v

    
def readNumericalSolution(filename, fieldName):
    """Converts a vtu file into a list of x-coordinates and field values."""
    f = vtktools.vtu(filename)
    x = f.GetLocations()[:,0]
    v = f.GetScalarField(fieldName)
    if not numpy.all(numpy.diff(x) > 0): 
        print ("Issue with numerical mesh: "
               "x-coordinates not monotonically increasing.")
        sys.exit()
    return x, v


def computeNorm(errors, gridSpacing, normIndex):
    """Computes a norm.  Uniform, linear elements are assumed, so the
    integral over the domain can be approximated by the trapezium rule
    (for the L1 norm at least).
    """
    E = numpy.array(errors)
    if normIndex == 1: 
        I = (E[0]/2 + sum(E[1:-1]) + E[-1]/2)*gridSpacing
    elif normIndex == 2:
        # does not make much difference, but previously I = (E[0]**2/2 +
        # sum(E[1:-1]**2) + E[-1]**2/2)*gridSpacing
        I = (E[0]**2 + sum(E[:-1]*E[1:]) + 2*sum(E[1:-1]**2) + \
                 E[-1]**2)*gridSpacing/3
    else: 
        print("Only l1 or l2 norms may be computed; i.e. normIndex "
              "must be 1 or 2.")
        sys.exit()
    return I**(1/float(normIndex))
   
   
def find(report, key):
    """Searches a report and returns the value associated with key
    (ID). """
    report.seek(0)
    v = 0.
    for l in report:
        if not l.strip(): continue
        if l.split()[0] == key: 
            v = float(l.split()[1])
            break
    return v
