#!/usr/bin/env python

""" 

 BL_Cubic_Interp.py
 
 The BL_Cubic_Interp class will read a solution field,
 as well as the hermite cubic interpolation coefficient
 then find the interpolated solution value at specified 
 x coordinate.
 
 The driver main will read a interpolation mesh and 
 then output the interpolated solution field
 
 The driver main is used via:
 
 BL_Cubic_Interp.py solution_file_name poly_coeff_file_name interpolation_mesh_file_name interpolation_solution_file_name
 
"""

import sys
import os
  
class BL_Cubic_Interp:
   
   def __init__(self, solution_file_name, poly_coeff_file_name):
      """ Initialise BL_Cubic_Interp class """
      
      # open two input files
      self.solution_file = open(solution_file_name,"r")
      self.poly_coeff_file = open(poly_coeff_file_name,"r")
      
      # read the solution file into lists
      self.read_solution_file()
      
      # read the x coord and the poly coeff set into lists
      self.read_poly_coeff_file()
            
      # if num_x_coord != num_poly_coeff_sets + 1 something wrong
      if self.num_x_coord != self.num_poly_coeff_sets + 1:
         
         print "Error for input file ",self.poly_coeff_file
         print "num_x_coord /= num_poly_coeff_sets + 1"
         print "num_x_coord ",self.num_x_coord
         print "num_poly_coeff_sets ",self.num_poly_coeff_sets
         
         sys.exit()
      
      # if number of x_coord not bigger than 1 then not going to work
      if self.num_x_coord < 2:
         
         print "Error as x coord of solution has less than 2 values"
         
         sys.exit()
      
   def read_solution_file(self):
      """ Read the solution input file """
      
      # initialise the input lists
      self.x_coord = []
      self.solution = []
      
      # go to start of file
      self.solution_file.seek(0)      
      
      # This assume a particular file format
                   
      # initialise counter
      self.num_x_coord = 0

      for line in self.solution_file:
         
         # skip blank lines (ie a line that cannot be stripped)
         if not line.strip():
            
            continue
         
         else:
            
            pass
         
         # count num_x_coord found
         self.num_x_coord = self.num_x_coord + 1
         
         self.x_coord.append(float(line.split()[0]))
         self.solution.append(float(line.split()[1]))
   
   def read_poly_coeff_file(self):
      """ Read the poly_coeff input file """
      
      # initialise the input lists
      self.poly_coeff = []
      
      # go to start of file
      self.poly_coeff_file.seek(0)      
      
      # This assume a particular file format
                   
      # initialise counter
      self.num_poly_coeff_sets = 0

      for line in self.poly_coeff_file:
         
         # skip blank lines (ie a line that cannot be stripped)
         if not line.strip():
            
            continue
         
         else:
            
            pass
         
         # count poly_coeff_sets found
         self.num_poly_coeff_sets = self.num_poly_coeff_sets + 1
         
         # first append an zerod list (lists within list)
         self.poly_coeff.append([0.0,0.0,0.0,0.0])

         # add the coeff into each local list
         self.poly_coeff[self.num_poly_coeff_sets-1][0] = float(line.split()[0])
         self.poly_coeff[self.num_poly_coeff_sets-1][1] = float(line.split()[1])
         self.poly_coeff[self.num_poly_coeff_sets-1][2] = float(line.split()[2])
         self.poly_coeff[self.num_poly_coeff_sets-1][3] = float(line.split()[3])
         
   def interp(self, x):
      """ Find the interpolated value at coordinate x """        
      
      # bound the interpolation to the known x_coord min and max values
      
      if x <= self.x_coord[0]:
         
         s = self.solution[0]
      
      elif x >= self.x_coord[-1]:
         
         s = self.solution[-1]
      
      else:
         
         # find which ele (or interval) this x resides in
         for ele in range(1,len(self.x_coord)):
            
            x_right = self.x_coord[ele]
            
            if x <= x_right:
               
               break
         
         # intialise s
         s = 0.0
         
         # interpolate using the cubic hermite polynomials
         for p in [1,2,3,4]:
            
            s = s + self.poly_coeff[ele-1][p-1] * (x - self.x_coord[ele-1])**(4-p)
                  
      return s
         
if __name__ == "__main__":
   
   try:

      solution_file_name = sys.argv[1]
      
      poly_coeff_file_name = sys.argv[2]
      
      interpolation_mesh_file_name = sys.argv[3]

      interpolation_solution_file_name = sys.argv[4]
            
   except:
    
     print "Usage: BL_Cubic_Interp.py solution_file_name poly_coeff_file_name interpolation_mesh_file_name interpolation_solution_file_name"
     
     sys.exit()
   
   print "Running BL_Cubic_Interp.py ",solution_file_name, poly_coeff_file_name, interpolation_mesh_file_name, interpolation_solution_file_name
   print "where poly_coeff_file_name is an output from JGomes Octave program"
   
   # initialise
   BL = BL_Cubic_Interp(solution_file_name, poly_coeff_file_name)
   
   # initialise the interpolated x coord and solution lists
   x_coord_interp = []
   solution_interp = []
   
   # open the file which has the solution interpolation mesh
   interpolation_mesh_file = open(interpolation_mesh_file_name,"r")
   
   # read the interpolated solution x coord and initialise interpolated solution
   for line in interpolation_mesh_file:
      
      # skip blank lines (ie a line that cannot be stripped)
      if not line.strip():
            
         continue
         
      else:
            
         pass
      
      x_coord_interp.append(float(line.split()[0]))         
       
      solution_interp.append(0.0)
   
   # find the interpolated solution for each of the required points
   for x_index in range(len(x_coord_interp)): 
    
     solution_interp[x_index] = BL.interp(x_coord_interp[x_index])
   
   # open the interpolated solution output file name
   interpolation_solution_file = open(interpolation_solution_file_name,"w")
   
   # write out the interpolated x_coord_interp solution_interp
   for x_index in range(len(x_coord_interp)): 
      
     interpolation_solution_file.write(str(x_coord_interp[x_index])+" "+str(solution_interp[x_index])+"\n") 
   
   interpolation_solution_file.close()
   
