#!/usr/bin/env python

""" 

 Linear_Interp_1D.py
 
 The Linear_Interp_1D class will read a 1d solution field,
 then find the interpolated solution value at specified 
 x coordinate.
 
 The driver main will read a interpolation mesh and 
 then output the interpolated solution field.
 
 The x coordinate assoicated with the solution field must 
 be ascending.
 
 The driver main is used via:
 
 Linear_Interp.py solution_file_name interpolation_mesh_file_name interpolation_solution_file_name
 
"""

import sys
import os
  
class Linear_Interp_1D:
   
   def __init__(self, solution_file_name):
      """ Initialise Linear_Interp_1D class """
      
      # open two input files
      self.solution_file = open(solution_file_name,"r")
      
      # read the solution file  and x coord into lists
      self.read_solution_file()
      
      # check that the coord mesh has ascending x values
      self.check_mesh_has_ascending_x_coord()
      
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
   
   def check_mesh_has_ascending_x_coord(self):
      """ Check that the mesh has an ascending x coordinate """
      
      for node in range(2,len(self.x_coord)):
         
         if self.x_coord[node] <= self.x_coord[node-1]:
            
            print "Issue with solution mesh: x coordinate not ascending for nodes:",node,node-1
            
            sys.exit()
            
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
         
         # Find the x left, s left and right
         x_left = self.x_coord[ele - 1]
         
         s_right = self.solution[ele]
         s_left  = self.solution[ele-1]
         
         # Find delta x
         delta_x = (x - x_left) / (x_right - x_left)
         
         # linear interpolate
         s = (1.0 - delta_x) * s_left + delta_x * s_right
                                    
      return s
         
if __name__ == "__main__":
   
   try:

      solution_file_name = sys.argv[1]
            
      interpolation_mesh_file_name = sys.argv[2]

      interpolation_solution_file_name = sys.argv[3]
            
   except:
    
     print "Usage: Linear_Interp_1D.py solution_file_name interpolation_mesh_file_name interpolation_solution_file_name"
     
     sys.exit()
   
   print "Running Linear_Interp_1D.py ",solution_file_name, interpolation_mesh_file_name, interpolation_solution_file_name
   
   # initialise
   LI = Linear_Interp_1D(solution_file_name)
   
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
    
     solution_interp[x_index] = LI.interp(x_coord_interp[x_index])
   
   # open the interpolated solution output file name
   interpolation_solution_file = open(interpolation_solution_file_name,"w")
   
   # write out the interpolated x_coord_interp solution_interp
   for x_index in range(len(x_coord_interp)): 
      
     interpolation_solution_file.write(str(x_coord_interp[x_index])+" "+str(solution_interp[x_index])+"\n") 
   
   interpolation_solution_file.close()
   
