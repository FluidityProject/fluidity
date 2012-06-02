#!/usr/bin/env python

""" 

 Compare_Numerical_To_Analytic_1d.py
 
 This script will read a numerical solution field
 (chosen as a keyword on initialisation) from a given vtu file and 
 also read the analytic solution from another given file, 
 then find the numerical solution on the analytic fine 
 point mesh, then compute the various errors and 
 output the results as simple text files that can 
 be plotted with xmgrace. The numerical solution is 
 considered as both FE and CV. The L2 norm error 
 is calculated assuming a uniform analytic solution 
 mesh via a simple approximation = sum(dx*error_node**2)
 where dx = domain_length / number_node_analytic_mesh.
  
 The x coordinate assoicated with the solution field and 
 the analytic solution be ascending.
 
 This is used via:
 
 Compare_Numerical_To_Analytic_1d.py solution_file_name analytic_file_name output_file_name vtu_field_name
 
"""

import sys
import os
import vtktools
  
class Compare_Numerical_To_Analytic_1d:
   
   def __init__(self, solution_file_name, analytic_file_name, output_file_name, vtu_field_name):
      """ Initialise Compare_Numerical_To_Analytic_1d class """
      
      # Set the vtu field name
      self.vtu_field_name = vtu_field_name
      
      # Open the output files
      self.fe_solution_analytic_mesh_txt_output_file = open(output_file_name + "_fe_solution_analytic_mesh.txt","w")
      self.cv_solution_analytic_mesh_txt_output_file = open(output_file_name + "_cv_solution_analytic_mesh.txt","w")
      self.fe_solution_solution_mesh_txt_output_file = open(output_file_name + "_fe_solution_solution_mesh.txt","w")
      self.cv_solution_solution_mesh_txt_output_file = open(output_file_name + "_cv_solution_solution_mesh.txt","w")
      self.fe_error_analytic_mesh_txt_output_file    = open(output_file_name + "_fe_error_analytic_mesh.txt","w")
      self.cv_error_analytic_mesh_txt_output_file    = open(output_file_name + "_cv_error_analytic_mesh.txt","w")
      self.error_norm_txt_output_file                = open(output_file_name + "_error_norm.txt","w")
      
      # open two input files
      self.solution_vtu_input_file = vtktools.vtu(solution_file_name)    
      self.analytic_txt_input_file = open(analytic_file_name,"r") 
      
      # read the solution vtu file into lists
      self.read_solution_vtu_input_file()
      
      # read the analytic file into lists
      self.read_analytic_txt_input_file()
      
      # check that the coord meshes have ascending x values
      self.check_meshes_have_ascending_x_coord()
      
      # if number of analytic x_coord not bigger than 1 then not going to work
      if self.number_analytic_x_coord < 2:
         
         print "Error as x coord of analytic solution has less than 2 values"
         
         sys.exit()

      # if number of solution x_coord not bigger than 1 then not going to work
      if self.number_solution_x_coord < 2:
         
         print "Error as x coord of numerical solution has less than 2 values"
         
         sys.exit()
      
      # initialise fields to calculate
      self.fe_solution_val_analytic_mesh = []
      self.cv_solution_val_analytic_mesh = []
      self.cv_solution_val_solution_mesh = []
      self.fe_error_val_analytic_mesh    = []
      self.cv_error_val_analytic_mesh    = []
      
   def read_solution_vtu_input_file(self):
      """ Read the solution vtu input file """
      
      # Read the coordinates from vtu
      xyz = self.solution_vtu_input_file.GetLocations()
      
      # Extract the x coordinate values and set in self
      self.solution_x_coord = xyz[:,0]
      
      # Deduce number of x coords
      self.number_solution_x_coord = len(self.solution_x_coord)
      
      # Extract the solution field values and set in self
      self.solution_val = self.solution_vtu_input_file.GetScalarField(self.vtu_field_name)
      
   def read_analytic_txt_input_file(self):
      """ Read the analytic txt input file """
      
      # initialise the input lists
      self.analytic_x_coord = []
      self.analytic_val     = []
      
      # go to start of file
      self.analytic_txt_input_file.seek(0)      
      
      # This assume a particular file format
                   
      # initialise counter
      self.number_analytic_x_coord = 0

      for line in self.analytic_txt_input_file:
         
         # skip blank lines (ie a line that cannot be stripped)
         if not line.strip():
            
            continue
         
         else:
            
            pass
         
         # count number of analytic x coord found
         self.number_analytic_x_coord = self.number_analytic_x_coord + 1
         
         self.analytic_x_coord.append(float(line.split()[0]))
         self.analytic_val.append(float(line.split()[1]))
   
   def check_meshes_have_ascending_x_coord(self):
      """ Check that the meshes have ascending x coordinates """
      
      for node in range(2,len(self.analytic_x_coord)):
         
         if self.analytic_x_coord[node] < self.analytic_x_coord[node-1]:
            
            print "Issue with analytic mesh: x coordinate not ascending for nodes:",node,node-1
            
            sys.exit()
      
      for node in range(2,len(self.solution_x_coord)):
         
         if self.solution_x_coord[node] < self.solution_x_coord[node-1]:
            
            print "Issue with solution mesh: x coordinate not ascending for nodes:",node,node-1
            
            sys.exit()
            
   def find_fe_and_cv_solutions_on_analytic_mesh(self):
      """ Find the FE and CV solution on the analytic coordinate mesh """
            
      for anode in range(self.number_analytic_x_coord):
         
         # find left (lower) bounding solution node x coord
         # and deduce the field values on the analytic mesh
                  
         if self.analytic_x_coord[anode] < self.solution_x_coord[0]:
                        
            self.fe_solution_val_analytic_mesh.append(self.solution_val[0])
            self.cv_solution_val_analytic_mesh.append(self.solution_val[0])
            
         elif self.analytic_x_coord[anode] > self.solution_x_coord[-1]:
                             
            self.fe_solution_val_analytic_mesh.append(self.solution_val[-1])
            self.cv_solution_val_analytic_mesh.append(self.solution_val[-1])
            
         else:
                        
            for snode in range(1,self.number_solution_x_coord):
                              
               if self.analytic_x_coord[anode] <= self.solution_x_coord[snode]:
                  
                  cv_surf_x = ((self.solution_x_coord[snode] - self.solution_x_coord[snode-1]) / 2) + \
                              self.solution_x_coord[snode-1]
                  
                  alpha = (self.analytic_x_coord[anode] - self.solution_x_coord[snode-1]) / \
                          (self.solution_x_coord[snode] - self.solution_x_coord[snode-1])
                                                      
                  fe_val = alpha*self.solution_val[snode] + (1.0 - alpha)*self.solution_val[snode-1]
                  
                  if self.analytic_x_coord[anode] <= cv_surf_x:
                     
                     cv_val = self.solution_val[snode-1]
                  
                  else:
                     
                     cv_val = self.solution_val[snode]
                  
                  self.fe_solution_val_analytic_mesh.append(fe_val)
                  
                  self.cv_solution_val_analytic_mesh.append(cv_val)
                  
                  break

   def find_fe_and_cv_solution_errors(self):
      """ Calculate the FE and CV solution errors """
      
      for node in range(self.number_analytic_x_coord):
         
         self.fe_error_val_analytic_mesh.append(abs(self.fe_solution_val_analytic_mesh[node] - self.analytic_val[node]))
         self.cv_error_val_analytic_mesh.append(abs(self.cv_solution_val_analytic_mesh[node] - self.analytic_val[node]))
      
      # Calculate the max error (Linf) for both FE and CV
      self.linf_fe_error = max(self.fe_error_val_analytic_mesh)
      self.linf_cv_error = max(self.cv_error_val_analytic_mesh)
      
      # Calculate the roungh mesh size per node
      dx = (self.analytic_x_coord[-1] - self.analytic_x_coord[0]) / float(self.number_analytic_x_coord)
            
      self.l2_fe_error = 0.0
      self.l2_cv_error = 0.0
      self.l1_fe_error = 0.0
      self.l1_cv_error = 0.0
      
      for node in range(self.number_analytic_x_coord):
         
         self.l2_fe_error = self.l2_fe_error + dx*(self.fe_error_val_analytic_mesh[node]**2.0) 
         self.l2_cv_error = self.l2_cv_error + dx*(self.cv_error_val_analytic_mesh[node]**2.0) 
         self.l1_fe_error = self.l1_fe_error + dx*(self.fe_error_val_analytic_mesh[node]) 
         self.l1_cv_error = self.l1_cv_error + dx*(self.cv_error_val_analytic_mesh[node]) 
               
   def output_fe_and_cv_solutions_and_errors_on_analytic_mesh(self):
      """ Output the FE and CV solutions on the analytic coordinate mesh """
      
      for node in range(self.number_analytic_x_coord):
         
         self.fe_solution_analytic_mesh_txt_output_file.write(str(self.analytic_x_coord[node]) + " " + str(self.fe_solution_val_analytic_mesh[node]) + "\n")
         self.cv_solution_analytic_mesh_txt_output_file.write(str(self.analytic_x_coord[node]) + " " + str(self.cv_solution_val_analytic_mesh[node]) + "\n")
         self.fe_error_analytic_mesh_txt_output_file.write(str(self.analytic_x_coord[node]) + " " + str(self.fe_error_val_analytic_mesh[node]) + "\n")
         self.cv_error_analytic_mesh_txt_output_file.write(str(self.analytic_x_coord[node]) + " " + str(self.cv_error_val_analytic_mesh[node]) + "\n")

   def output_fe_and_cv_solutions_on_solution_mesh(self):
      """ Output the FE and CV solutions on the solution coordinate mesh """
      
      self.cv_solution_solution_mesh_txt_output_file.write(str(self.solution_x_coord[0]) + " " + str(self.solution_val[0]) + "\n")      
      
      for node in range(self.number_solution_x_coord):
         
         self.fe_solution_solution_mesh_txt_output_file.write(str(self.solution_x_coord[node]) + " " + str(self.solution_val[node]) + "\n")
         
         cv_surf_x = ((self.solution_x_coord[node] - self.solution_x_coord[node-1]) / 2) + \
                     self.solution_x_coord[node-1]
         
         if node > 0:
            
            self.cv_solution_solution_mesh_txt_output_file.write(str(cv_surf_x) + " " + str(self.solution_val[node-1]) + "\n")
            self.cv_solution_solution_mesh_txt_output_file.write(str(cv_surf_x) + " " + str(self.solution_val[node]) + "\n")
            self.cv_solution_solution_mesh_txt_output_file.write(str(self.solution_x_coord[node]) + " " + str(self.solution_val[node]) + "\n") 
   
   def output_error_norms(self):
      """ Ouput error norms """
      
      self.error_norm_txt_output_file.write("FE Linf error: "+str(CNTA1D.linf_fe_error)+"\n")
      self.error_norm_txt_output_file.write("FE L2 error: "+str(CNTA1D.l2_fe_error)+"\n")
      self.error_norm_txt_output_file.write("FE L1 error: "+str(CNTA1D.l1_fe_error)+"\n")
      self.error_norm_txt_output_file.write("CV Linf error: "+str(CNTA1D.linf_cv_error)+"\n")
      self.error_norm_txt_output_file.write("CV L2 error: "+str(CNTA1D.l2_cv_error)+"\n")
      self.error_norm_txt_output_file.write("CV L1 error: "+str(CNTA1D.l1_cv_error)+"\n")
         
   def close_files(self):
      """ Close opened files """
      
      self.fe_solution_analytic_mesh_txt_output_file.close()
      self.cv_solution_analytic_mesh_txt_output_file.close()
      self.fe_solution_solution_mesh_txt_output_file.close()
      self.cv_solution_solution_mesh_txt_output_file.close()
      self.fe_error_analytic_mesh_txt_output_file.close()
      self.cv_error_analytic_mesh_txt_output_file.close()      
      self.error_norm_txt_output_file.close()
      
      self.analytic_txt_input_file.close()
     
if __name__ == "__main__":
   
   try:

      solution_file_name = sys.argv[1]
            
      analytic_file_name = sys.argv[2]

      output_file_name   = sys.argv[3]

      vtu_field_name     = sys.argv[4]
            
   except:
    
     print "Usage: Compare_Numerical_To_Analytic_1d.py solution_file_name analytic_file_name output_file_name vtu_field_name"
     
     sys.exit()
   
   print "Running Compare_Numerical_To_Analytic_1d.py ",solution_file_name, analytic_file_name, output_file_name, vtu_field_name
   
   # initialise
   CNTA1D = Compare_Numerical_To_Analytic_1d(solution_file_name,analytic_file_name,output_file_name,vtu_field_name)
      
   # find the fe and cv solution fields on the analytic mesh
   CNTA1D.find_fe_and_cv_solutions_on_analytic_mesh()
   
   # find the fe and cv errors
   CNTA1D.find_fe_and_cv_solution_errors()
   
   # output the error norms to a file
   CNTA1D.output_error_norms()   
   
   # output the fe and cv solution and error fields on the analytic mesh
   CNTA1D.output_fe_and_cv_solutions_and_errors_on_analytic_mesh()  
   
   # output the fe and cv solution fields on the solution mesh
   CNTA1D.output_fe_and_cv_solutions_on_solution_mesh()  
   
   # close opened files
   CNTA1D.close_files()  
