# some of these functions below where copied from
# svn+ssh://amcg.ese.ic.ac.uk/svn/beans/trunk/tools/Beans_Tools.py

import sys
import os
import re

# --------------------------------------------------------------------------------------------------------------------
def compare_two_files_dump_seris(file_series_name1,file_series_name2,dump_start,dump_end,rel_tolerance,start1,end1,start2,end2):
    # compare two files seris dumps where each dump  contain two columns list of numbers 
    # to within a relative tolerance    
    # the pass number(1=all pass, 0=any fail)
    
    # the dump series should be called file_series_name1.d.1 etc
    
    print "Compare two files dump seris"
    
    pass_number = 1
    
    file_series_name3 = str(file_series_name1)
    file_series_name4 = str(file_series_name2)    
    print "Input files series:", file_series_name3,file_series_name4
    
    # loop from the start to the end dump and compare the two files
    dump_difference = dump_end - dump_start
    maxdiff_alldump = 0.0
    maxdiff_dumpnum = 0
    for dump_num in range(0,dump_difference+1):
              
       filein_1 = open(file_series_name3+".d."+str(dump_num+1),"r")    
       filein_2 = open(file_series_name4+".d."+str(dump_num+1),"r")    
    
       # count the number of lines in file
       all_lines_1 = filein_1.readlines()
       all_lines_2 = filein_2.readlines()
       number_lines_1 = len(all_lines_1)
       number_lines_2 = len(all_lines_2)    
    
       if number_lines_1 == number_lines_2:
          number_lines = number_lines_1  
       else:
          pass_number = 0
          number_lines = 0
          print "Error, two files have different number of lines",number_lines_1,number_lines_2
    
       # go back to start of file
       filein_1.seek(0)
       filein_2.seek(0)
    
       # Read in files line by line ======================================
       # start with zero length lists then append in results on the go
       value_1 = []     
       value_2 = [] 
    
       # python list starts from 0 --> [0,1,2,3,...]
       # then range(1,5) = [1,2,3,4] so dont include -1 at end      
       for i in range(0,number_lines):                       

          line_1 = filein_1.readline()
          line_2 = filein_2.readline()
       
          line_split_1 = line_1.split()
          line_split_2 = line_2.split()
    
          value_1.append( float(line_split_1[1]) )    
          value_2.append( float(line_split_2[1]) )
        
       # close the file
       filein_1.close()
       filein_2.close()
        
       # Pass or Fail Test ===============================================
       # Now compare 
       # pass_number = 1 --> pass
       # pass_number = 0 --> fail
    
       difference1 = end1 - start1
       difference2 = end2 - start2
        
       if difference1 == difference2:
          difference = difference1  
       else:
          pass_number = 0
          difference = 0
          print "Error, two seperate ranges for the files asked for"
    
       maxdiff = 0.0
       # python list starts from 0 --> [0,1,2,3,...]
       # then range(1,5) = [1,2,3,4] so dont include -1 at end  
       for i in range(0,difference):
        
          j = i - 1 + start1
          k = i - 1 + start2
       
          if abs(value_2[j]) > 0.0:
             maxdiff = max(maxdiff,abs(value_1[j] - value_2[k])/abs(value_2[j]))
            
          if abs(value_1[j] - value_2[k]) > rel_tolerance*abs(value_2[j]):   
             pass_number = 0
             print "*****Failed on value:",j+1,k+1
             print "     value_1,value_2:",value_1[j],value_2[k] 
             print "     For dump number:",dump_num + 1    
       
       if maxdiff > maxdiff_alldump:
          maxdiff_dumpnum = dump_num + 1
          
       maxdiff_alldump = max(maxdiff_alldump,maxdiff)
       
       if pass_number == 0:
          exit

    print "Relative tolerance, Max relative difference over all dumps and associated dump number",rel_tolerance,maxdiff_alldump,maxdiff_dumpnum
    print "Finish Compare two file dump series: pass_number(1=pass,0=fail)",pass_number  
                                    
    return pass_number

# --------------------------------------------------------------------------------------------------------------------

def compare_two_files(file_name1,file_name2,rel_tolerance,start1,end1,start2,end2):
    # compare two files that contain two columns list of numbers 
    # to within a relative tolerance    
    # the pass number(1=all pass, 0=any fail)
    
    print "Compare two files"
    
    pass_number = 1
    
    file_name3 = str(file_name1)
    file_name4 = str(file_name2)    
    print "Input files:", file_name3,file_name4
    
    filein_1 = open(file_name3,"r")    
    filein_2 = open(file_name4,"r")    
    
    # count the number of lines in file
    all_lines_1 = filein_1.readlines()
    all_lines_2 = filein_2.readlines()
    number_lines_1 = len(all_lines_1)
    number_lines_2 = len(all_lines_2)    
    
    if number_lines_1 == number_lines_2:
       number_lines = number_lines_1  
    else:
       pass_number = 0
       number_lines = 0
       print "Error, two files have different number of lines",number_lines_1,number_lines_2
    
    # go back to start of file
    filein_1.seek(0)
    filein_2.seek(0)
    
    # Read in files line by line ======================================
    # start with zero length lists then append in results on the go
    value_1 = []     
    value_2 = [] 
    
    # python list starts from 0 --> [0,1,2,3,...]
    # then range(1,5) = [1,2,3,4] so dont include -1 at end      
    for i in range(0,number_lines):                       

       line_1 = filein_1.readline()
       line_2 = filein_2.readline()
       
       line_split_1 = line_1.split()
       line_split_2 = line_2.split()
    
       value_1.append( float(line_split_1[1]) )    
       value_2.append( float(line_split_2[1]) )
        
    # close the file
    filein_1.close()
    filein_2.close()
        
    # Pass or Fail Test ===============================================
    # Now compare 
    # pass_number = 1 --> pass
    # pass_number = 0 --> fail
    
    difference1 = end1 - start1
    difference2 = end2 - start2
        
    if difference1 == difference2:
       difference = difference1  
    else:
       pass_number = 0
       difference = 0
       print "Error, two seperate ranges for the files asked for"
    
    maxdiff = 0.0
    # python list starts from 0 --> [0,1,2,3,...]
    # then range(1,5) = [1,2,3,4] so dont include -1 at end  
    for i in range(0,difference):
        
       j = i - 1 + start1
       k = i - 1 + start2
       
       if abs(value_2[j]) > 0.0:
          maxdiff = max(maxdiff,abs(value_1[j] - value_2[k])/abs(value_2[j]))

       if abs(value_1[j] - value_2[k]) > rel_tolerance*abs(value_2[j]):   
          pass_number = 0
          print "*****Failed on value:",j+1,k+1
          print "     value_1,value_2:",value_1[j],value_2[k]
    
    print "Relative tolerance, Max relative difference",rel_tolerance,maxdiff
    print "Finish Compare two file: pass_number(1=pass,0=fail)",pass_number            
                     
    return pass_number

# --------------------------------------------------------------------------------------------------------------------

def divide_variables(numerator,denominator):
   # divide two numbers, return value
   
   value = numerator / denominator
   
   return value


# --------------------------------------------------------------------------------------------------------------------

def compare_variables(reference, current, error, zerotol=1.0e-14):
    # This takes in an array for a particular variable
    # (e.g. kinetic energy) containing the values of that
    # variable at the timesteps. It compares the output
    # of the current run against a reference run,
    # checking that the relative error is within the bound
    # specified.    

    assert len(reference) == len(current)
    relerrs = []

    for i in range(len(current)):
       if abs(reference[i]) > zerotol: # decide if reference[i] is "0.0" or not
           diff = abs(reference[i] - current[i])
           relerr = (diff / reference[i])
           relerrs.append(relerr)
       else:
           relerrs.append(abs(current[i])) # not really a relative error but however

    maxerr = max(relerrs)
    print "Asserting max relative error is smaller than", error
    print "max relative error: %s; index: %s" % (maxerr, relerrs.index(maxerr))
    assert maxerr < error

# --------------------------------------------------------------------------------------------------------------------

def compare_variable(reference, current, error, zerotol=1.0e-14):
    # Compares current value with reference value. Relative error
    # should be smaller than 'error'. If the reference is within
    # zerotol however, an absolute error will be used.

    diff = abs(reference - current)
    if abs(reference) > zerotol: # decide if reference is "0.0" or not
       relerr = (diff / reference)
       print "Asserting relative error is smaller than", error
       print "relative error: " + `relerr`
       assert relerr < error
    else:
       print "Asserting absolute error is smaller than", error
       print "absolute error: "+ `diff`
       assert diff < error

# --------------------------------------------------------------------------------------------------------------------

def diff_squared(value1,value2):
    # find the square of the difference between two values and return answer
    diff = value1 - value2
    answer = diff*diff
    
    return answer

# --------------------------------------------------------------------------------------------------------------------

def squared(value):
    # find the square of the value and return answer
    answer = value*value
    
    return answer

# --------------------------------------------------------------------------------------------------------------------

def diff_absolute(value1,value2):
    # find the absolute difference between two values and return answer
    answer = abs(value1 - value2)
    
    return answer

# --------------------------------------------------------------------------------------------------------------------

def diff_relative(value1,value2):
    # find the relative difference between two values and return answer
    
    if abs(value1) > 0.0:
       
       answer = abs(value1 - value2) / abs(value1)
    
    else:
       
       answer = 0.0
    
    return answer

# --------------------------------------------------------------------------------------------------------------------

def add(value1,value2):
    # add two values together and return the answer
    answer = value1 + value2
    
    return answer

# --------------------------------------------------------------------------------------------------------------------

def add_absolute(value1,value2):
    # add two absolute values together and return the answer
    answer = abs(value1) + abs(value2)
    
    return answer

# --------------------------------------------------------------------------------------------------------------------

def sqrt(value):
    # find the sqrt of value and return answer
    answer = value**(0.5)
    
    return answer

# --------------------------------------------------------------------------------------------------------------------

def times(value1,value2):
    # times value1 and value2 together and return the answer
    answer = value1*value2

    return answer

# --------------------------------------------------------------------------------------------------------------------

def max_value(value1,value2):
    # find max of value1 and value2 and return the answer
    answer = max(value1,value2)

    return answer

# --------------------------------------------------------------------------------------------------------------------

def print_value(value,value_name):
    # print the value name and value
    
    value_name1 = str(value_name)
    
    print value_name1,value
    
    return 

# --------------------------------------------------------------------------------------------------------------------
