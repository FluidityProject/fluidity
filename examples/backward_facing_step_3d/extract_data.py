#!/usr/bin/env python

import glob
import sys
import os
import vtktools
import numpy

pwd=os.getenv('PWD')
print(pwd)
mydir=pwd.split('backward_facing_step_3d')[0]+'backward_facing_step_3d'
print(mydir)

def restresseslemoin(off):
  # off is the plot offset

  # get profiles from Le&Moin graphs. x=4
  Le_uu4 = open(str(mydir)+'/Le-profiles/Le-profile1-uu-x4.dat', 'r').readlines()
  rs_uu4 = [off+float(line.split()[0]) for line in Le_uu4]
  y_uu4 = [float(line.split()[1]) for line in Le_uu4]
  Le_vv4 = open(str(mydir)+'/Le-profiles/Le-profile1-vv-x4.dat', 'r').readlines()
  rs_vv4 = [float(line.split()[0]) for line in Le_vv4]
  y_vv4 = [float(line.split()[1]) for line in Le_vv4]
  Le_uv4 = open(str(mydir)+'/Le-profiles/Le-profile1-uv-x4.dat', 'r').readlines()
  rs_uv4 = [float(line.split()[0]) for line in Le_uv4]
  y_uv4 = [float(line.split()[1]) for line in Le_uv4]

  # get profiles from Le&Moin graphs. x=6
  Le_uu6 = open(str(mydir)+'/Le-profiles/Le-profile1-uu-x6.dat', 'r').readlines()
  rs_uu6 = [off+float(line.split()[0]) for line in Le_uu6]
  y_uu6 = [float(line.split()[1]) for line in Le_uu6]
  Le_vv6 = open(str(mydir)+'/Le-profiles/Le-profile1-vv-x6.dat', 'r').readlines()
  rs_vv6 = [float(line.split()[0]) for line in Le_vv6]
  y_vv6 = [float(line.split()[1]) for line in Le_vv6]
  Le_uv6 = open(str(mydir)+'/Le-profiles/Le-profile1-uv-x6.dat', 'r').readlines()
  rs_uv6 = [float(line.split()[0]) for line in Le_uv6]
  y_uv6 = [float(line.split()[1]) for line in Le_uv6]

  # get profiles from Le&Moin graphs. x=10
  Le_uu10 = open(str(mydir)+'/Le-profiles/Le-profile1-uu-x10.dat', 'r').readlines()
  rs_uu10 = [off+float(line.split()[0]) for line in Le_uu10]
  y_uu10 = [float(line.split()[1]) for line in Le_uu10]
  Le_vv10 = open(str(mydir)+'/Le-profiles/Le-profile1-vv-x10.dat', 'r').readlines()
  rs_vv10 = [float(line.split()[0]) for line in Le_vv10]
  y_vv10 = [float(line.split()[1]) for line in Le_vv10]
  Le_uv10 = open(str(mydir)+'/Le-profiles/Le-profile1-uv-x10.dat', 'r').readlines()
  rs_uv10 = [float(line.split()[0]) for line in Le_uv10]
  y_uv10 = [float(line.split()[1]) for line in Le_uv10]

  # get profiles from Le&Moin graphs. x=19
  Le_uu19 = open(str(mydir)+'/Le-profiles/Le-profile1-uu-x19.dat', 'r').readlines()
  rs_uu19 = [off+float(line.split()[0]) for line in Le_uu19]
  y_uu19 = [float(line.split()[1]) for line in Le_uu19]
  Le_vv19 = open(str(mydir)+'/Le-profiles/Le-profile1-vv-x19.dat', 'r').readlines()
  rs_vv19 = [float(line.split()[0]) for line in Le_vv19]
  y_vv19 = [float(line.split()[1]) for line in Le_vv19]
  Le_uv19 = open(str(mydir)+'/Le-profiles/Le-profile1-uv-x19.dat', 'r').readlines()
  rs_uv19 = [float(line.split()[0]) for line in Le_uv19]
  y_uv19 = [float(line.split()[1]) for line in Le_uv19]

  return y_uu4,rs_uu4,y_vv4,rs_vv4,y_uv4,rs_uv4,y_uu6,rs_uu6,y_vv6,rs_vv6,y_uv6,rs_uv6,y_uu10,rs_uu10,y_vv10,rs_vv10,y_uv10,rs_uv10,y_uu19,rs_uu19,y_vv19,rs_vv19,y_uv19,rs_uv19

def velocityprofileslemoin():

  Le = open(str(mydir)+'/Le-profiles/Le-profile1-U-x4.dat', 'r').readlines()
  Le_u4 = [float(line.split()[0]) for line in Le]
  Le_y4 = [float(line.split()[1]) for line in Le]
  jd = open(str(mydir)+'/Le-profiles/JD-profile1-U-x4.dat', 'r').readlines()
  jd_u4 = [float(line.split()[0]) for line in jd]
  jd_y4 = [float(line.split()[1]) for line in jd]
  Le = open(str(mydir)+'/Le-profiles/Le-profile1-U-x6.dat', 'r').readlines()
  Le_u6 = [float(line.split()[0]) for line in Le]
  Le_y6 = [float(line.split()[1]) for line in Le]
  jd = open(str(mydir)+'/Le-profiles/JD-profile1-U-x6.dat', 'r').readlines()
  jd_u6 = [float(line.split()[0]) for line in jd]
  jd_y6 = [float(line.split()[1]) for line in jd]
  Le = open(str(mydir)+'/Le-profiles/Le-profile1-U-x10.dat', 'r').readlines()
  Le_u10 = [float(line.split()[0]) for line in Le]
  Le_y10 = [float(line.split()[1]) for line in Le]
  jd = open(str(mydir)+'/Le-profiles/JD-profile1-U-x10.dat', 'r').readlines()
  jd_u10 = [float(line.split()[0]) for line in jd]
  jd_y10 = [float(line.split()[1]) for line in jd]
  Le = open(str(mydir)+'/Le-profiles/Le-profile1-U-x19.dat', 'r').readlines()
  Le_u19 = [float(line.split()[0]) for line in Le]
  Le_y19 = [float(line.split()[1]) for line in Le]
  jd = open(str(mydir)+'/Le-profiles/JD-profile1-U-x19.dat', 'r').readlines()
  jd_u19 = [float(line.split()[0]) for line in jd]
  jd_y19 = [float(line.split()[1]) for line in jd]

  return Le_y4,Le_u4,jd_y4,jd_u4,Le_y6,Le_u6,jd_y6,jd_u6,Le_y10,Le_u10,jd_y10,jd_u10,Le_y19,Le_u19,jd_y19,jd_u19

def ercoftacvelocityprofiles():

  y4=[];U4=[];y6=[];U6=[];y10=[];U10=[];y19=[];U19=[]
  datafile = open(str(mydir)+'/Ercoftac-test31-BFS/BFS-SEM-ERCOFTAC360-table.dat', 'r')
  for line in range(1): datafile.readline()
  for line in datafile: y4.append(float(line.split()[0])); U4.append(float(line.split()[1]))
  datafile = open(str(mydir)+'/Ercoftac-test31-BFS/BFS-SEM-ERCOFTAC411-table.dat', 'r')
  for line in range(1): datafile.readline()
  for line in datafile: y6.append(float(line.split()[0])); U6.append(float(line.split()[1]))
  datafile = open(str(mydir)+'/Ercoftac-test31-BFS/BFS-SEM-ERCOFTAC513-table.dat', 'r')
  for line in range(1): datafile.readline()
  for line in datafile: y10.append(float(line.split()[0])); U10.append(float(line.split()[1]))
  datafile = open(str(mydir)+'/Ercoftac-test31-BFS/BFS-SEM-ERCOFTAC744-table.dat', 'r')
  for line in range(1): datafile.readline()  
  for line in datafile: y19.append(float(line.split()[0])); U19.append(float(line.split()[1]))

  return y4,U4,y6,U6,y10,U10,y19,U19

def ercoftacrestressprofiles():

  datafile = open(str(mydir)+'/Ercoftac-test31-BFS/BFS-SEM-ERCOFTAC181-table.dat', 'r')
  print("reading in data from file: BFS-SEM-ERCOFTAC181-table.dat")
  # ignore header line
  for line in range(1):
    datafile.readline()
  
  y=[];U=[];uu=[];vv=[];uv=[]
  for line in datafile:
    y.append(1.0+float(line.split()[0]))
    U.append(0.1*float(line.split()[1]))
    uu.append(float(line.split()[2]))
    vv.append(float(line.split()[3]))
    uv.append(float(line.split()[4]))

  return y,U,uu,vv,uv

def panjwaniprofiles(variable):

  off=0.05 if variable=='uu' else 0.0
  inv=-1.0 if variable=='uv' else 1.0
  # get data from relevant files
  pj = open(str(mydir)+'/panjwani-profiles/panjwani-dynles-'+str(variable)+'-4.8.dat', 'r').readlines()
  pj_u1 = [off+inv*float(line.split()[0]) for line in pj]
  pj_y1 = [float(line.split()[1]) for line in pj]
  pj = open(str(mydir)+'/panjwani-profiles/panjwani-dynles-'+str(variable)+'-7.3.dat', 'r').readlines()
  pj_u2 = [off+inv*float(line.split()[0]) for line in pj]
  pj_y2 = [float(line.split()[1]) for line in pj]
  pj = open(str(mydir)+'/panjwani-profiles/panjwani-dynles-'+str(variable)+'-12.2.dat', 'r').readlines()
  pj_u3 = [off+inv*float(line.split()[0]) for line in pj]
  pj_y3 = [float(line.split()[1]) for line in pj]
  pj = open(str(mydir)+'/panjwani-profiles/panjwani-dynles-'+str(variable)+'-18.2.dat', 'r').readlines()
  pj_u4 = [off+inv*float(line.split()[0]) for line in pj]
  pj_y4 = [float(line.split()[1]) for line in pj]

  return pj_u1,pj_y1,pj_u2,pj_y2,pj_u3,pj_y3,pj_u4,pj_y4



