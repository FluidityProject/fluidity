/* Copyright (C) 2004-2006 by Gerard Gorman
   Copyright (C) 2006- Imperial College London and others.
   
   Please see the AUTHORS file in the main source directory for a full list
   of copyright holders.
   
   Dr Gerard J Gorman
   Applied Modelling and Computation Group
   Department of Earth Science and Engineering
   Imperial College London
   
   g.gorman@imperial.ac.uk
   
   This library is free software; you can redistribute it and/or
   modify it under the terms of the GNU Lesser General Public
   License as published by the Free Software Foundation; either
   version 2.1 of the License.
   
   This library is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
   Lesser General Public License for more details.
   
   You should have received a copy of the GNU Lesser General Public
   License along with this library; if not, write to the Free Software
   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307
   USA
*/

extern "C"{
  void vtpopen(char *outName, int *len1, char *vtkTitle, int *len2){}
  void vtpwritepoints(int *NPoints,
		       float *x, float *y, float *z){}
  void vtpwritepointsd(int *NNodes, int *NElems, 
			double *x, double *y, double *z,
			int *enlist, int *elementTypes, int *elementSizes){}
  void vtpstartn(){}
  void vtpwriteghostlevels(int *ghost_levels){}
  void vtpwriteisn(int *vect, char *name, int *len){}
  void vtpwritefsn(float *vect, char *name, int *len){}
  void vtpwritedsn(double *vect, char *name, int *len){}
  void vtpwritefvn(float *vx, float *vy, float *vz,
		      char *name, int *len){}
  void vtpwritedvn(double *vx, double *vy, double *vz,
		      char *name, int *len){}
  void vtpwriteftn(float *v1, float *v2, float *v3, 
          float *v4, float *v5, float *v6, 
          float *v7, float *v8, float *v9, 
          char *name, int *len){}
  void vtpwriteftc(float *v1, float *v2, float *v3, 
          float *v4, float *v5, float *v6, 
          float *v7, float *v8, float *v9, 
          char *name, int *len){}
  
  void vtpwritedtn(double *v1, double *v2, double *v3, 
          double *v4, double *v5, double *v6, 
          double *v7, double *v8, double *v9, 
          char *name, int *len){}
  void vtpwritedtc(double *v1, double *v2, double *v3, 
          double *v4, double *v5, double *v6, 
          double *v7, double *v8, double *v9, 
          char *name, int *len){}
  
  void vtpwriteftn2(float *v1, float *v2, float *v3, 
		       float *v4, float *v5, float *v6, 
		       char *name, int *len){}
  
  void vtpstartc(){}
  void vtpwriteisc(int *vect, char *name, int *len){}
  void vtpwritefsc(float *vect, char *name, int *len){}
  void vtpwritedsc(double *vect, char *name, int *len){}
  void vtpwritefvc(float *vx, float *vy, float *vz,
          char *name, int *len){}
  void vtpwritedvc(double *vx, double *vy, double *vz, 
          char *name, int *len){}
  void vtpclose(){}
  void vtppclose(int *rank, int *npartitions){}
  void vtpsetactivescalars(char* name){}
  void vtpsetactivevectors(char* name){}
  void vtpsetactivetensors(char* name){}

}
