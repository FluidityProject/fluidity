/* Copyright (C) 2006 Imperial College London and others.

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
#ifndef PARMETIS_H
#define PARMETIS_H

typedef int idxtype;

void ParMETIS_V3_AdaptiveRepart(idxtype *vtxdist, idxtype *xadj, idxtype *adjncy,
                                idxtype *vwgt, idxtype *vsize, idxtype *adjwgt, 
                                int *wgtflag, int *numflag, int *ncon, int *nparts,
                                float *tpwgts, float *ubvec, float *itr, int *options,
                                int *edgecut, idxtype *part, MPI_Comm *comm);

void ParMETIS_V3_PartKway(idxtype *, idxtype *, idxtype *, idxtype *, idxtype *,
                          int *, int *, int *, int *, float *, float *,
                          int *, int *, idxtype *, MPI_Comm *);
#endif
