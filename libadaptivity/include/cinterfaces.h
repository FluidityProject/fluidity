/* Copyright (C) 2006 Imperial College London.

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

#ifndef CINTERFACES_H
#define CINTERFACES_H

// Methods
extern "C" {
#define adptvy_fc F77_FUNC(adptvy, ADPTVY) 
  void adptvy_fc(int intBuffer[], const int *intBuffer_size, void *floatBuffer, const int *floatBuffer_size,
                 const int *Geom3D, const int *SRFGMY, const int *USEQ,
                 const int *NNodes, const int *NElements, const int *NSElements, const int *mxnods,
                 const int *sizeENList, const int ENLBasePtr[], const int ENList[], const int ElementRegion[],
                 const int *CLCGMY, const int *sizeSENList, const int SENLBasePtr[], const int SENList[], const int surfID[],
                 const int PRDNDS[], const int *NPRDND,
                 const void *X, const void *Y, const void *Z,
                 const int *ref_NNodes, const int *ref_NElements, const int *ref_sizeENList, const int *ref_ENList, const int *ref_ENLBasePtr,
                 const void *ref_X, const void *ref_Y, const void *ref_Z,
                 const void *Metric, const void *fields, const int *NFREE, const int *TOTFRE, const int *NFIELD,
                 const int *XPCTEL, 
                 int *newNNodes, int *newNElements, int *newNSElements,
                 int *NWSZEN, int *NWSZSN, int *NWSZNN, int *NWNDLC, int *NWSROW,
                 int *NWENLB, int *NWENLS, int *NWSNLB, int *NWSNLS, int *NWSFID,
                 int *NWELRG, int *NWNODX, int *NWNODY, int *NWNODZ,
                 int *NEWMTX, int *NEWFLD,
                 int *ADPBIG, int *ADPNOD,
                 const void *DOTOP, const void *minChance, const int *AdaptIter, const int AdaptOpts[], const int *TWOSTG, const int *TOGTHR,
                 int Gather[], int Scatter[], int *NGather, int *NScatter, int *NPrivateNodes,
                 int ATOSEN[], int ATOREC[], const int *NProcs, const int *debug_level, const int *dbg, const int *chcnsy);
  
#define get_predicted_nelements_fc F77_FUNC(get_predicted_nelements, GET_EXPECTED_NELEMENTS)
  int get_predicted_nelements_fc(const void *, const void *, const void *, const void *, 
                                 const int *, const int *, const int *, const int *);
  
#define adaptmem_fc F77_FUNC(adaptmem, ADAPTMEM)
  int adaptmem_fc(const int *, const int *, const int *, const int *, const int *, 
                  const int *, const int *, const int *, const int *,const  int *, const int *);
  
#define adapt_false_fc F77_FUNC(adapt_false, ADAPT_FALSE)
  int adapt_false_fc();

#define adapt_true_fc F77_FUNC(adapt_true, ADAPT_TRUE)
  int adapt_true_fc();
}

#define LOGICAL_FALSE (adapt_false_fc())
#define LOGICAL_TRUE  (adapt_true_fc())

#endif
