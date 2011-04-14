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

#include "MetricTensor.h"
using namespace std;

bool MetricTensor::verbose = false;

MetricTensor::MetricTensor(int dimension){
  if(verbose)
    cout<<"MetricTensor::MetricTensor("<<dimension<<")\n";

  dim = dimension;
  if((dim!=2)&&(dim!=3)){
    cerr<<"ERROR ("<<__FILE__<<", "<<__LINE__<<"): Unsupported dimension "<<dim<<endl;
    exit(-1);
  }

  t_size = (dim*dim - (dim*(dim-1))/2);
  metric.resize(t_size);
}

MetricTensor::MetricTensor(double m11,
                           double m12, double m22){
  if(verbose)
    cout<<"MetricTensor::MetricTensor("<<m11<<",\n"
        <<"                           "<<m12<<", "<<m22<<")\n";
  
  dim = 2;
  t_size = 3;
  metric.resize(t_size);
  metric[0] = m11;
  metric[1] = m12;
  metric[2] = m22;
}

MetricTensor::MetricTensor(double m11,
                           double m12, double m22,
                           double m13, double m23, double m33){
  if(verbose)
    cout<<"MetricTensor::MetricTensor("<<m11<<",\n"
        <<"                           "<<m12<<", "<<m22<<",\n"
        <<"                           "<<m13<<", "<<m23<<", "<<m33<<")\n";

  dim = 3;
  t_size = 6;
  metric.resize(t_size);
  metric[0] = m11;
  metric[1] = m12; metric[2] = m22;
  metric[3] = m13; metric[4] = m23; metric[5] = m33;
}

MetricTensor::MetricTensor(int dimension, const double *M){
  if(verbose){
    cout<<"MetricTensor::MetricTensor("<<dimension<<", const double *M)\n";
    
    cout<<"M =\n";
    for(int i=0;i<dimension;i++){
      for(int j=0;j<dimension;j++){
        cout<<M[i*dimension+j]<<",\t";
      }
      cout<<"\n";
    }
  }
  
  dim = dimension;
  
  if(dim==2){
    t_size = 3;
    metric.resize(t_size);
    metric[0] = M[0];
    metric[1] = M[2]; metric[2] = M[3];
  }else if(dim==3){
    t_size = 6;
    metric.resize(t_size);
    metric[0] = M[0];
    metric[1] = M[3]; metric[2] = M[4];
    metric[3] = M[6]; metric[4] = M[7]; metric[5] = M[8];
  }else{
    cerr<<"ERROR: unexpected dimenbsion = "<<dim<<endl;
    exit(-1);
  }
}

MetricTensor::MetricTensor(const MetricTensor& A){
  if(verbose)
    cout<<"MetricTensor::MetricTensor(const MetricTensor& A)\n";
  *this = A;
}

const MetricTensor& MetricTensor::operator=(const MetricTensor &A){
  dim = A.dim;
  t_size = A.t_size;
  metric = A.metric;
  return *this;
}

std::ostream &operator<<(std::ostream& out, const MetricTensor &A){
  for(size_t i=0;i<A.dim;i++){
    for(size_t j=0;j<A.dim;j++){
      out<<A.metric[A.lookup(i,j)]<<"\t";
    }
    out<<"\n";
  }
  return out;
}

double MetricTensor::aspect_ratio() const{
  if(verbose)
    cout<<"double MetricTensor::aspect_ratio() const\n";

  double D[dim], V[dim*dim];
  eigen_decomp(D, V);
  double hmin =1/sqrt(D[0]);
  double hmax =1/sqrt(D[0]);
  for(size_t i=1;i<dim;i++){
    hmin = min(hmin, 1/sqrt(D[i]));
    hmax = max(hmax, 1/sqrt(D[i]));
  }
  
  return hmax/hmin;
}

double MetricTensor::average_length() const{
  if(verbose)
    cout<<"double MetricTensor::average_ev() const\n";

  double D[dim], V[dim*dim];
  eigen_decomp(D, V);
  
  double sum = D[0];
  for(size_t i=1;i<dim;i++)
    sum+=D[i];
  double average = sum/dim;

  return sqrt(1.0/average);
}

int MetricTensor::blas_spev(char jobz, char uplo, int N, const double ap[],  
			     double eigenvalues[],  double eigenvectors[],
			     int ldz,  double work[]) const{
  
  int info;
  dspev_(&jobz, &uplo, &N, ap, eigenvalues, eigenvectors, &ldz, work, &info);
  return info;
}

int MetricTensor::circumscribe(const MetricTensor& M){
  if(verbose)
    cout<<"int MetricTensor::circumscribe(const MetricTensor& M)\n";

  map(M); 
  double D[dim], V[dim*dim];
  eigen_decomp(D, V);
  
  for(size_t i=0;i<dim;i++)
    D[i] = std::max((double)1.0, D[i]);
  
  eigen_undecomp(D, V);
  unmap(M); 
  return 0;
}

double MetricTensor::cofactor(int i, int j) const{
  if(verbose)
    cout<<"double MetricTensor::cofactor(int i, int j) const\n";

  double cofactor;
  if(dim==2){
    int sign=1;
    if((i+j)%2)
      sign = -1;
    
    int ii = (i+1)%2;
    int jj = (j+1)%2;
    
    cofactor = sign*metric[lookup(ii, jj)];
  }else{ // if(dim==3)
    int sign=1;
    if((i+j)%2)
      sign = -1;
    
    double c[4];
    int ii=0;
    for(int k=0;k<3;k++){
      if(i!=k){
        int jj=0;
        for(int l=0;l<3;l++)
          if(j!=l){
            c[ii*2+jj++] = metric[lookup(k, l)];
          }
        ii++;
      }
    }
    
    cofactor = sign*(c[0]*c[3] - c[1]*c[2]);
  }

  return cofactor;
}

double MetricTensor::det() const{
  if(verbose)
    cout<<"double MetricTensor::det() const\n";

  double d=0.0;
  if(dim==2){
    d = metric[0]*metric[2] - metric[1]*metric[1];
  }else if(dim==3){
    d =
      metric[0]*(metric[2]*metric[5] - metric[4]*metric[4]) -
      metric[1]*(metric[1]*metric[5] - metric[3]*metric[4]) +
      metric[3]*(metric[1]*metric[4] - metric[3]*metric[2]);
  }else{
    cerr<<"ERROR: unexpected dimenbsion = "<<dim<<endl;
    exit(-1);
  }

  return d;
}

int MetricTensor::eigen_decomp(double *eigenvalues, double *eigenvectors) const{
  if(verbose)
    cout<<"int MetricTensor::eigen_decomp(double *eigenvalues, double *eigenvectors) const\n";

  char jobz = 'V';
  char uplo = 'U';
  double ap[t_size];
  for(size_t i=0;i<t_size;i++)
    ap[i] = metric[i];
  double work[dim*dim];
  int info = blas_spev(jobz, uplo, dim, ap, eigenvalues, eigenvectors, dim, work);

  if(info<0){
    cerr<<"Failed in eigenvalue decomposition. Argument "<<abs(info)<<" was foobar\n";
  }else if(info>0){
    cerr<<"Failed in eigenvalue decomposition. The algorithm  failed  to  converge; "
        <<info<<" off-diagonal elements of an intermediate tridiagonal form did not "
        <<"converge to zero.\n";
    for(size_t i=0;i<t_size;i++)
      cerr<<"metric = "<<metric[i]<<"\t";
    cerr<<endl;
  }

  for(size_t i=0;i<dim;i++)
    eigenvalues[i] = fabs(eigenvalues[i]);

  return info;
}

int MetricTensor::eigen_decomp(const double *T, double *eigenvalues, double *eigenvectors) const{
  if(verbose)
    cout<<"int MetricTensor::eigen_decomp(const double *T, double *eigenvalues, double *eigenvectors) const\n";

  char jobz = 'V';
  char uplo = 'U';
  double ap[t_size];

  for(size_t i=0;i<dim;i++){
    for(size_t j=0;j<=i;j++){
      ap[lookup(i, j)] = T[i*dim+j];
    }
  }
  double work[dim*dim];
  int info = blas_spev(jobz, uplo, dim, ap, eigenvalues, eigenvectors, dim, work);

  if(info<0){
    cerr<<"Failed in eigenvalue decomposition. Argument "<<abs(info)<<" was foobar\n";
  }else if(info>0){
    cerr<<"Failed in eigenvalue decomposition. The algorithm  failed  to  converge; "
        <<info<<" off-diagonal elements of an intermediate tridiagonal form did not "
        <<"converge to zero.\n";
    for(size_t i=0;i<t_size;i++)
      cerr<<"metric = "<<metric[i]<<"\t";
    cerr<<endl;
  }

  for(size_t i=0;i<dim;i++)
    eigenvalues[i] = fabs(eigenvalues[i]);

  return info;
}

int MetricTensor::eigen_undecomp(const double *D, const double *V){
  if(verbose)
    cout<<"int MetricTensor::eigen_undecomp(const double *D, const double *V)\n";

  // Insure eigenvalues are positive
  double eigenvalues[dim];
  for(size_t i=0;i<dim;i++)
    eigenvalues[i] = fabs(D[i]);
  
  for(size_t i=0;i<dim;i++)
    for(size_t j=0;j<dim;j++){
      int ii = lookup(i,j);
      metric[ii] = 0.0;
      for(size_t k=0;k<dim;k++)
        metric[ii]+=eigenvalues[k]*V[k*dim+i]*V[k*dim+j];
    }

  return 0;
}

int MetricTensor::eigen_undecomp(const double *D, const double *V, double *T) const{
  if(verbose)
    cout<<"int MetricTensor::eigen_undecomp(const double *D, const double *V)\n";

  // Insure eigenvalues are positive
  double eigenvalues[dim];
  for(size_t i=0;i<dim;i++)
    eigenvalues[i] = fabs(D[i]);
  
  for(size_t i=0;i<dim;i++)
    for(size_t j=0;j<dim;j++){
      T[i*dim+j] = 0.0;
      for(size_t k=0;k<dim;k++)
        T[i*dim+j]+=eigenvalues[k]*V[k*dim+i]*V[k*dim+j];
    }

  return 0;
}

MetricTensor MetricTensor::dot(const MetricTensor &M2){
  if(verbose)
    cout<<"MetricTensor MetricTensor::dot(const MetricTensor &M2)\n";

  if(dim==2){
    double m11 = 0.0;
    for(int i=0;i<2;i++)
      m11+=metric[lookup(i,0)]*M2.metric[lookup(0,i)];
    
    double m12 = 0.0;
    for(int i=0;i<2;i++)
      m12+=metric[lookup(i,0)]*M2.metric[lookup(1,i)];
    
    double m22 = 0.0;
    for(int i=0;i<2;i++)
      m22+=metric[lookup(i,1)]*M2.metric[lookup(1,i)];
    
    return MetricTensor(m11,
                        m12, m22);
  }else{ // if(dim==3)
    double m11 = 0.0;
    for(int i=0;i<3;i++)
      m11+=metric[lookup(i,0)]*M2.metric[lookup(0,i)];
    
    double m12 = 0.0;
    for(int i=0;i<3;i++)
      m12+=metric[lookup(i,0)]*M2.metric[lookup(1,i)];
    
    double m13 = 0.0;
    for(int i=0;i<3;i++)
      m13+=metric[lookup(i,0)]*M2.metric[lookup(2,i)];
    
    double m22 = 0.0;
    for(int i=0;i<3;i++)
      m22+=metric[lookup(i,1)]*M2.metric[lookup(1,i)];
    
    double m23 = 0.0;
    for(int i=0;i<3;i++)
      m23+=metric[lookup(i,1)]*M2.metric[lookup(2,i)];
    
    double m33 = 0.0;
    for(int i=0;i<3;i++)
      m33+=metric[lookup(i,2)]*M2.metric[lookup(2,i)];
    
    return MetricTensor(m11,
                        m12, m22,
                        m13, m23, m33); 
  }
}

void MetricTensor::get_metric(double M[]) const{
  if(verbose)
    cout<<"void MetricTensor::get_metric(double M[]) const\n";

  for(size_t i=0;i<t_size;i++)
    M[i] = metric[i];
  return;
}

void MetricTensor::get_metric2(double M[]) const{
  if(verbose)
    cout<<"void MetricTensor::get_metric2(double M[]) const\n";

  for(size_t i=0;i<dim;i++)
    for(size_t j=0;j<dim;j++)
      M[i*dim+j] = metric[lookup(i, j)];
  return;
}

int MetricTensor::inscribe(const MetricTensor& M){
  if(verbose)
    cout<<"int MetricTensor::inscribe(const MetricTensor& M)\n";

  map(M); // Map self to the same space where M is the unit sphere

  double D[dim], V[dim*dim];
  eigen_decomp(D, V);
  
  
  for(size_t i=0;i<dim;i++)
    D[i] = std::min((double)1.0, D[i]);
  
  eigen_undecomp(D, V);
  
  unmap(M);
  return 0;
}

int MetricTensor::limit_max_size(const double *max_len){
  if(verbose)
    cout<<"int MetricTensor::limit_max_size(const double *max_len)\n";

  double D[dim], V[dim*dim];
  eigen_decomp(max_len, D, V);
  for(size_t i=0;i<dim;i++)
    D[i] = 1.0/(D[i]*D[i]);
  
  double T[dim*dim];
  eigen_undecomp(D, V, T);
    
  circumscribe(MetricTensor(dim, T));

  return 0;
}

int MetricTensor::limit_min_size(const double *min_len){
  if(verbose)
    cout<<"int MetricTensor::limit_min_size(const double *min_len)\n";

  double D[dim], V[dim*dim];

  eigen_decomp(min_len, D, V);
  for(size_t i=0;i<dim;i++)
    D[i] = 1.0/(D[i]*D[i]);

  double T[dim*dim];
  eigen_undecomp(D, V, T);

  inscribe(MetricTensor(dim, T));
  
  return 0;
}

MetricTensor MetricTensor::inv() const{
  if(verbose)
    cout<<"MetricTensor MetricTensor::inv() const\n";

  double d = det();
  if(fabs(d)<FLT_SMALL){
    vector<double> components(t_size, 0.0);
    return MetricTensor(dim, &(components[0])); 
  }else{
    vector<double> components(t_size);
    for(size_t i=0, loc=0;i<dim;i++)
      for(size_t j=0;j<=i;j++,loc++)
        components[loc] = cofactor(i, j)/d;
    
    return MetricTensor(dim, &(components[0]));
  }
}

int MetricTensor::inv(double &v00, double &v01,
                      double &v10, double &v11) const{
  if(verbose)
    cout<<"int MetricTensor::inv(double &v00, double &v01,double &v10, double &v11) const\n";

  double d = v00*v11 - v10*v01;
  
  if(fabs(d)<FLT_SMALL)
    return -1;
  
  double _v00 = v00, _v01 = v01;
  double _v10 = v10, _v11 = v11;

  v00 =  _v11/d; v01 = -_v10/d;
  v10 = -_v01/d; v11 =  _v00/d;

  return 0;
}

int MetricTensor::inv(double &v00, double &v01, double &v02,
                      double &v10, double &v11, double &v12,
                      double &v20, double &v21, double &v22) const{
  if(verbose)
    cout<<"int MetricTensor::inv(double &v00, double &v01, double &v02, "
        <<"double &v10, double &v11, double &v12, "
        <<"double &v20, double &v21, double &v22) const\n";

  double d = 
    v00*(v11*v22 - v12*v21) -
    v01*(v10*v22 - v12*v20) +
    v02*(v10*v21 - v11*v20);
  
  if(fabs(d)<FLT_SMALL)
    return -1;
  
  double _v00 = v00, _v01 = v01, _v02 = v02;
  double _v10 = v10, _v11 = v11, _v12 = v12;
  double _v20 = v20, _v21 = v21, _v22 = v22;

  v00 =  (_v11*_v22 - _v12*_v21)/d; v01 = -(_v10*_v22 - _v12*_v20)/d; v02 =  (_v10*_v21 - _v11*_v20)/d;
  v10 = -(_v01*_v22 - _v02*_v21)/d; v11 =  (_v00*_v22 - _v02*_v20)/d; v12 = -(_v00*_v21 - _v01*_v20)/d;
  v20 =  (_v01*_v12 - _v02*_v11)/d; v21 = -(_v00*_v12 - _v02*_v10)/d; v22 =  (_v00*_v11 - _v01*_v10)/d;

  return 0;
}

bool MetricTensor::is_null() const{
  return fabs(det())<FLT_SMALL;
}

int MetricTensor::limit_aspect(double MaximumAspectRatio){
  if(verbose)
    cout<<"int MetricTensor::limit_aspect(double MaximumAspectRatio)\n";

  double D[dim], V[dim*dim];
  eigen_decomp(D, V);

  double max_D = D[0];
  for(size_t i=1;i<dim;i++)
    max_D = std::max(max_D, D[i]);

  double min_D = max_D/(MaximumAspectRatio*MaximumAspectRatio);

  for(size_t i=0;i<dim;i++)
    D[i] = std::max(min_D, D[i]);

  eigen_undecomp(D, V);

  return 0;
}

int MetricTensor::map(const MetricTensor &M){
  if(verbose)
    cout<<"int MetricTensor::map(const MetricTensor &M)\n";

  double D[dim], V[dim*dim];
  M.eigen_decomp(D, V);
  for(size_t i=0;i<dim;i++)
    D[i] = sqrt(D[i]);
  
  double F[dim*dim];
  if(dim==2){
    F[0] = D[0]*V[0]; F[2] = D[0]*V[2]; 
    F[1] = D[1]*V[1]; F[3] = D[1]*V[3];
    
    if(inv(F[0], F[2],
           F[1], F[3]))
      return -1;
  }else{
    F[0] = D[0]*V[0]; F[3] = D[0]*V[1]; F[6] = D[0]*V[2];
    F[1] = D[1]*V[3]; F[4] = D[1]*V[4]; F[7] = D[1]*V[5];
    F[2] = D[2]*V[6]; F[5] = D[2]*V[7]; F[8] = D[2]*V[8];
    
    if(inv(F[0], F[3], F[6],
           F[1], F[4], F[7],
           F[2], F[5], F[8]))
      return -1;
  }
  
  double M1[dim*dim];
  for(size_t i=0;i<dim;i++)
    for(size_t j=0;j<dim;j++)
      M1[i*dim+j] = metric[lookup(i, j)];
  
  double C[dim*dim];
  MatrixDotMatrix(F,false,M1,false,C);
  MatrixDotMatrix(C,false,F,true,M1);
  
  for(size_t i=0;i<dim;i++)
    for(size_t j=0;j<=i;j++)
      metric[lookup(i, j)] = M1[i*dim+j];

  return(0);
}

int MetricTensor::unmap(const MetricTensor &M){
  if(verbose)
    cout<<"int MetricTensor::unmap(const MetricTensor &M)\n";

  double D[dim], V[dim*dim];
  M.eigen_decomp(D, V);
  for(unsigned i=0;i<dim;i++) 
    D[i] = 1.0/sqrt(D[i]);

  double F[dim*dim];
  if(dim==2){
    F[0] = D[0]*V[0]; F[2] = D[0]*V[2]; 
    F[1] = D[1]*V[1]; F[3] = D[1]*V[3];
    
    inv(F[0], F[2],
        F[1], F[3]);
  }else{
    F[0] = D[0]*V[0]; F[3] = D[0]*V[1]; F[6] = D[0]*V[2];
    F[1] = D[1]*V[3]; F[4] = D[1]*V[4]; F[7] = D[1]*V[5];
    F[2] = D[2]*V[6]; F[5] = D[2]*V[7]; F[8] = D[2]*V[8];
    
    inv(F[0], F[3], F[6],
        F[1], F[4], F[7],
        F[2], F[5], F[8]);
  }

  double M1[dim*dim];
  for(size_t i=0;i<dim;i++)
    for(size_t j=0;j<dim;j++)
      M1[i*dim+j] = metric[lookup(i, j)];
  
  double C[dim*dim];
  MatrixDotMatrix(F,true,M1,false,C);
  MatrixDotMatrix(C,false,F,false,M1);
  
  for(size_t i=0;i<dim;i++)
    for(size_t j=0;j<=i;j++)
      metric[lookup(i, j)] = M1[i*dim+j];
  
  return(0);
}

int MetricTensor::MatrixDotMatrix(double *A, bool aT, double *B, bool bT, double *C) const{
  if(verbose)
    cout<<"int MetricTensor::MatrixDotMatrix(double *A, bool aT, double *B, bool bT, double *C) const\n";

  char TRANSA='N';
  if(aT) TRANSA='T';
  
  char TRANSB='N';
  if(bT) TRANSB='T';
  
  int M=dim, N=dim, K=dim;
  int LDA=dim, LDB=dim, LDC=dim;
  
  double ALPHA = 1.0;
  double BETA=0.0;
  
  dgemm_(&TRANSA, &TRANSB, &M, &N, &K, &ALPHA, A, &LDA, B, &LDB, &BETA, C, &LDC);
  
  return 0;
}

void MetricTensor::print_eigen(double *D, double *V) const{
  cout<<"eigenvalues = ";
  for(size_t i=0;i<dim;i++)
    cout<<D[i]<<"\t";
  cout<<endl;
  cout<<"eigenvectors = "<<endl;
  for(size_t i=0;i<dim;i++){
    for(size_t j=0;j<dim;j++)
      cout<<V[i*dim+j]<<"\t";
    cout<<endl;
  }
}

void MetricTensor::scale(double s){
  for(size_t i=0;i<dim;i++)
    metric[i]*=s;
}

int MetricTensor::superimpose(const MetricTensor& M){
  return circumscribe(M);
}


int MetricTensor::lookup(size_t i, size_t j) const{
  assert(i>=0); assert(i<dim);
  assert(j>=0); assert(j<dim);
  if(dim==2){
    switch(i){
    case(0):
      switch(j){
      case(0): return 0;
      case(1): return 1;
      }
    case(1):
      switch(j){
      case(0): return 1;
      case(1): return 2;
      }
    }
  }else{
    switch(i){
    case(0):
      switch(j){
      case(0): return 0;
      case(1): return 1;
      case(2): return 3;
      }
    case(1):
      switch(j){
      case(0): return 1;
      case(1): return 2;
      case(2): return 4;
      }
    case(2): return 3+j;
    } 
  }
  return -1;
}

void MetricTensor::verbose_off(){
  verbose = false;
}

void MetricTensor::verbose_on(){
  verbose = true;
}

int MetricTensor::write_vtk(const char *name) const{
#ifdef HAVE_VTK
  vtkPolyData *dataSet = vtkPolyData::New();
  dataSet->Allocate();
  
  vtkPoints *newPts = vtkPoints::New();
  newPts->InsertNextPoint(0.0, 0.0, 0.0);
  newPts->InsertNextPoint(1.0, 1.0, 1.0);
  dataSet->SetPoints(newPts);
  dataSet->Update();
  newPts->Delete();

  vtkIdType Cell[]={0,1};
  dataSet->InsertNextCell(VTK_POLY_VERTEX, 2, Cell);

  vtkDoubleArray *tensor = vtkDoubleArray::New();
  tensor->SetName("Metric tensor");
  tensor->SetNumberOfComponents(9);
  if(dim==2){
    tensor->InsertNextTuple9(metric[0], metric[1], 0.0,
                             metric[1], metric[2], 0.0,
                             0.0,       0.0,       0.0);
  }else{
    tensor->InsertNextTuple9(metric[0], metric[1], metric[3],
                             metric[1], metric[2], metric[4],
                             metric[3], metric[4], metric[5]);
  }
  tensor->InsertNextTuple9(FLT_MIN, 0.0,     0.0,
			   0.0,     FLT_MIN, 0.0,
			   0.0,     0.0,     FLT_MIN);
  dataSet->GetPointData()->AddArray(tensor);
  dataSet->GetPointData()->SetActiveAttribute("Metric tensor", vtkDataSetAttributes::TENSORS);
  dataSet->Update();
  tensor->Delete();
  
  vtkXMLPolyDataWriter *writer= vtkXMLPolyDataWriter::New();
  vtkZLibDataCompressor* compressor = vtkZLibDataCompressor::New();
  writer->SetCompressor(compressor);
  compressor->Delete();
  
  std::string fname(name);
  fname.append(".vtp");
  writer->SetFileName(fname.c_str());
  writer->SetInput(dataSet);
  writer->Write();
  writer->Delete();
  dataSet->Delete();
#endif
  return 0;
}
