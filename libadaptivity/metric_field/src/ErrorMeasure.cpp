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

#include "ErrorMeasure.h"

using namespace std;

ErrorMeasure::ErrorMeasure(){
  verbose_off();
  
 ug = NULL;
}

ErrorMeasure::~ErrorMeasure(){
}

void ErrorMeasure::add_field(string field, double error, bool relative, double sigma){
  if(verbose)
    cout<<"void ErrorMeasure::add_field(...)\n";
  
  if(ug==NULL){
    cerr<<"ERROR ("<<__FILE__<<", "<<__LINE__<<"): no mesh has been provided\n";
  }
  
  assert(ug->GetPointData()->GetArray(field.c_str())!=NULL);
  get_hessian(field);
  string hessian_name(field+"_Hessian");

  double bbox[9];
  ug->GetBounds(bbox);
  double dbbox[dim];
  for(size_t i=0;i<dim;i++)
    dbbox[i] = bbox[i*2+1]-bbox[i*2];
  
  if(ug->GetPointData()->GetArray("metric")==NULL){
    vtkDoubleArray *array = vtkDoubleArray::New();
    array->SetName("metric");
    array->SetNumberOfComponents(dim*dim);
    array->SetNumberOfTuples(ug->GetNumberOfPoints());
    
    for(int i=0;i<ug->GetNumberOfPoints();i++){
      if(dim==2){
        array->SetTuple9(i,
                         1/(dbbox[0]*dbbox[0]), 0.0,                   0.0,
                         0.0,                   1/(dbbox[1]*dbbox[1]), 0.0,
                         0.0,                   0.0,                   1.0);
        
      }else if(dim==3){
        array->SetTuple9(i,
                         1/(dbbox[0]*dbbox[0]), 0.0,                   0.0,
                         0.0,                   1/(dbbox[1]*dbbox[1]), 0.0,
                         0.0,                   0.0,                   1/(dbbox[2]*dbbox[2]));
      }else{
        cerr<<"ERROR: unexpected dimension = "<<dim<<endl;
        exit(-1);
      }
    }
    ug->GetPointData()->AddArray(array);
    array->Delete();
  }
  vtkDataArray *m = ug->GetPointData()->GetArray("metric");

  // Convert the Hessian to the required metric tensor
  vector<double> H(dim*dim);
  vtkDataArray *hessian_array = ug->GetPointData()->GetArray(hessian_name.c_str());
  if(hessian_array==NULL){
    cerr<<"ERROR ("<<__FILE__<<", "<<__LINE__<<"): Hessian (\""<<hessian_name<<"\") is missing!\n"
        <<"Data arrays are:\n";
    for(int i=0;i<ug->GetPointData()->GetNumberOfArrays();i++){
      cerr<<ug->GetPointData()->GetArrayName(i)<<endl;
    }
    
    exit(-1);
  }
  
  for(int i=0;i<ug->GetNumberOfPoints();i++){
    hessian_array->GetTuple(i, &(H[0]));
    MetricTensor metric(dim, &(H[0]));
    
    if(!metric.is_null()){
      // Make relative if required
      double local_error;
      if(relative){
        double val = ug->GetPointData()->GetArray(field.c_str())->GetTuple1(i);
        local_error = error*max(sigma, fabs(val));
      }else{
        local_error = error;
      }
      
      metric.scale(1.0/local_error);
      
      // Superimpose metrics
      m->GetTuple(i, &(H[0]));
      
      MetricTensor m1(metric), m2(dim, &(H[0]));

      m1.superimpose(m2);
      if(m1.is_null()){
        m1 = metric;
        m2.superimpose(m1);
        metric = m2;
      }else{
        metric = m1;
      }
         
      // Copy back metric
      metric.get_metric2(&(H[0]));
      m->SetTuple(i, &(H[0]));
    }
  }
  
  ug->GetPointData()->RemoveArray(hessian_name.c_str());;
  return;
}

void ErrorMeasure::apply_gradation(double gradation){
  if(verbose)
    cout<<"void ErrorMeasure::apply_gradation()\n";
  
  // Form NNlist.
  deque< set<size_t> > NNList(ug->GetNumberOfPoints());
  for(int e=0;e<ug->GetNumberOfCells();e++){
    if(ug->GetCell(e)->GetCellType()==VTK_TRIANGLE){
      vtkTriangle *tri = (vtkTriangle *)ug->GetCell(e);
      for(size_t i=0;i<3;i++){
        for(size_t j=i+1;j<3;j++){
          NNList[tri->GetPointId(i)].insert(tri->GetPointId(j));
          NNList[tri->GetPointId(j)].insert(tri->GetPointId(i));
        }
      }
    }else if(ug->GetCell(0)->GetCellType()==VTK_TETRA){
      vtkTetra *tetra = (vtkTetra *)ug->GetCell(e);
      for(size_t i=0;i<4;i++){
        for(size_t j=i+1;j<4;j++){
          NNList[tetra->GetPointId(i)].insert(tetra->GetPointId(j));
          NNList[tetra->GetPointId(j)].insert(tetra->GetPointId(i));
        }
      }
    }else{
      cerr<<"ERROR: unsupported cell type "
          <<ug->GetCell(0)->GetCellType()<<endl;
    }
  }
  
  vtkDataArray *m = ug->GetPointData()->GetArray("metric");
  double log_gradation = log(gradation);
  diagnostics();

  // This is used to ensure we don't revisit parts of the mesh that
  // are known to have converged.
  set<int> hits;
  for(size_t cnt=0; cnt<10; cnt++){
    multimap<double, size_t> ordered_edges;
    if(cnt==0){
      // Iterate over everything.
      for(int n=0;n<ug->GetNumberOfPoints();n++){
        double l = ug->GetPointData()->GetArray("mean_desired_lengths")->GetTuple1(n);
        ordered_edges.insert(pair<double, size_t>(l, n));
      }
    }else{
      // Iterate over only nodes which were modified in the previous iteration.
      for(set<int>::const_iterator n=hits.begin();n!=hits.end();++n){
        double l = ug->GetPointData()->GetArray("mean_desired_lengths")->GetTuple1(*n);
        ordered_edges.insert(pair<double, size_t>(l, *n));
      }
      hits.clear();
    }

    for(multimap<double, size_t>::const_iterator n=ordered_edges.begin();n!=ordered_edges.end(); n++){
      // Used to ensure that the front cannot go back on itself.
      set<size_t> swept;
      
      // Start the new front
      set<size_t> front;
      front.insert(n->second);
      
      while(!front.empty()){
        size_t p=*(front.begin());
        front.erase(p);
        swept.insert(p);
        
        vector<double> Tp(dim*dim);
        vector<double> Dp(dim), Vp(dim*dim);
        
        vector<double> Tq(dim*dim);
        vector<double> Dq(dim), Vq(dim*dim);
        
        for(set<size_t>::const_iterator it=NNList[p].begin(); it!=NNList[p].end();it++){
          size_t q=*it;

          if(swept.count(q))
            continue;
          else
            swept.insert(q);
          
          m->GetTuple(p, &(Tp[0]));
          MetricTensor Mp(dim, &(Tp[0]));        
          Mp.eigen_decomp(&(Dp[0]), &(Vp[0]));
          
          m->GetTuple(q, &(Tq[0]));
          MetricTensor Mq(dim, &(Tq[0]));      
          Mq.eigen_decomp(&(Dq[0]), &(Vq[0]));
          
          // Pair the eigenvectors between p and q by minimising the angle between them.
          vector<int> pairs(dim, -1);
          vector<bool> paired(dim, false);
          for(size_t d=0;d<dim;d++){
            vector<double> angle(dim);
            for(size_t k=0;k<dim;k++){
              if(paired[k])
                continue;
              angle[k] = Vp[d*dim]*Vq[k*dim];
              for(size_t l=1;l<dim;l++)
                angle[k] += Vp[d*dim+l]*Vq[k*dim+l];
              angle[k] = acos(fabs(angle[k]));
            }
            
            size_t r=0;
            for(;r<dim;r++){
              if(!paired[r]){
                pairs[d] = r;
                break;
              }
            }
            r++;
            
            for(;r<dim;r++){
              if(angle[pairs[d]]<angle[r]){
                pairs[d] = r;
              }
            }
            
            paired[pairs[d]] = true;
            
            assert(pairs[d]!=-1);
          }
          
          // Resize eigenvalues if necessary
          double Lpq=length(p, q);
          double dh=Lpq*log_gradation;
          bool add_p=false, add_q=false;
          for(size_t k=0;k<dim;k++){
            double hp = 1.0/sqrt(Dp[k]);
            double hq = 1.0/sqrt(Dq[pairs[k]]);
            double gamma = exp(fabs(hp - hq)/Lpq);
            
            if(isinf(gamma))
              gamma = DBL_MAX;
            if(gamma>(1.05*gradation)){
              if(hp>hq){
                hp = hq + dh;
                Dp[k] = 1.0/(hp*hp);
                add_p = true;
              }else{
                hq = hp + dh;
                Dq[pairs[k]] = 1.0/(hq*hq);
                add_q = true;
              }
            }
          }

          // Reform metrics if modified
          if(add_p){
            front.insert(p);
            
            Mp.eigen_undecomp(&(Dp[0]), &(Vp[0]));
            Mp.get_metric2(&(Tp[0]));
            m->SetTuple(p, &(Tp[0]));
            hits.insert(p);
          }
          if(add_q){
            front.insert(q);
            
            Mq.eigen_undecomp(&(Dq[0]), &(Vq[0]));
            Mq.get_metric2(&(Tq[0]));
            m->SetTuple(q, &(Tq[0]));
            hits.insert(p);
          }
        }
      }
    }
    
    // Refresh diagnostics since we have changed the metric.
    diagnostics();

    if(hits.empty())
      break;
  }

  return;
}

int ErrorMeasure::blas_spev(char jobz, char uplo, int N, const double ap[],
                            double eigenvalues[],  double eigenvectors[],  int ldz,  double work[]){
  int info;
  dspev_(&jobz, &uplo, &N, ap, eigenvalues, eigenvectors, &ldz, work, &info);
  return info;
}

int ErrorMeasure::blas_sgemm(char TRANSA, char TRANSB, int N, double A[], double B[], double C[]){
  double ALPHA=1.0;
  double BETA=1.0;
  dgemm_(&TRANSA, &TRANSB, &N, &N, &N, &ALPHA, A, &N, B, &N, &BETA, C, &N);
  return 0;
}
 
 
/*
  int ErrorMeasure::get_expected_nelements(){
  if(verbose)
  cout<<"int ErrorMeasure::get_expected_nelements()\n";
  // foofar - fix
  return ExpectedNElements;
  }
*/

void ErrorMeasure::diagnostics(){
  if(verbose)
    cout<<"void ErrorMeasure::diagnostics()\n";

  vtkDataArray *m = ug->GetPointData()->GetArray("metric");
  if(m==NULL)
    return;
  
  if(ug->GetPointData()->GetArray("mean_desired_lengths")==NULL){
    vtkDataArray *array = vtkDoubleArray::New();
    array->SetName("mean_desired_lengths");
    array->SetNumberOfComponents(1);
    array->SetNumberOfTuples(ug->GetNumberOfPoints());
    ug->GetPointData()->AddArray(array);
    array->Delete();
  }
  vtkDataArray *mean_lengths = ug->GetPointData()->GetArray("mean_desired_lengths");

  if(ug->GetPointData()->GetArray("desired_lengths")==NULL){
    vtkDataArray *array = vtkDoubleArray::New();
    array->SetName("desired_lengths");
    array->SetNumberOfComponents(9);
    array->SetNumberOfTuples(ug->GetNumberOfPoints());
    ug->GetPointData()->AddArray(array);
    array->Delete();
  }
  vtkDataArray *lengths = ug->GetPointData()->GetArray("desired_lengths");

  double H[dim*dim], D[dim], V[dim*dim];
  for(int i=0;i<ug->GetNumberOfPoints();i++){
    m->GetTuple(i, H);
    MetricTensor metric(dim, H);
    mean_lengths->SetTuple1(i, metric.average_length());
    
    metric.eigen_decomp(D, V);
    for(size_t j=0;j<dim;j++)
      D[j] = 1.0/sqrt(D[j]);
    metric.eigen_undecomp(D, V);
    metric.get_metric2(H);
    lengths->SetTuple(i, H);
  }

  return;
}

void ErrorMeasure::get_hessian(string field){
  if(verbose)
    cout<<"void ErrorMeasure::get_hessian("<<field<<")\n";
  
  // First derivative
  grad(field);
  
  // Second derivatives
  grad(field+"_dx");
  grad(field+"_dy");
  grad(field+"_dz");

  // Put together the hessian
  vtkDoubleArray *H = vtkDoubleArray::New();
  H->SetName(string(field+"_Hessian").c_str());
  H->SetNumberOfComponents(dim*dim);
  H->SetNumberOfTuples(ug->GetNumberOfPoints());
  
  for(int i=0;i<ug->GetNumberOfPoints();i++){
    double ddxx = ug->GetPointData()->GetArray(string(field+"_dx_dx").c_str())->GetTuple1(i);
    double ddxy = ug->GetPointData()->GetArray(string(field+"_dx_dy").c_str())->GetTuple1(i);
    double ddxz = ug->GetPointData()->GetArray(string(field+"_dx_dz").c_str())->GetTuple1(i);

    double ddyx = ug->GetPointData()->GetArray(string(field+"_dy_dx").c_str())->GetTuple1(i);
    double ddyy = ug->GetPointData()->GetArray(string(field+"_dy_dy").c_str())->GetTuple1(i);
    double ddyz = ug->GetPointData()->GetArray(string(field+"_dy_dz").c_str())->GetTuple1(i);

    double ddzx = ug->GetPointData()->GetArray(string(field+"_dz_dx").c_str())->GetTuple1(i);
    double ddzy = ug->GetPointData()->GetArray(string(field+"_dz_dy").c_str())->GetTuple1(i);
    double ddzz = ug->GetPointData()->GetArray(string(field+"_dz_dz").c_str())->GetTuple1(i);

    H->InsertTuple9(i,
                    ddxx,          (ddxy+ddyx)/2, (ddxz+ddzx)/2,
                    (ddxy+ddyx)/2, ddyy,          (ddyz+ddzy)/2,
                    (ddxz+ddzx)/2, (ddyz+ddzy)/2, ddzz        );
  }
  ug->GetPointData()->AddArray(H);
  ug->Update();
  H->Delete();

  ug->GetPointData()->RemoveArray(string(field+"_dx").c_str());
  ug->GetPointData()->RemoveArray(string(field+"_dy").c_str());
  ug->GetPointData()->RemoveArray(string(field+"_dz").c_str());

  ug->GetPointData()->RemoveArray(string(field+"_dx_dx").c_str());
  ug->GetPointData()->RemoveArray(string(field+"_dx_dy").c_str());
  ug->GetPointData()->RemoveArray(string(field+"_dx_dz").c_str());
  ug->GetPointData()->RemoveArray(string(field+"_dy_dx").c_str());
  ug->GetPointData()->RemoveArray(string(field+"_dy_dy").c_str());
  ug->GetPointData()->RemoveArray(string(field+"_dy_dz").c_str());
  ug->GetPointData()->RemoveArray(string(field+"_dz_dx").c_str());
  ug->GetPointData()->RemoveArray(string(field+"_dz_dy").c_str());
  ug->GetPointData()->RemoveArray(string(field+"_dz_dz").c_str());
  ug->Update();

  return;
}

inline double ErrorMeasure::length(size_t n0, size_t n1) const{
  assert(n0!=n1);

  double r0[3], r1[3];
  ug->GetPoints()->GetPoint(n0, r0);
  ug->GetPoints()->GetPoint(n1, r1);
  
  double l = 0.0;
  for(size_t i=0;i<dim;i++)
    l += ((r0[i]-r1[i])*(r0[i]-r1[i]));

  return sqrt(l);
}

int ErrorMeasure::grad(std::string field){
  if(verbose)
    cout<<"vtkUnstructuredGrid *ErrorMeasure::grad("<<field<<")\n";

  ug->GetPointData()->SetActiveScalars(field.c_str());
  
  vtkCellDerivatives *derivatives = vtkCellDerivatives::New();
  derivatives->SetInput(ug);
  derivatives->Update();

  vtkCellDataToPointData *cell2point = vtkCellDataToPointData::New();
  cell2point->SetInput(derivatives->GetUnstructuredGridOutput());
  cell2point->PassCellDataOff();
  cell2point->Update();

  vtkExtractVectorComponents *components = vtkExtractVectorComponents::New();
  components->ExtractToFieldDataOn();
  components->SetInput(cell2point->GetUnstructuredGridOutput());
  components->Update();

  vtkUnstructuredGrid *dug = components->GetUnstructuredGridOutput();

  if(dug->GetPointData()->GetArray("Vorticity")!=NULL){
    dug->GetPointData()->RemoveArray("Vorticity");
    dug->GetPointData()->GetArray("Vorticity-x")->SetName(string(field+"_dx").c_str());
    dug->GetPointData()->GetArray("Vorticity-y")->SetName(string(field+"_dy").c_str());
    dug->GetPointData()->GetArray("Vorticity-z")->SetName(string(field+"_dz").c_str());
  }else{ // Assume this is a newer version of VTK where "Vorticity" is now "ScalarGradient".
    dug->GetPointData()->RemoveArray("ScalarGradient");
    dug->GetPointData()->GetArray("ScalarGradient-x")->SetName(string(field+"_dx").c_str());
    dug->GetPointData()->GetArray("ScalarGradient-y")->SetName(string(field+"_dy").c_str());
    dug->GetPointData()->GetArray("ScalarGradient-z")->SetName(string(field+"_dz").c_str());
  }

  ug->GetPointData()->AddArray(dug->GetPointData()->GetArray(string(field+"_dx").c_str()));
  ug->GetPointData()->AddArray(dug->GetPointData()->GetArray(string(field+"_dy").c_str()));
  ug->GetPointData()->AddArray(dug->GetPointData()->GetArray(string(field+"_dz").c_str()));
  
  derivatives->Delete();
  cell2point->Delete();
  components->Delete();
  
  return 0;
}

void ErrorMeasure::set_max_length(double _max_len){
  if(verbose)
    cout<<"void ErrorMeasure::set_max_length("<<_max_len<<")\n";
  
  double max_len[dim*dim];
  for(size_t i=0;i<dim;i++){
    for(size_t j=0;j<dim;j++){
      if(i==j){
        max_len[i*dim+j] = _max_len;
      }else{
        max_len[i*dim+j] = 0.0;
      }
    }
  }
  
  set_max_length(max_len, 1);

  return;
}
 
void ErrorMeasure::set_max_length(double *max_len, int len){
  if(verbose)
    cout<<"void ErrorMeasure::set_max_length(double *max_len, int len)\n";
  
  if(ug->GetPointData()->GetArray("metric")==NULL)
    return;
  
  assert((len==1)||(len==ug->GetNumberOfPoints()));
  
  vtkDataArray *m = ug->GetPointData()->GetArray("metric");
  
  double H[dim*dim];
  for(int i=0;i<ug->GetNumberOfPoints();i++){
    m->GetTuple(i, H);
    MetricTensor metric(dim, H);
    
    // Limit maximum edge length
    if(len==1)
      metric.limit_max_size(max_len);
    else
      metric.limit_max_size(max_len+dim*dim*i);

    metric.get_metric2(H);
    m->SetTuple(i, H);
  }

  return;
}

void ErrorMeasure::set_max_nodes(int max_nodes){
  if(verbose)
    cout<<"void ErrorMeasure::set_max_nodes("<<max_nodes<<")\n";

  if(ug->GetPointData()->GetArray("metric")==NULL)
    return;
  
  vtkDataArray *m = ug->GetPointData()->GetArray("metric");
  
  int NNodes = ug->GetNumberOfPoints();
  int NElements = ug->GetNumberOfCells();
  vector<double> Metric(9*NNodes), X(NNodes), Y(NNodes), Z(NNodes);
  for(int i=0; i<NNodes; i++){
    double r[3];
    ug->GetPoints()->GetPoint(i, r);
    X[i] = r[0];
    Y[i] = r[1];
    Z[i] = r[2];
    
    m->GetTuple(i, &(Metric[i*9]));
  }
  
  vector<int> ENList(NElements*4);
  for(int i=0;i<NElements;i++){
    assert(ug->GetCell(i)->GetCellType()==VTK_TETRA);
    vtkTetra *tetra = (vtkTetra *)ug->GetCell(i);
    
    for(int j=3;j>=0;j--){
      ENList[i*4+j] = tetra->GetPointId(j) + 1;
    }
  }
  
  int nloc = 4;

  // expected_elements could be a big number so we'll use a double to
  // avoid an integer overflow.
  double expected_elements = get_predicted_nelements_fc(&(Metric[0]), &(X[0]), &(Y[0]), &(Z[0]),
                                                        &(ENList[0]), &NNodes, &NElements, &nloc);  
  double eles_per_node = (double)NElements/NNodes;
  double max_elements = eles_per_node * max_nodes;
  
  if(expected_elements>max_elements){
    double scale = pow(((1.0 / expected_elements) * max_elements), 2.0/dim);

    for(int i=0;i<NNodes;i++){
      for(int j=0;j<9;j++)
        Metric[i*9+j] *= scale;
      
      m->SetTuple(i, &(Metric[i*9]));
    }
  }

  expected_elements = get_predicted_nelements_fc(&(Metric[0]), &(X[0]), &(Y[0]), &(Z[0]),
                                                 &(ENList[0]), &NNodes, &NElements, &nloc);
  return;
}

void ErrorMeasure::set_min_length(double _min_len){
  if(verbose)
    cout<<"void ErrorMeasure::set_min_length("<<_min_len<<")\n";
  
  double min_len[dim*dim];
  for(size_t i=0;i<dim;i++){
    for(size_t j=0;j<dim;j++){
      if(i==j){
        min_len[i*dim+j] = _min_len;
      }else{
        min_len[i*dim+j] = 0.0;
      }
    }
  }
  
  set_min_length(min_len, 1);
  
  return;
}
 
void ErrorMeasure::set_min_length(double *min_len, int len){
  if(verbose)
    cout<<"void ErrorMeasure::set_min_length(double *min_len, int len)\n";
  
  if(ug->GetPointData()->GetArray("metric")==NULL)
    return;
  
  assert((len==1)||(len==ug->GetNumberOfPoints()));
  
  vtkDataArray *m = ug->GetPointData()->GetArray("metric");
  
  double H[dim*dim];
  for(int i=0;i<ug->GetNumberOfPoints();i++){
    m->GetTuple(i, H);
    MetricTensor metric(dim, H);
    
    // Limit minimum edge length
    if(len==1)
      metric.limit_min_size(min_len);
    else
      metric.limit_min_size(min_len+dim*dim*i);
    
    metric.get_metric2(H);
    m->SetTuple(i, H);
  }

  return;
}

void ErrorMeasure::set_input(vtkUnstructuredGrid *_ug){
  if(verbose)
    cout<<"void ErrorMeasure::set_input( ... )\n";
  ug = _ug;

  if(ug->GetCell(0)->GetCellType()==VTK_TRIANGLE){
    dim = 2;
  }else if(ug->GetCell(0)->GetCellType()==VTK_TETRA){
    dim = 3;
  }else{
    cerr<<"ERROR: unsupported cell type "
        <<ug->GetCell(0)->GetCellType()<<endl;
  }

  return;
}

void ErrorMeasure::verbose_off(){
  verbose=false;
}

void ErrorMeasure::verbose_on(){
  verbose=true;
}
