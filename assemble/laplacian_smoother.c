void lap_smoother(int dimension, int num_nodes, int num_elements, int num_surf_elements, int * connectivity, double * phys_mesh, double* smooth_mesh, double * comp_mesh, int * surf_connectivity) {

//The smoothing algorithm
// Just a straight copy for now

  int i;
  for(i=0;i<num_nodes; i++) {
    smooth_mesh[i] = phys_mesh[i];
  }

}
