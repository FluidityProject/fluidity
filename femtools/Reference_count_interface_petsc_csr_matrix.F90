  interface addref
     module procedure addref_petsc_csr_matrix
  end interface

  interface incref
     module procedure incref_petsc_csr_matrix
  end interface

  interface decref
     module procedure decref_petsc_csr_matrix
  end interface  

  interface has_references
     module procedure has_references_petsc_csr_matrix
  end interface
 
