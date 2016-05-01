  interface addref
     module procedure addref_dynamic_csr_matrix
  end interface

  interface incref
     module procedure incref_dynamic_csr_matrix
  end interface

  interface decref
     module procedure decref_dynamic_csr_matrix
  end interface  

  interface has_references
     module procedure has_references_dynamic_csr_matrix
  end interface
 
