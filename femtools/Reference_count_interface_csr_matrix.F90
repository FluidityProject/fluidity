  interface addref
     module procedure addref_csr_matrix
  end interface

  interface incref
     module procedure incref_csr_matrix
  end interface

  interface decref
     module procedure decref_csr_matrix
  end interface  

  interface has_references
     module procedure has_references_csr_matrix
  end interface
 
