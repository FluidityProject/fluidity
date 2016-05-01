  interface addref
     module procedure addref_quadrature_type
  end interface

  interface incref
     module procedure incref_quadrature_type
  end interface

  interface decref
     module procedure decref_quadrature_type
  end interface  

  interface has_references
     module procedure has_references_quadrature_type
  end interface
 
