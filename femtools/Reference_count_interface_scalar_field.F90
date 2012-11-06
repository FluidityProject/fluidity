  interface addref
     module procedure addref_scalar_field
  end interface

  interface incref
     module procedure incref_scalar_field
  end interface

  interface decref
     module procedure decref_scalar_field
  end interface  

  interface has_references
     module procedure has_references_scalar_field
  end interface
 
