  interface addref
     module procedure addref_vector_field
  end interface

  interface incref
     module procedure incref_vector_field
  end interface

  interface decref
     module procedure decref_vector_field
  end interface  

  interface has_references
     module procedure has_references_vector_field
  end interface
 
