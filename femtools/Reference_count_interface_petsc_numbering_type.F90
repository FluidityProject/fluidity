  interface addref
     module procedure addref_petsc_numbering_type
  end interface

  interface incref
     module procedure incref_petsc_numbering_type
  end interface

  interface decref
     module procedure decref_petsc_numbering_type
  end interface  

  interface has_references
     module procedure has_references_petsc_numbering_type
  end interface
 
