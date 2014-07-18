  interface addref
     module procedure addref_tensor_field
  end interface

  interface incref
     module procedure incref_tensor_field
  end interface

  interface decref
     module procedure decref_tensor_field
  end interface  

  interface has_references
     module procedure has_references_tensor_field
  end interface
 
