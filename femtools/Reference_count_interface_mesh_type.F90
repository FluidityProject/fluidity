  interface addref
     module procedure addref_mesh_type
  end interface

  interface incref
     module procedure incref_mesh_type
  end interface

  interface decref
     module procedure decref_mesh_type
  end interface  

  interface has_references
     module procedure has_references_mesh_type
  end interface
 
