  interface addref
     module procedure addref_halo_type
  end interface

  interface incref
     module procedure incref_halo_type
  end interface

  interface decref
     module procedure decref_halo_type
  end interface  

  interface has_references
     module procedure has_references_halo_type
  end interface
 
