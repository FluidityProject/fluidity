GFORTRAN module version '6' created from Boundary_Conditions.F90 on Fri Nov  2 16:08:23 2012
MD5:bd9df04a57c49d89e5e67630a680aa6c -- If you edit this, you'll get what you deserve.

(() () () () () () () () () () () () () () () () () () () () () () ()
() () () ())

()

(('add_boundary_condition_surface_elements' 'boundary_conditions' 2 3) (
'add_boundary_condition' 'boundary_conditions' 4 5) (
'apply_dirichlet_conditions' 'boundary_conditions' 6 7 8 9 10 11 12) (
'extract_scalar_surface_field' 'boundary_conditions' 13 14) (
'extract_surface_field' 'boundary_conditions' 15 16 17 18) (
'get_boundary_condition' 'boundary_conditions' 19 20 21 22) (
'get_entire_boundary_condition' 'boundary_conditions' 23 24) (
'get_boundary_condition_nodes' 'boundary_conditions' 25 26) (
'has_boundary_condition' 'boundary_conditions' 27 28) (
'has_boundary_condition_name' 'boundary_conditions' 29 30) (
'has_surface_field' 'boundary_conditions' 31 32) (
'get_boundary_condition_count' 'boundary_conditions' 33 34) (
'insert_surface_field' 'boundary_conditions' 35 36 37 38 39 40) (
'remove_boundary_condition' 'boundary_conditions' 41 42) (
'set_dirichlet_consistent' 'boundary_conditions' 43 44 45) (
'set_reference_node' 'boundary_conditions' 46 47))

(('mpi_fortran_argv_null' 48 0 0 'mpi_fortran_argv_null') (
'mpi_fortran_errcodes_ignore' 49 0 0 'mpi_fortran_errcodes_ignore') (
'mpi_fortran_in_place' 50 0 0 'mpi_fortran_in_place') (
'mpi_fortran_status_ignore' 51 0 0 'mpi_fortran_status_ignore') (
'mpi_fortran_statuses_ignore' 52 0 0 'mpi_fortran_statuses_ignore') (
'petscfortran1' 53 0 0 'petscfortran1') ('petscfortran10' 54 0 0
'petscfortran10') ('petscfortran2' 55 0 0 'petscfortran2') (
'petscfortran3' 56 0 0 'petscfortran3') ('mpi_fortran_argvs_null' 57 0 0
'mpi_fortran_argvs_null') ('mpi_fortran_bottom' 58 0 0
'mpi_fortran_bottom') ('petscfortran4' 59 0 0 'petscfortran4') (
'petscfortran5' 60 0 0 'petscfortran5') ('petscfortran6' 61 0 0
'petscfortran6') ('petscfortran7' 62 0 0 'petscfortran7') (
'petscfortran8' 63 0 0 'petscfortran8') ('petscfortran9' 64 0 0
'petscfortran9'))

()

()

(65 'derive_collapsed_bcs' 'boundary_conditions' 'derive_collapsed_bcs'
1 ((PROCEDURE UNKNOWN-INTENT MODULE-PROC DECL UNKNOWN 0 0 SUBROUTINE
ALWAYS_EXPLICIT) (UNKNOWN 0 0 0 UNKNOWN ()) 66 0 (67 68 69) () 0 () () ()
0 0)
70 'get_dg_surface_mesh' 'boundary_conditions' 'get_dg_surface_mesh' 1 (
(PROCEDURE UNKNOWN-INTENT MODULE-PROC DECL UNKNOWN 0 0 POINTER FUNCTION
ALWAYS_EXPLICIT) (DERIVED 71 0 0 DERIVED ()) 72 0 (73) () 70 () () () 0
0)
74 'get_periodic_boundary_condition' 'boundary_conditions'
'get_periodic_boundary_condition' 1 ((PROCEDURE UNKNOWN-INTENT
MODULE-PROC DECL UNKNOWN 0 0 SUBROUTINE ALWAYS_EXPLICIT) (UNKNOWN 0 0 0
UNKNOWN ()) 75 0 (76 77) () 0 () () () 0 0)
78 'has_scalar_surface_field' 'boundary_conditions'
'has_scalar_surface_field' 1 ((PROCEDURE UNKNOWN-INTENT MODULE-PROC DECL
UNKNOWN 0 0 FUNCTION) (LOGICAL 4 0 0 LOGICAL ()) 79 0 (80 81 82) () 78 ()
() () 0 0)
45 'set_dirichlet_consistent' 'boundary_conditions'
'set_dirichlet_consistent' 1 ((PROCEDURE UNKNOWN-INTENT MODULE-PROC DECL
UNKNOWN 0 0 SUBROUTINE GENERIC ALWAYS_EXPLICIT) (UNKNOWN 0 0 0 UNKNOWN ())
83 0 (84) () 0 () () () 0 0)
71 'mesh_type' 'fields_data_types' 'mesh_type' 1 ((DERIVED
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 POINTER_COMP) (UNKNOWN 0
0 0 UNKNOWN ()) 0 0 () () 0 ((85 'ndglno' (INTEGER 4 0 0 INTEGER ()) (1
0 DEFERRED () ()) (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN
UNKNOWN 0 0 DIMENSION POINTER) UNKNOWN-ACCESS ()) (86 'wrapped' (
LOGICAL 4 0 0 LOGICAL ()) () (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0) UNKNOWN-ACCESS (CONSTANT (LOGICAL 4 0 0 LOGICAL ())
0 1)) (87 'shape' (DERIVED 88 0 0 DERIVED ()) () (UNKNOWN-FL
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0) UNKNOWN-ACCESS (
STRUCTURE (DERIVED 88 0 0 DERIVED ()) 0 ((() ()) (() ()) (() ()) (() ())
((NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0) ()) ((NULL (UNKNOWN 0 0 0 UNKNOWN ())
0) ()) ((NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0) ()) ((NULL (UNKNOWN 0 0 0
UNKNOWN ()) 0) ()) ((NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0) ()) ((NULL (
UNKNOWN 0 0 0 UNKNOWN ()) 0) ()) ((NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0) ())
((STRUCTURE (DERIVED 89 0 0 DERIVED ()) 0 ((() ()) (() ()) (() ()) (() ())
((NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0) ()) ((NULL (UNKNOWN 0 0 0 UNKNOWN ())
0) ()) (() ()) ((NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0) ()) (() ())) ()) ())
((NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0) ()) ((NULL (UNKNOWN 0 0 0 UNKNOWN ())
0) ()) ((NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0) ()) ((NULL (UNKNOWN 0 0 0
UNKNOWN ()) 0) ()) (() ())) ())) (90 'elements' (INTEGER 4 0 0 INTEGER ())
() (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0)
UNKNOWN-ACCESS ()) (91 'nodes' (INTEGER 4 0 0 INTEGER ()) () (
UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0)
UNKNOWN-ACCESS ()) (92 'name' (CHARACTER 1 0 0 CHARACTER ((CONSTANT (
INTEGER 4 0 0 INTEGER ()) 0 '101'))) () (UNKNOWN-FL UNKNOWN-INTENT
UNKNOWN-PROC UNKNOWN UNKNOWN 0 0) UNKNOWN-ACCESS ()) (93 'option_path' (
CHARACTER 1 0 0 CHARACTER ((CONSTANT (INTEGER 4 0 0 INTEGER ()) 0 '8192')))
() (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0)
UNKNOWN-ACCESS (CONSTANT (CHARACTER 1 0 0 CHARACTER (())) 0 8192
'/uninitialised_path/                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            '))
(94 'continuity' (INTEGER 4 0 0 INTEGER ()) () (UNKNOWN-FL
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0) UNKNOWN-ACCESS (
CONSTANT (INTEGER 4 0 0 INTEGER ()) 0 '0')) (95 'refcount' (DERIVED 96 0
0 DERIVED ()) () (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN
0 0 POINTER) UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0)) (97
'faces' (DERIVED 98 0 0 DERIVED ()) () (UNKNOWN-FL UNKNOWN-INTENT
UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 POINTER) UNKNOWN-ACCESS (NULL (UNKNOWN
0 0 0 UNKNOWN ()) 0)) (99 'subdomain_mesh' (DERIVED 100 0 0 DERIVED ())
() (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 POINTER)
UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0)) (101 'adj_lists' (
DERIVED 102 0 0 DERIVED ()) () (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0 POINTER) UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN
()) 0)) (103 'columns' (INTEGER 4 0 0 INTEGER ()) (1 0 DEFERRED () ()) (
UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 DIMENSION
POINTER) UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0)) (104
'element_columns' (INTEGER 4 0 0 INTEGER ()) (1 0 DEFERRED () ()) (
UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 DIMENSION
POINTER) UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0)) (105
'region_ids' (INTEGER 4 0 0 INTEGER ()) (1 0 DEFERRED () ()) (
UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 DIMENSION
POINTER) UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0)) (106 'halos'
(DERIVED 107 0 0 DERIVED ()) (1 0 DEFERRED () ()) (UNKNOWN-FL
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 DIMENSION POINTER)
UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0)) (108 'element_halos'
(DERIVED 107 0 0 DERIVED ()) (1 0 DEFERRED () ()) (UNKNOWN-FL
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 DIMENSION POINTER)
UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0)) (109 'periodic' (
LOGICAL 4 0 0 LOGICAL ()) () (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0) UNKNOWN-ACCESS (CONSTANT (LOGICAL 4 0 0 LOGICAL ())
0 0))) PUBLIC (() () () ()) () 0 0 2572791)
48 'mpi_argv_null' 'mpi_interfaces' 'mpi_argv_null' 1 ((VARIABLE
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 DIMENSION IN_COMMON) (
CHARACTER 1 0 0 CHARACTER ((CONSTANT (INTEGER 4 0 0 INTEGER ()) 0 '1')))
0 0 () (1 0 EXPLICIT (CONSTANT (INTEGER 4 0 0 INTEGER ()) 0 '1') (
CONSTANT (INTEGER 4 0 0 INTEGER ()) 0 '1')) 0 () () () 0 0)
56 'petsc_null' 'petscsysdef' 'petsc_null' 1 ((VARIABLE UNKNOWN-INTENT
UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 IN_COMMON) (INTEGER 4 0 0 INTEGER ()) 0
0 () () 0 () () () 0 0)
59 'petsc_null_scalar' 'petscsysdef' 'petsc_null_scalar' 1 ((VARIABLE
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 IN_COMMON) (REAL 8 0 0
REAL ()) 0 0 () () 0 () () () 0 0)
64 'petsc_comm_world' 'petscsysdef' 'petsc_comm_world' 1 ((VARIABLE
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 IN_COMMON) (INTEGER 4 0
0 INTEGER ()) 0 0 () () 0 () () () 0 0)
63 'petsc_null_object' 'petscsysdef' 'petsc_null_object' 1 ((VARIABLE
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 IN_COMMON) (INTEGER 8 0
0 INTEGER ()) 0 0 () () 0 () () () 0 0)
62 'petsc_null_truth' 'petscsysdef' 'petsc_null_truth' 1 ((VARIABLE
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 IN_COMMON) (LOGICAL 4 0
0 LOGICAL ()) 0 0 () () 0 () () () 0 0)
61 'petsc_null_real' 'petscsysdef' 'petsc_null_real' 1 ((VARIABLE
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 IN_COMMON) (REAL 8 0 0
REAL ()) 0 0 () () 0 () () () 0 0)
58 'mpi_bottom' 'mpi_interfaces' 'mpi_bottom' 1 ((VARIABLE
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 IN_COMMON) (INTEGER 4 0
0 INTEGER ()) 0 0 () () 0 () () () 0 0)
57 'mpi_argvs_null' 'mpi_interfaces' 'mpi_argvs_null' 1 ((VARIABLE
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 IN_COMMON) (REAL 8 0 0
REAL ()) 0 0 () () 0 () () () 0 0)
60 'petsc_null_double' 'petscsysdef' 'petsc_null_double' 1 ((VARIABLE
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 IN_COMMON) (REAL 8 0 0
REAL ()) 0 0 () () 0 () () () 0 0)
52 'mpi_statuses_ignore' 'mpi_interfaces' 'mpi_statuses_ignore' 1 ((
VARIABLE UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 IN_COMMON) (
REAL 8 0 0 REAL ()) 0 0 () () 0 () () () 0 0)
51 'mpi_status_ignore' 'mpi_interfaces' 'mpi_status_ignore' 1 ((
VARIABLE UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 DIMENSION
IN_COMMON) (INTEGER 4 0 0 INTEGER ()) 0 0 () (1 0 EXPLICIT (CONSTANT (
INTEGER 4 0 0 INTEGER ()) 0 '1') (CONSTANT (INTEGER 4 0 0 INTEGER ()) 0
'5')) 0 () () () 0 0)
53 'petsc_null_character' 'petscsysdef' 'petsc_null_character' 1 ((
VARIABLE UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN IMPLICIT-SAVE 0 0 IN_COMMON)
(CHARACTER 1 0 0 CHARACTER ((CONSTANT (INTEGER 4 0 0 INTEGER ()) 0 '80')))
0 0 () () 0 () () () 0 0)
50 'mpi_in_place' 'mpi_interfaces' 'mpi_in_place' 1 ((VARIABLE
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 IN_COMMON) (INTEGER 4 0
0 INTEGER ()) 0 0 () () 0 () () () 0 0)
49 'mpi_errcodes_ignore' 'mpi_interfaces' 'mpi_errcodes_ignore' 1 ((
VARIABLE UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 DIMENSION
IN_COMMON) (INTEGER 4 0 0 INTEGER ()) 0 0 () (1 0 EXPLICIT (CONSTANT (
INTEGER 4 0 0 INTEGER ()) 0 '1') (CONSTANT (INTEGER 4 0 0 INTEGER ()) 0
'1')) 0 () () () 0 0)
55 'petsc_null_integer' 'petscsysdef' 'petsc_null_integer' 1 ((VARIABLE
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 IN_COMMON) (INTEGER 4 0
0 INTEGER ()) 0 0 () () 0 () () () 0 0)
88 'element_type' 'elements' 'element_type' 1 ((DERIVED UNKNOWN-INTENT
UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 POINTER_COMP) (UNKNOWN 0 0 0 UNKNOWN ())
0 0 () () 0 ((110 'dim' (INTEGER 4 0 0 INTEGER ()) () (UNKNOWN-FL
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0) UNKNOWN-ACCESS ()) (
111 'loc' (INTEGER 4 0 0 INTEGER ()) () (UNKNOWN-FL UNKNOWN-INTENT
UNKNOWN-PROC UNKNOWN UNKNOWN 0 0) UNKNOWN-ACCESS ()) (112 'ngi' (
INTEGER 4 0 0 INTEGER ()) () (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0) UNKNOWN-ACCESS ()) (113 'degree' (INTEGER 4 0 0
INTEGER ()) () (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0
0) UNKNOWN-ACCESS ()) (114 'n' (REAL 8 0 0 REAL ()) (2 0 DEFERRED () ()
() ()) (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DIMENSION POINTER) UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0)) (
115 'dn' (REAL 8 0 0 REAL ()) (3 0 DEFERRED () () () () () ()) (
UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 DIMENSION
POINTER) UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0)) (116 'n_s'
(REAL 8 0 0 REAL ()) (3 0 DEFERRED () () () () () ()) (UNKNOWN-FL
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 DIMENSION POINTER)
UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0)) (117 'dn_s' (REAL 8
0 0 REAL ()) (4 0 DEFERRED () () () () () () () ()) (UNKNOWN-FL
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 DIMENSION POINTER)
UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0)) (118 'spoly' (
DERIVED 119 0 0 DERIVED ()) (2 0 DEFERRED () () () ()) (UNKNOWN-FL
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 DIMENSION POINTER)
UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0)) (120 'dspoly' (
DERIVED 119 0 0 DERIVED ()) (2 0 DEFERRED () () () ()) (UNKNOWN-FL
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 DIMENSION POINTER)
UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0)) (121 'numbering' (
DERIVED 122 0 0 DERIVED ()) () (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0 POINTER) UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN
()) 0)) (123 'quadrature' (DERIVED 89 0 0 DERIVED ()) () (UNKNOWN-FL
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0) UNKNOWN-ACCESS (
STRUCTURE (DERIVED 89 0 0 DERIVED ()) 0 ((() ()) (() ()) (() ()) (() ())
((NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0) ()) ((NULL (UNKNOWN 0 0 0 UNKNOWN ())
0) ()) (() ()) ((NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0) ()) (() ())) ())) (
124 'surface_quadrature' (DERIVED 89 0 0 DERIVED ()) () (UNKNOWN-FL
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 POINTER) UNKNOWN-ACCESS
(NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0)) (125 'superconvergence' (DERIVED
126 0 0 DERIVED ()) () (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN
UNKNOWN 0 0 POINTER) UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0))
(127 'constraints' (DERIVED 128 0 0 DERIVED ()) () (UNKNOWN-FL
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 POINTER) UNKNOWN-ACCESS
(NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0)) (129 'refcount' (DERIVED 96 0 0
DERIVED ()) () (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0
0 POINTER) UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0)) (130 'name'
(CHARACTER 1 0 0 CHARACTER ((CONSTANT (INTEGER 4 0 0 INTEGER ()) 0 '0')))
() (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0)
UNKNOWN-ACCESS ())) PUBLIC (() () () ()) () 0 0 79461029)
107 'halo_type' 'halo_data_types' 'halo_type' 1 ((DERIVED UNKNOWN-INTENT
UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 POINTER_COMP) (UNKNOWN 0 0 0 UNKNOWN ())
0 0 () () 0 ((131 'name' (CHARACTER 1 0 0 CHARACTER ((CONSTANT (INTEGER
4 0 0 INTEGER ()) 0 '101'))) () (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0) UNKNOWN-ACCESS ()) (132 'refcount' (DERIVED 96 0 0
DERIVED ()) () (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0
0 POINTER) UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0)) (133
'data_type' (INTEGER 4 0 0 INTEGER ()) () (UNKNOWN-FL UNKNOWN-INTENT
UNKNOWN-PROC UNKNOWN UNKNOWN 0 0) UNKNOWN-ACCESS (CONSTANT (INTEGER 4 0
0 INTEGER ()) 0 '0')) (134 'ordering_scheme' (INTEGER 4 0 0 INTEGER ())
() (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0)
UNKNOWN-ACCESS (CONSTANT (INTEGER 4 0 0 INTEGER ()) 0 '0')) (135
'communicator' (INTEGER 4 0 0 INTEGER ()) () (UNKNOWN-FL UNKNOWN-INTENT
UNKNOWN-PROC UNKNOWN UNKNOWN 0 0) UNKNOWN-ACCESS ()) (136 'nprocs' (
INTEGER 4 0 0 INTEGER ()) () (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0) UNKNOWN-ACCESS (CONSTANT (INTEGER 4 0 0 INTEGER ())
0 '0')) (137 'sends' (DERIVED 138 0 0 DERIVED ()) (1 0 DEFERRED () ()) (
UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 DIMENSION
POINTER) UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0)) (139
'receives' (DERIVED 138 0 0 DERIVED ()) (1 0 DEFERRED () ()) (
UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 DIMENSION
POINTER) UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0)) (140
'nowned_nodes' (INTEGER 4 0 0 INTEGER ()) () (UNKNOWN-FL UNKNOWN-INTENT
UNKNOWN-PROC UNKNOWN UNKNOWN 0 0) UNKNOWN-ACCESS (CONSTANT (INTEGER 4 0
0 INTEGER ()) 0 '-1')) (141 'owners' (INTEGER 4 0 0 INTEGER ()) (1 0
DEFERRED () ()) (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN
0 0 DIMENSION POINTER) UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0))
(142 'unn_count' (INTEGER 4 0 0 INTEGER ()) () (UNKNOWN-FL
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0) UNKNOWN-ACCESS (
CONSTANT (INTEGER 4 0 0 INTEGER ()) 0 '-1')) (143 'owned_nodes_unn_base'
(INTEGER 4 0 0 INTEGER ()) (1 0 DEFERRED () ()) (UNKNOWN-FL
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 DIMENSION POINTER)
UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0)) (144
'my_owned_nodes_unn_base' (INTEGER 4 0 0 INTEGER ()) () (UNKNOWN-FL
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0) UNKNOWN-ACCESS (
CONSTANT (INTEGER 4 0 0 INTEGER ()) 0 '-1')) (145 'receives_gnn_to_unn'
(INTEGER 4 0 0 INTEGER ()) (1 0 DEFERRED () ()) (UNKNOWN-FL
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 DIMENSION POINTER)
UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0)) (146 'gnn_to_unn' (
INTEGER 4 0 0 INTEGER ()) (1 0 DEFERRED () ()) (UNKNOWN-FL
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 DIMENSION POINTER)
UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0))) PUBLIC (() () () ())
() 0 0 63994053)
20 'get_scalar_boundary_condition_by_name' 'boundary_conditions'
'get_scalar_boundary_condition_by_name' 1 ((PROCEDURE UNKNOWN-INTENT
MODULE-PROC DECL UNKNOWN 0 0 SUBROUTINE ALWAYS_EXPLICIT) (UNKNOWN 0 0 0
UNKNOWN ()) 147 0 (148 149 150 151 152 153 154) () 0 () () () 0 0)
19 'get_vector_boundary_condition_by_name' 'boundary_conditions'
'get_vector_boundary_condition_by_name' 1 ((PROCEDURE UNKNOWN-INTENT
MODULE-PROC DECL UNKNOWN 0 0 SUBROUTINE ALWAYS_EXPLICIT) (UNKNOWN 0 0 0
UNKNOWN ()) 155 0 (156 157 158 159 160 161 162 163) () 0 () () () 0 0)
40 'insert_scalar_surface_field' 'boundary_conditions'
'insert_scalar_surface_field' 1 ((PROCEDURE UNKNOWN-INTENT MODULE-PROC
DECL UNKNOWN 0 0 SUBROUTINE) (UNKNOWN 0 0 0 UNKNOWN ()) 164 0 (165 166
167) () 0 () () () 0 0)
39 'insert_vector_surface_field' 'boundary_conditions'
'insert_vector_surface_field' 1 ((PROCEDURE UNKNOWN-INTENT MODULE-PROC
DECL UNKNOWN 0 0 SUBROUTINE) (UNKNOWN 0 0 0 UNKNOWN ()) 168 0 (169 170
171) () 0 () () () 0 0)
38 'insert_scalar_surface_field_by_name' 'boundary_conditions'
'insert_scalar_surface_field_by_name' 1 ((PROCEDURE UNKNOWN-INTENT
MODULE-PROC DECL UNKNOWN 0 0 SUBROUTINE) (UNKNOWN 0 0 0 UNKNOWN ()) 172
0 (173 174 175) () 0 () () () 0 0)
37 'insert_vector_surface_field_by_name' 'boundary_conditions'
'insert_vector_surface_field_by_name' 1 ((PROCEDURE UNKNOWN-INTENT
MODULE-PROC DECL UNKNOWN 0 0 SUBROUTINE) (UNKNOWN 0 0 0 UNKNOWN ()) 176
0 (177 178 179) () 0 () () () 0 0)
36 'insert_vector_scalar_surface_field' 'boundary_conditions'
'insert_vector_scalar_surface_field' 1 ((PROCEDURE UNKNOWN-INTENT
MODULE-PROC DECL UNKNOWN 0 0 SUBROUTINE) (UNKNOWN 0 0 0 UNKNOWN ()) 180
0 (181 182 183) () 0 () () () 0 0)
35 'insert_vector_scalar_surface_field_by_name' 'boundary_conditions'
'insert_vector_scalar_surface_field_by_name' 1 ((PROCEDURE
UNKNOWN-INTENT MODULE-PROC DECL UNKNOWN 0 0 SUBROUTINE) (UNKNOWN 0 0 0
UNKNOWN ()) 184 0 (185 186 187) () 0 () () () 0 0)
18 'extract_scalar_surface_field_by_number' 'boundary_conditions'
'extract_scalar_surface_field_by_number' 1 ((PROCEDURE UNKNOWN-INTENT
MODULE-PROC DECL UNKNOWN 0 0 POINTER FUNCTION ALWAYS_EXPLICIT) (DERIVED
188 0 0 DERIVED ()) 189 0 (190 191 192) () 193 () () () 0 0)
17 'extract_vector_surface_field_by_number' 'boundary_conditions'
'extract_vector_surface_field_by_number' 1 ((PROCEDURE UNKNOWN-INTENT
MODULE-PROC DECL UNKNOWN 0 0 POINTER FUNCTION ALWAYS_EXPLICIT) (DERIVED
194 0 0 DERIVED ()) 195 0 (196 197 198) () 199 () () () 0 0)
16 'extract_scalar_surface_field_by_name' 'boundary_conditions'
'extract_scalar_surface_field_by_name' 1 ((PROCEDURE UNKNOWN-INTENT
MODULE-PROC DECL UNKNOWN 0 0 POINTER FUNCTION ALWAYS_EXPLICIT) (DERIVED
188 0 0 DERIVED ()) 200 0 (201 202 203) () 204 () () () 0 0)
15 'extract_vector_surface_field_by_name' 'boundary_conditions'
'extract_vector_surface_field_by_name' 1 ((PROCEDURE UNKNOWN-INTENT
MODULE-PROC DECL UNKNOWN 0 0 POINTER FUNCTION ALWAYS_EXPLICIT) (DERIVED
194 0 0 DERIVED ()) 205 0 (206 207 208) () 209 () () () 0 0)
14 'extract_vector_scalar_surface_field' 'boundary_conditions'
'extract_vector_scalar_surface_field' 1 ((PROCEDURE UNKNOWN-INTENT
MODULE-PROC DECL UNKNOWN 0 0 POINTER FUNCTION ALWAYS_EXPLICIT) (DERIVED
188 0 0 DERIVED ()) 210 0 (211 212 213) () 214 () () () 0 0)
13 'extract_vector_scalar_surface_field_by_name' 'boundary_conditions'
'extract_vector_scalar_surface_field_by_name' 1 ((PROCEDURE
UNKNOWN-INTENT MODULE-PROC DECL UNKNOWN 0 0 POINTER FUNCTION
ALWAYS_EXPLICIT) (DERIVED 188 0 0 DERIVED ()) 215 0 (216 217 218) () 219
() () () 0 0)
32 'has_scalar_surface_field_specific' 'boundary_conditions'
'has_scalar_surface_field_specific' 1 ((PROCEDURE UNKNOWN-INTENT
MODULE-PROC DECL UNKNOWN 0 0 FUNCTION) (LOGICAL 4 0 0 LOGICAL ()) 220 0
(221 222 223) () 32 () () () 0 0)
31 'has_vector_surface_field_specific' 'boundary_conditions'
'has_vector_surface_field_specific' 1 ((PROCEDURE UNKNOWN-INTENT
MODULE-PROC DECL UNKNOWN 0 0 FUNCTION) (LOGICAL 4 0 0 LOGICAL ()) 224 0
(225 226 227) () 31 () () () 0 0)
34 'get_scalar_boundary_condition_count' 'boundary_conditions'
'get_scalar_boundary_condition_count' 1 ((PROCEDURE UNKNOWN-INTENT
MODULE-PROC DECL UNKNOWN 0 0 FUNCTION IMPLICIT_PURE) (INTEGER 4 0 0
INTEGER ()) 228 0 (229) () 34 () () () 0 0)
33 'get_vector_boundary_condition_count' 'boundary_conditions'
'get_vector_boundary_condition_count' 1 ((PROCEDURE UNKNOWN-INTENT
MODULE-PROC DECL UNKNOWN 0 0 FUNCTION IMPLICIT_PURE) (INTEGER 4 0 0
INTEGER ()) 230 0 (231) () 33 () () () 0 0)
24 'get_entire_scalar_boundary_condition' 'boundary_conditions'
'get_entire_scalar_boundary_condition' 1 ((PROCEDURE UNKNOWN-INTENT
MODULE-PROC DECL UNKNOWN 0 0 SUBROUTINE ALWAYS_EXPLICIT) (UNKNOWN 0 0 0
UNKNOWN ()) 232 0 (233 234 235 236 237 238) () 0 () () () 0 0)
23 'get_entire_vector_boundary_condition' 'boundary_conditions'
'get_entire_vector_boundary_condition' 1 ((PROCEDURE UNKNOWN-INTENT
MODULE-PROC DECL UNKNOWN 0 0 SUBROUTINE ALWAYS_EXPLICIT) (UNKNOWN 0 0 0
UNKNOWN ()) 239 0 (240 241 242 243 244) () 0 () () () 0 0)
26 'get_scalar_boundary_condition_nodes' 'boundary_conditions'
'get_scalar_boundary_condition_nodes' 1 ((PROCEDURE UNKNOWN-INTENT
MODULE-PROC DECL UNKNOWN 0 0 SUBROUTINE ALWAYS_EXPLICIT) (UNKNOWN 0 0 0
UNKNOWN ()) 245 0 (246 247 248) () 0 () () () 0 0)
25 'get_vector_boundary_condition_nodes' 'boundary_conditions'
'get_vector_boundary_condition_nodes' 1 ((PROCEDURE UNKNOWN-INTENT
MODULE-PROC DECL UNKNOWN 0 0 SUBROUTINE ALWAYS_EXPLICIT) (UNKNOWN 0 0 0
UNKNOWN ()) 249 0 (250 251 252) () 0 () () () 0 0)
47 'set_reference_node_scalar' 'boundary_conditions'
'set_reference_node_scalar' 1 ((PROCEDURE UNKNOWN-INTENT MODULE-PROC
DECL UNKNOWN 0 0 SUBROUTINE ALWAYS_EXPLICIT) (UNKNOWN 0 0 0 UNKNOWN ())
253 0 (254 255 256 257 258) () 0 () () () 0 0)
46 'set_reference_node_vector_petsc' 'boundary_conditions'
'set_reference_node_vector_petsc' 1 ((PROCEDURE UNKNOWN-INTENT
MODULE-PROC DECL UNKNOWN 0 0 SUBROUTINE ALWAYS_EXPLICIT) (UNKNOWN 0 0 0
UNKNOWN ()) 259 0 (260 261 262 263 264) () 0 () () () 0 0)
28 'has_boundary_condition_scalar' 'boundary_conditions'
'has_boundary_condition_scalar' 1 ((PROCEDURE UNKNOWN-INTENT MODULE-PROC
DECL UNKNOWN 0 0 FUNCTION IMPLICIT_PURE) (LOGICAL 4 0 0 LOGICAL ()) 265
0 (266 267) () 28 () () () 0 0)
27 'has_boundary_condition_vector' 'boundary_conditions'
'has_boundary_condition_vector' 1 ((PROCEDURE UNKNOWN-INTENT MODULE-PROC
DECL UNKNOWN 0 0 FUNCTION IMPLICIT_PURE) (LOGICAL 4 0 0 LOGICAL ()) 268
0 (269 270) () 27 () () () 0 0)
30 'has_boundary_condition_name_scalar' 'boundary_conditions'
'has_boundary_condition_name_scalar' 1 ((PROCEDURE UNKNOWN-INTENT
MODULE-PROC DECL UNKNOWN 0 0 FUNCTION IMPLICIT_PURE) (LOGICAL 4 0 0
LOGICAL ()) 271 0 (272 273) () 30 () () () 0 0)
29 'has_boundary_condition_name_vector' 'boundary_conditions'
'has_boundary_condition_name_vector' 1 ((PROCEDURE UNKNOWN-INTENT
MODULE-PROC DECL UNKNOWN 0 0 FUNCTION IMPLICIT_PURE) (LOGICAL 4 0 0
LOGICAL ()) 274 0 (275 276) () 29 () () () 0 0)
44 'set_dirichlet_consistent_scalar' 'boundary_conditions'
'set_dirichlet_consistent_scalar' 1 ((PROCEDURE UNKNOWN-INTENT
MODULE-PROC DECL UNKNOWN 0 0 SUBROUTINE) (UNKNOWN 0 0 0 UNKNOWN ()) 277
0 (278) () 0 () () () 0 0)
43 'set_dirichlet_consistent_vector' 'boundary_conditions'
'set_dirichlet_consistent_vector' 1 ((PROCEDURE UNKNOWN-INTENT
MODULE-PROC DECL UNKNOWN 0 0 SUBROUTINE) (UNKNOWN 0 0 0 UNKNOWN ()) 279
0 (280) () 0 () () () 0 0)
12 'apply_dirichlet_conditions_scalar' 'boundary_conditions'
'apply_dirichlet_conditions_scalar' 1 ((PROCEDURE UNKNOWN-INTENT
MODULE-PROC DECL UNKNOWN 0 0 SUBROUTINE ALWAYS_EXPLICIT) (UNKNOWN 0 0 0
UNKNOWN ()) 281 0 (282 283 284 285) () 0 () () () 0 0)
11 'apply_dirichlet_conditions_scalar_lumped' 'boundary_conditions'
'apply_dirichlet_conditions_scalar_lumped' 1 ((PROCEDURE UNKNOWN-INTENT
MODULE-PROC DECL UNKNOWN 0 0 SUBROUTINE) (UNKNOWN 0 0 0 UNKNOWN ()) 286
0 (287 288 289 290) () 0 () () () 0 0)
10 'apply_dirichlet_conditions_vector' 'boundary_conditions'
'apply_dirichlet_conditions_vector' 1 ((PROCEDURE UNKNOWN-INTENT
MODULE-PROC DECL UNKNOWN 0 0 SUBROUTINE ALWAYS_EXPLICIT) (UNKNOWN 0 0 0
UNKNOWN ()) 291 0 (292 293 294 295) () 0 () () () 0 0)
9 'apply_dirichlet_conditions_vector_petsc_csr' 'boundary_conditions'
'apply_dirichlet_conditions_vector_petsc_csr' 1 ((PROCEDURE
UNKNOWN-INTENT MODULE-PROC DECL UNKNOWN 0 0 SUBROUTINE ALWAYS_EXPLICIT)
(UNKNOWN 0 0 0 UNKNOWN ()) 296 0 (297 298 299 300) () 0 () () () 0 0)
8 'apply_dirichlet_conditions_vector_component' 'boundary_conditions'
'apply_dirichlet_conditions_vector_component' 1 ((PROCEDURE
UNKNOWN-INTENT MODULE-PROC DECL UNKNOWN 0 0 SUBROUTINE ALWAYS_EXPLICIT)
(UNKNOWN 0 0 0 UNKNOWN ()) 301 0 (302 303 304 305 306) () 0 () () () 0 0)
7 'apply_dirichlet_conditions_vector_lumped' 'boundary_conditions'
'apply_dirichlet_conditions_vector_lumped' 1 ((PROCEDURE UNKNOWN-INTENT
MODULE-PROC DECL UNKNOWN 0 0 SUBROUTINE ALWAYS_EXPLICIT) (UNKNOWN 0 0 0
UNKNOWN ()) 307 0 (308 309 310 311) () 0 () () () 0 0)
6 'apply_dirichlet_conditions_vector_component_lumped'
'boundary_conditions' 'apply_dirichlet_conditions_vector_component_lumped'
1 ((PROCEDURE UNKNOWN-INTENT MODULE-PROC DECL UNKNOWN 0 0 SUBROUTINE
ALWAYS_EXPLICIT) (UNKNOWN 0 0 0 UNKNOWN ()) 312 0 (313 314 315 316 317)
() 0 () () () 0 0)
165 'field' '' 'field' 164 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0
0 DUMMY) (DERIVED 188 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
166 'n' '' 'n' 164 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 DUMMY)
(INTEGER 4 0 0 INTEGER ()) 0 0 () () 0 () () () 0 0)
167 'surface_field' '' 'surface_field' 164 ((VARIABLE IN UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0 DUMMY) (DERIVED 188 0 0 DERIVED ()) 0 0 () () 0 () ()
() 0 0)
169 'field' '' 'field' 168 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0
0 DUMMY) (DERIVED 194 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
170 'n' '' 'n' 168 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 DUMMY)
(INTEGER 4 0 0 INTEGER ()) 0 0 () () 0 () () () 0 0)
171 'surface_field' '' 'surface_field' 168 ((VARIABLE IN UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0 DUMMY) (DERIVED 194 0 0 DERIVED ()) 0 0 () () 0 () ()
() 0 0)
181 'field' '' 'field' 180 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0
0 DUMMY) (DERIVED 194 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
182 'n' '' 'n' 180 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 DUMMY)
(INTEGER 4 0 0 INTEGER ()) 0 0 () () 0 () () () 0 0)
102 'adjacency_cache' 'fields_data_types' 'adjacency_cache' 1 ((DERIVED
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 POINTER_COMP) (UNKNOWN 0
0 0 UNKNOWN ()) 0 0 () () 0 ((318 'nnlist' (DERIVED 319 0 0 DERIVED ())
() (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 POINTER)
UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0)) (320 'nelist' (
DERIVED 319 0 0 DERIVED ()) () (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0 POINTER) UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN
()) 0)) (321 'eelist' (DERIVED 319 0 0 DERIVED ()) () (UNKNOWN-FL
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 POINTER) UNKNOWN-ACCESS
(NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0))) PUBLIC (() () () ()) () 0 0
37419158)
98 'mesh_faces' 'fields_data_types' 'mesh_faces' 1 ((DERIVED
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 POINTER_COMP) (UNKNOWN 0
0 0 UNKNOWN ()) 0 0 () () 0 ((322 'shape' (DERIVED 88 0 0 DERIVED ()) ()
(UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 POINTER)
UNKNOWN-ACCESS ()) (323 'face_list' (DERIVED 324 0 0 DERIVED ()) () (
UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0)
UNKNOWN-ACCESS (STRUCTURE (DERIVED 324 0 0 DERIVED ()) 0 (((STRUCTURE (
DERIVED 319 0 0 DERIVED ()) 0 (((NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0) ())
((NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0) ()) ((NULL (UNKNOWN 0 0 0 UNKNOWN ())
0) ()) (() ()) ((NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0) ()) ((NULL (UNKNOWN
0 0 0 UNKNOWN ()) 0) ()) ((NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0) ()) ((
CONSTANT (CHARACTER 1 0 0 CHARACTER (())) 0 101
'                                                                                                     ')
()) ((CONSTANT (LOGICAL 4 0 0 LOGICAL ()) 0 0) ()) ((CONSTANT (LOGICAL 4
0 0 LOGICAL ()) 0 0) ())) ()) ()) ((NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0) ())
((NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0) ()) ((CONSTANT (LOGICAL 4 0 0
LOGICAL ()) 0 0) ()) ((CONSTANT (LOGICAL 4 0 0 LOGICAL ()) 0 0) ()) (()
()) ((NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0) ()) ((NULL (UNKNOWN 0 0 0
UNKNOWN ()) 0) ()) ((CONSTANT (CHARACTER 1 0 0 CHARACTER (())) 0 101
'                                                                                                     ')
())) ())) (325 'face_lno' (INTEGER 4 0 0 INTEGER ()) (1 0 DEFERRED () ())
(UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 DIMENSION
POINTER) UNKNOWN-ACCESS ()) (326 'surface_mesh' (DERIVED 71 0 0 DERIVED
()) () (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0)
UNKNOWN-ACCESS (STRUCTURE (DERIVED 71 0 0 DERIVED ()) 0 ((() ()) ((
CONSTANT (LOGICAL 4 0 0 LOGICAL ()) 0 1) ()) ((STRUCTURE (DERIVED 88 0 0
DERIVED ()) 0 ((() ()) (() ()) (() ()) (() ()) ((NULL (UNKNOWN 0 0 0
UNKNOWN ()) 0) ()) ((NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0) ()) ((NULL (
UNKNOWN 0 0 0 UNKNOWN ()) 0) ()) ((NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0) ())
((NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0) ()) ((NULL (UNKNOWN 0 0 0 UNKNOWN ())
0) ()) ((NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0) ()) ((STRUCTURE (DERIVED 89
0 0 DERIVED ()) 0 ((() ()) (() ()) (() ()) (() ()) ((NULL (UNKNOWN 0 0 0
UNKNOWN ()) 0) ()) ((NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0) ()) (() ()) ((
NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0) ()) (() ())) ()) ()) ((NULL (UNKNOWN
0 0 0 UNKNOWN ()) 0) ()) ((NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0) ()) ((
NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0) ()) ((NULL (UNKNOWN 0 0 0 UNKNOWN ())
0) ()) (() ())) ()) ()) (() ()) (() ()) (() ()) ((CONSTANT (CHARACTER 1
0 0 CHARACTER (())) 0 8192
'/uninitialised_path/                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            ')
()) ((CONSTANT (INTEGER 4 0 0 INTEGER ()) 0 '0') ()) ((NULL (UNKNOWN 0 0
0 UNKNOWN ()) 0) ()) ((NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0) ()) ((NULL (
UNKNOWN 0 0 0 UNKNOWN ()) 0) ()) ((NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0) ())
((NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0) ()) ((NULL (UNKNOWN 0 0 0 UNKNOWN ())
0) ()) ((NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0) ()) ((NULL (UNKNOWN 0 0 0
UNKNOWN ()) 0) ()) ((NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0) ()) ((CONSTANT (
LOGICAL 4 0 0 LOGICAL ()) 0 0) ())) ())) (327 'surface_node_list' (
INTEGER 4 0 0 INTEGER ()) (1 0 DEFERRED () ()) (UNKNOWN-FL
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 DIMENSION POINTER)
UNKNOWN-ACCESS ()) (328 'face_element_list' (INTEGER 4 0 0 INTEGER ()) (
1 0 DEFERRED () ()) (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN
UNKNOWN 0 0 DIMENSION POINTER) UNKNOWN-ACCESS ()) (329 'boundary_ids' (
INTEGER 4 0 0 INTEGER ()) (1 0 DEFERRED () ()) (UNKNOWN-FL
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 DIMENSION POINTER)
UNKNOWN-ACCESS ()) (330 'coplanar_ids' (INTEGER 4 0 0 INTEGER ()) (1 0
DEFERRED () ()) (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN
0 0 DIMENSION POINTER) UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0))
(331 'dg_surface_mesh' (DERIVED 71 0 0 DERIVED ()) () (UNKNOWN-FL
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 POINTER) UNKNOWN-ACCESS
(NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0)) (332 'has_internal_boundaries' (
LOGICAL 4 0 0 LOGICAL ()) () (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0) UNKNOWN-ACCESS (CONSTANT (LOGICAL 4 0 0 LOGICAL ())
0 0))) PUBLIC (() () () ()) () 0 0 11936185)
100 'mesh_subdomain_mesh' 'fields_data_types' 'mesh_subdomain_mesh' 1 (
(DERIVED UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 POINTER_COMP) (
UNKNOWN 0 0 0 UNKNOWN ()) 0 0 () () 0 ((333 'element_list' (INTEGER 4 0
0 INTEGER ()) (1 0 DEFERRED () ()) (UNKNOWN-FL UNKNOWN-INTENT
UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 DIMENSION POINTER) UNKNOWN-ACCESS ()) (
334 'node_list' (INTEGER 4 0 0 INTEGER ()) (1 0 DEFERRED () ()) (
UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 DIMENSION
POINTER) UNKNOWN-ACCESS ())) PUBLIC (() () () ()) () 0 0 85302181)
188 'scalar_field' 'fields_data_types' 'scalar_field' 1 ((DERIVED
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 POINTER_COMP) (UNKNOWN 0
0 0 UNKNOWN ()) 0 0 () () 0 ((335 'val' (REAL 8 0 0 REAL ()) (1 0
DEFERRED () ()) (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN
0 0 DIMENSION POINTER) UNKNOWN-ACCESS ()) (336 'val_stride' (INTEGER 4 0
0 INTEGER ()) () (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN
0 0) UNKNOWN-ACCESS (CONSTANT (INTEGER 4 0 0 INTEGER ()) 0 '1')) (337
'wrapped' (LOGICAL 4 0 0 LOGICAL ()) () (UNKNOWN-FL UNKNOWN-INTENT
UNKNOWN-PROC UNKNOWN UNKNOWN 0 0) UNKNOWN-ACCESS (CONSTANT (LOGICAL 4 0
0 LOGICAL ()) 0 1)) (338 'field_type' (INTEGER 4 0 0 INTEGER ()) () (
UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0)
UNKNOWN-ACCESS (CONSTANT (INTEGER 4 0 0 INTEGER ()) 0 '0')) (339 'bc' (
DERIVED 340 0 0 DERIVED ()) () (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0 POINTER) UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN
()) 0)) (341 'name' (CHARACTER 1 0 0 CHARACTER ((CONSTANT (INTEGER 4 0 0
INTEGER ()) 0 '101'))) () (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0) UNKNOWN-ACCESS ()) (342 'option_path' (CHARACTER 1
0 0 CHARACTER ((CONSTANT (INTEGER 4 0 0 INTEGER ()) 0 '8192'))) () (
UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0)
UNKNOWN-ACCESS (CONSTANT (CHARACTER 1 0 0 CHARACTER (())) 0 8192
'/uninitialised_path/                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            '))
(343 'mesh' (DERIVED 71 0 0 DERIVED ()) () (UNKNOWN-FL UNKNOWN-INTENT
UNKNOWN-PROC UNKNOWN UNKNOWN 0 0) UNKNOWN-ACCESS (STRUCTURE (DERIVED 71
0 0 DERIVED ()) 0 ((() ()) ((CONSTANT (LOGICAL 4 0 0 LOGICAL ()) 0 1) ())
((STRUCTURE (DERIVED 88 0 0 DERIVED ()) 0 ((() ()) (() ()) (() ()) (() ())
((NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0) ()) ((NULL (UNKNOWN 0 0 0 UNKNOWN ())
0) ()) ((NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0) ()) ((NULL (UNKNOWN 0 0 0
UNKNOWN ()) 0) ()) ((NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0) ()) ((NULL (
UNKNOWN 0 0 0 UNKNOWN ()) 0) ()) ((NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0) ())
((STRUCTURE (DERIVED 89 0 0 DERIVED ()) 0 ((() ()) (() ()) (() ()) (() ())
((NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0) ()) ((NULL (UNKNOWN 0 0 0 UNKNOWN ())
0) ()) (() ()) ((NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0) ()) (() ())) ()) ())
((NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0) ()) ((NULL (UNKNOWN 0 0 0 UNKNOWN ())
0) ()) ((NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0) ()) ((NULL (UNKNOWN 0 0 0
UNKNOWN ()) 0) ()) (() ())) ()) ()) (() ()) (() ()) (() ()) ((CONSTANT (
CHARACTER 1 0 0 CHARACTER (())) 0 8192
'/uninitialised_path/                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            ')
()) ((CONSTANT (INTEGER 4 0 0 INTEGER ()) 0 '0') ()) ((NULL (UNKNOWN 0 0
0 UNKNOWN ()) 0) ()) ((NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0) ()) ((NULL (
UNKNOWN 0 0 0 UNKNOWN ()) 0) ()) ((NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0) ())
((NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0) ()) ((NULL (UNKNOWN 0 0 0 UNKNOWN ())
0) ()) ((NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0) ()) ((NULL (UNKNOWN 0 0 0
UNKNOWN ()) 0) ()) ((NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0) ()) ((CONSTANT (
LOGICAL 4 0 0 LOGICAL ()) 0 0) ())) ())) (344 'refcount' (DERIVED 96 0 0
DERIVED ()) () (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0
0 POINTER) UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0)) (345
'aliased' (LOGICAL 4 0 0 LOGICAL ()) () (UNKNOWN-FL UNKNOWN-INTENT
UNKNOWN-PROC UNKNOWN UNKNOWN 0 0) UNKNOWN-ACCESS (CONSTANT (LOGICAL 4 0
0 LOGICAL ()) 0 0)) (346 'py_locweight' (REAL 8 0 0 REAL ()) (2 0
DEFERRED () () () ()) (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN
UNKNOWN 0 0 DIMENSION POINTER) UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0
UNKNOWN ()) 0)) (347 'py_func' (CHARACTER 1 0 0 CHARACTER ((CONSTANT (
INTEGER 4 0 0 INTEGER ()) 0 '8192'))) () (UNKNOWN-FL UNKNOWN-INTENT
UNKNOWN-PROC UNKNOWN UNKNOWN 0 0) UNKNOWN-ACCESS ()) (348 'py_positions'
(DERIVED 194 0 0 DERIVED ()) () (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0 POINTER) UNKNOWN-ACCESS ()) (349
'py_positions_same_mesh' (LOGICAL 4 0 0 LOGICAL ()) () (UNKNOWN-FL
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0) UNKNOWN-ACCESS ()) (
350 'py_dim' (INTEGER 4 0 0 INTEGER ()) () (UNKNOWN-FL UNKNOWN-INTENT
UNKNOWN-PROC UNKNOWN UNKNOWN 0 0) UNKNOWN-ACCESS ()) (351
'py_positions_shape' (DERIVED 88 0 0 DERIVED ()) () (UNKNOWN-FL
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 POINTER) UNKNOWN-ACCESS
(NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0))) PUBLIC (() () () ()) () 0 0
64912956)
194 'vector_field' 'fields_data_types' 'vector_field' 1 ((DERIVED
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 POINTER_COMP) (UNKNOWN 0
0 0 UNKNOWN ()) 0 0 () () 0 ((352 'val' (REAL 8 0 0 REAL ()) (2 0
DEFERRED () () () ()) (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN
UNKNOWN 0 0 DIMENSION POINTER) UNKNOWN-ACCESS ()) (353 'wrapped' (
LOGICAL 4 0 0 LOGICAL ()) () (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0) UNKNOWN-ACCESS (CONSTANT (LOGICAL 4 0 0 LOGICAL ())
0 1)) (354 'field_type' (INTEGER 4 0 0 INTEGER ()) () (UNKNOWN-FL
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0) UNKNOWN-ACCESS (
CONSTANT (INTEGER 4 0 0 INTEGER ()) 0 '0')) (355 'bc' (DERIVED 356 0 0
DERIVED ()) () (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0
0 POINTER) UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0)) (357 'name'
(CHARACTER 1 0 0 CHARACTER ((CONSTANT (INTEGER 4 0 0 INTEGER ()) 0 '101')))
() (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0)
UNKNOWN-ACCESS ()) (358 'dim' (INTEGER 4 0 0 INTEGER ()) () (UNKNOWN-FL
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0) UNKNOWN-ACCESS ()) (
359 'option_path' (CHARACTER 1 0 0 CHARACTER ((CONSTANT (INTEGER 4 0 0
INTEGER ()) 0 '8192'))) () (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0) UNKNOWN-ACCESS (CONSTANT (CHARACTER 1 0 0 CHARACTER
(())) 0 8192
'/uninitialised_path/                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            '))
(360 'mesh' (DERIVED 71 0 0 DERIVED ()) () (UNKNOWN-FL UNKNOWN-INTENT
UNKNOWN-PROC UNKNOWN UNKNOWN 0 0) UNKNOWN-ACCESS (STRUCTURE (DERIVED 71
0 0 DERIVED ()) 0 ((() ()) ((CONSTANT (LOGICAL 4 0 0 LOGICAL ()) 0 1) ())
((STRUCTURE (DERIVED 88 0 0 DERIVED ()) 0 ((() ()) (() ()) (() ()) (() ())
((NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0) ()) ((NULL (UNKNOWN 0 0 0 UNKNOWN ())
0) ()) ((NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0) ()) ((NULL (UNKNOWN 0 0 0
UNKNOWN ()) 0) ()) ((NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0) ()) ((NULL (
UNKNOWN 0 0 0 UNKNOWN ()) 0) ()) ((NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0) ())
((STRUCTURE (DERIVED 89 0 0 DERIVED ()) 0 ((() ()) (() ()) (() ()) (() ())
((NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0) ()) ((NULL (UNKNOWN 0 0 0 UNKNOWN ())
0) ()) (() ()) ((NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0) ()) (() ())) ()) ())
((NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0) ()) ((NULL (UNKNOWN 0 0 0 UNKNOWN ())
0) ()) ((NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0) ()) ((NULL (UNKNOWN 0 0 0
UNKNOWN ()) 0) ()) (() ())) ()) ()) (() ()) (() ()) (() ()) ((CONSTANT (
CHARACTER 1 0 0 CHARACTER (())) 0 8192
'/uninitialised_path/                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            ')
()) ((CONSTANT (INTEGER 4 0 0 INTEGER ()) 0 '0') ()) ((NULL (UNKNOWN 0 0
0 UNKNOWN ()) 0) ()) ((NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0) ()) ((NULL (
UNKNOWN 0 0 0 UNKNOWN ()) 0) ()) ((NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0) ())
((NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0) ()) ((NULL (UNKNOWN 0 0 0 UNKNOWN ())
0) ()) ((NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0) ()) ((NULL (UNKNOWN 0 0 0
UNKNOWN ()) 0) ()) ((NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0) ()) ((CONSTANT (
LOGICAL 4 0 0 LOGICAL ()) 0 0) ())) ())) (361 'refcount' (DERIVED 96 0 0
DERIVED ()) () (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0
0 POINTER) UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0)) (362
'aliased' (LOGICAL 4 0 0 LOGICAL ()) () (UNKNOWN-FL UNKNOWN-INTENT
UNKNOWN-PROC UNKNOWN UNKNOWN 0 0) UNKNOWN-ACCESS (CONSTANT (LOGICAL 4 0
0 LOGICAL ()) 0 0)) (363 'picker' (DERIVED 364 0 0 DERIVED ()) () (
UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 POINTER)
UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0))) PUBLIC (() () () ())
() 0 0 43598963)
324 'csr_matrix' 'sparse_tools' 'csr_matrix' 1 ((DERIVED UNKNOWN-INTENT
UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 POINTER_COMP) (UNKNOWN 0 0 0 UNKNOWN ())
0 0 () () 0 ((365 'sparsity' (DERIVED 319 0 0 DERIVED ()) () (
UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0)
UNKNOWN-ACCESS (STRUCTURE (DERIVED 319 0 0 DERIVED ()) 0 (((NULL (
UNKNOWN 0 0 0 UNKNOWN ()) 0) ()) ((NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0) ())
((NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0) ()) (() ()) ((NULL (UNKNOWN 0 0 0
UNKNOWN ()) 0) ()) ((NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0) ()) ((NULL (
UNKNOWN 0 0 0 UNKNOWN ()) 0) ()) ((CONSTANT (CHARACTER 1 0 0 CHARACTER (
())) 0 101
'                                                                                                     ')
()) ((CONSTANT (LOGICAL 4 0 0 LOGICAL ()) 0 0) ()) ((CONSTANT (LOGICAL 4
0 0 LOGICAL ()) 0 0) ())) ())) (366 'val' (REAL 8 0 0 REAL ()) (1 0
DEFERRED () ()) (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN
0 0 DIMENSION POINTER) UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0))
(367 'ival' (INTEGER 4 0 0 INTEGER ()) (1 0 DEFERRED () ()) (UNKNOWN-FL
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 DIMENSION POINTER)
UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0)) (368 'clone' (
LOGICAL 4 0 0 LOGICAL ()) () (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0) UNKNOWN-ACCESS (CONSTANT (LOGICAL 4 0 0 LOGICAL ())
0 0)) (369 'external_val' (LOGICAL 4 0 0 LOGICAL ()) () (UNKNOWN-FL
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0) UNKNOWN-ACCESS (
CONSTANT (LOGICAL 4 0 0 LOGICAL ()) 0 0)) (370 'inactive' (DERIVED 371 0
0 DERIVED ()) () (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN
0 0 POINTER) UNKNOWN-ACCESS ()) (372 'ksp' (INTEGER 8 0 0 INTEGER ()) ()
(UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 POINTER)
UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0)) (373 'refcount' (
DERIVED 96 0 0 DERIVED ()) () (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0 POINTER) UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN
()) 0)) (374 'name' (CHARACTER 1 0 0 CHARACTER ((CONSTANT (INTEGER 4 0 0
INTEGER ()) 0 '101'))) () (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0) UNKNOWN-ACCESS (CONSTANT (CHARACTER 1 0 0 CHARACTER
(())) 0 101
'                                                                                                     ')))
PUBLIC (() () () ()) () 0 0 66371681)
371 'logical_array_ptr' 'sparse_tools' 'logical_array_ptr' 1 ((DERIVED
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 POINTER_COMP) (UNKNOWN 0
0 0 UNKNOWN ()) 0 0 () () 0 ((375 'ptr' (LOGICAL 4 0 0 LOGICAL ()) (1 0
DEFERRED () ()) (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN
0 0 DIMENSION POINTER) UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0)))
PUBLIC (() () () ()) () 0 0 30974511)
119 'polynomial' 'polynomials' 'polynomial' 1 ((DERIVED UNKNOWN-INTENT
UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 POINTER_COMP) (UNKNOWN 0 0 0 UNKNOWN ())
0 0 () () 0 ((376 'coefs' (REAL 8 0 0 REAL ()) (1 0 DEFERRED () ()) (
UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 DIMENSION
POINTER) UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0)) (377 'degree'
(INTEGER 4 0 0 INTEGER ()) () (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0) UNKNOWN-ACCESS (CONSTANT (INTEGER 4 0 0 INTEGER ())
0 '-1'))) PUBLIC (() () () ()) () 0 0 87989236)
54 'petsc_comm_self' 'petscsysdef' 'petsc_comm_self' 1 ((VARIABLE
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 IN_COMMON) (INTEGER 4 0
0 INTEGER ()) 0 0 () () 0 () () () 0 0)
319 'csr_sparsity' 'sparse_tools' 'csr_sparsity' 1 ((DERIVED
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 POINTER_COMP) (UNKNOWN 0
0 0 UNKNOWN ()) 0 0 () () 0 ((378 'findrm' (INTEGER 4 0 0 INTEGER ()) (
1 0 DEFERRED () ()) (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN
UNKNOWN 0 0 DIMENSION POINTER) UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0
UNKNOWN ()) 0)) (379 'centrm' (INTEGER 4 0 0 INTEGER ()) (1 0 DEFERRED ()
()) (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DIMENSION POINTER) UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0)) (
380 'colm' (INTEGER 4 0 0 INTEGER ()) (1 0 DEFERRED () ()) (UNKNOWN-FL
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 DIMENSION POINTER)
UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0)) (381 'columns' (
INTEGER 4 0 0 INTEGER ()) () (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0) UNKNOWN-ACCESS ()) (382 'row_halo' (DERIVED 107 0 0
DERIVED ()) () (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0
0 POINTER) UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0)) (383
'column_halo' (DERIVED 107 0 0 DERIVED ()) () (UNKNOWN-FL UNKNOWN-INTENT
UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 POINTER) UNKNOWN-ACCESS (NULL (UNKNOWN
0 0 0 UNKNOWN ()) 0)) (384 'refcount' (DERIVED 96 0 0 DERIVED ()) () (
UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 POINTER)
UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0)) (385 'name' (
CHARACTER 1 0 0 CHARACTER ((CONSTANT (INTEGER 4 0 0 INTEGER ()) 0 '101')))
() (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0)
UNKNOWN-ACCESS (CONSTANT (CHARACTER 1 0 0 CHARACTER (())) 0 101
'                                                                                                     '))
(386 'wrapped' (LOGICAL 4 0 0 LOGICAL ()) () (UNKNOWN-FL UNKNOWN-INTENT
UNKNOWN-PROC UNKNOWN UNKNOWN 0 0) UNKNOWN-ACCESS (CONSTANT (LOGICAL 4 0
0 LOGICAL ()) 0 0)) (387 'sorted_rows' (LOGICAL 4 0 0 LOGICAL ()) () (
UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0)
UNKNOWN-ACCESS (CONSTANT (LOGICAL 4 0 0 LOGICAL ()) 0 0))) PUBLIC (() ()
() ()) () 0 0 78882345)
89 'quadrature_type' 'quadrature' 'quadrature_type' 1 ((DERIVED
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 POINTER_COMP) (UNKNOWN 0
0 0 UNKNOWN ()) 0 0 () () 0 ((388 'dim' (INTEGER 4 0 0 INTEGER ()) () (
UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0)
UNKNOWN-ACCESS ()) (389 'degree' (INTEGER 4 0 0 INTEGER ()) () (
UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0)
UNKNOWN-ACCESS ()) (390 'vertices' (INTEGER 4 0 0 INTEGER ()) () (
UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0)
UNKNOWN-ACCESS ()) (391 'ngi' (INTEGER 4 0 0 INTEGER ()) () (UNKNOWN-FL
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0) UNKNOWN-ACCESS ()) (
392 'weight' (REAL 8 0 0 REAL ()) (1 0 DEFERRED () ()) (UNKNOWN-FL
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 DIMENSION POINTER)
UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0)) (393 'l' (REAL 8 0 0
REAL ()) (2 0 DEFERRED () () () ()) (UNKNOWN-FL UNKNOWN-INTENT
UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 DIMENSION POINTER) UNKNOWN-ACCESS (
NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0)) (394 'name' (CHARACTER 1 0 0
CHARACTER ((CONSTANT (INTEGER 4 0 0 INTEGER ()) 0 '0'))) () (UNKNOWN-FL
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0) UNKNOWN-ACCESS ()) (
395 'refcount' (DERIVED 96 0 0 DERIVED ()) () (UNKNOWN-FL UNKNOWN-INTENT
UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 POINTER) UNKNOWN-ACCESS (NULL (UNKNOWN
0 0 0 UNKNOWN ()) 0)) (396 'family' (INTEGER 4 0 0 INTEGER ()) () (
UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0)
UNKNOWN-ACCESS ())) PUBLIC (() () () ()) () 0 0 59837722)
122 'ele_numbering_type' 'element_numbering' 'ele_numbering_type' 1 ((
DERIVED UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 POINTER_COMP) (
UNKNOWN 0 0 0 UNKNOWN ()) 0 0 () () 0 ((397 'faces' (INTEGER 4 0 0
INTEGER ()) () (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0
0) UNKNOWN-ACCESS ()) (398 'vertices' (INTEGER 4 0 0 INTEGER ()) () (
UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0)
UNKNOWN-ACCESS ()) (399 'edges' (INTEGER 4 0 0 INTEGER ()) () (
UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0)
UNKNOWN-ACCESS ()) (400 'boundaries' (INTEGER 4 0 0 INTEGER ()) () (
UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0)
UNKNOWN-ACCESS ()) (401 'degree' (INTEGER 4 0 0 INTEGER ()) () (
UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0)
UNKNOWN-ACCESS ()) (402 'dimension' (INTEGER 4 0 0 INTEGER ()) () (
UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0)
UNKNOWN-ACCESS ()) (403 'nodes' (INTEGER 4 0 0 INTEGER ()) () (
UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0)
UNKNOWN-ACCESS ()) (404 'type' (INTEGER 4 0 0 INTEGER ()) () (
UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0)
UNKNOWN-ACCESS (CONSTANT (INTEGER 4 0 0 INTEGER ()) 0 '1')) (405 'family'
(INTEGER 4 0 0 INTEGER ()) () (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0) UNKNOWN-ACCESS ()) (406 'count2number' (INTEGER 4 0
0 INTEGER ()) (3 0 DEFERRED () () () () () ()) (UNKNOWN-FL
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 DIMENSION POINTER)
UNKNOWN-ACCESS ()) (407 'number2count' (INTEGER 4 0 0 INTEGER ()) (2 0
DEFERRED () () () ()) (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN
UNKNOWN 0 0 DIMENSION POINTER) UNKNOWN-ACCESS ()) (408 'boundary_coord'
(INTEGER 4 0 0 INTEGER ()) (1 0 DEFERRED () ()) (UNKNOWN-FL
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 DIMENSION POINTER)
UNKNOWN-ACCESS ()) (409 'boundary_val' (INTEGER 4 0 0 INTEGER ()) (1 0
DEFERRED () ()) (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN
0 0 DIMENSION POINTER) UNKNOWN-ACCESS ())) PUBLIC (() () () ()) () 0 0
96431082)
126 'superconvergence_type' 'elements' 'superconvergence_type' 1 ((
DERIVED UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 POINTER_COMP) (
UNKNOWN 0 0 0 UNKNOWN ()) 0 0 () () 0 ((410 'nsp' (INTEGER 4 0 0 INTEGER
()) () (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0)
UNKNOWN-ACCESS ()) (411 'l' (REAL 8 0 0 REAL ()) (2 0 DEFERRED () () ()
()) (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DIMENSION POINTER) UNKNOWN-ACCESS ()) (412 'n' (REAL 8 0 0 REAL ()) (2 0
DEFERRED () () () ()) (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN
UNKNOWN 0 0 DIMENSION POINTER) UNKNOWN-ACCESS ()) (413 'dn' (REAL 8 0 0
REAL ()) (3 0 DEFERRED () () () () () ()) (UNKNOWN-FL UNKNOWN-INTENT
UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 DIMENSION POINTER) UNKNOWN-ACCESS ()))
PUBLIC (() () () ()) () 0 0 18282395)
128 'constraints_type' 'elements' 'constraints_type' 1 ((DERIVED
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 POINTER_COMP) (UNKNOWN 0
0 0 UNKNOWN ()) 0 0 () () 0 ((414 'type' (INTEGER 4 0 0 INTEGER ()) () (
UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0)
UNKNOWN-ACCESS ()) (415 'dim' (INTEGER 4 0 0 INTEGER ()) () (UNKNOWN-FL
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0) UNKNOWN-ACCESS ()) (
416 'degree' (INTEGER 4 0 0 INTEGER ()) () (UNKNOWN-FL UNKNOWN-INTENT
UNKNOWN-PROC UNKNOWN UNKNOWN 0 0) UNKNOWN-ACCESS ()) (417 'loc' (
INTEGER 4 0 0 INTEGER ()) () (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0) UNKNOWN-ACCESS ()) (418 'n_constraints' (INTEGER 4
0 0 INTEGER ()) () (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN
UNKNOWN 0 0) UNKNOWN-ACCESS ()) (419 'orthogonal' (REAL 8 0 0 REAL ()) (
3 0 DEFERRED () () () () () ()) (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0 DIMENSION POINTER) UNKNOWN-ACCESS (NULL (UNKNOWN 0 0
0 UNKNOWN ()) 0))) PUBLIC (() () () ()) () 0 0 57548971)
2 'add_vector_boundary_condition_surface_elements' 'boundary_conditions'
'add_vector_boundary_condition_surface_elements' 1 ((PROCEDURE
UNKNOWN-INTENT MODULE-PROC DECL UNKNOWN 0 0 SUBROUTINE ALWAYS_EXPLICIT)
(UNKNOWN 0 0 0 UNKNOWN ()) 420 0 (421 422 423 424 425 426) () 0 () () ()
0 0)
22 'get_scalar_boundary_condition_by_number' 'boundary_conditions'
'get_scalar_boundary_condition_by_number' 1 ((PROCEDURE UNKNOWN-INTENT
MODULE-PROC DECL UNKNOWN 0 0 SUBROUTINE ALWAYS_EXPLICIT) (UNKNOWN 0 0 0
UNKNOWN ()) 427 0 (428 429 430 431 432 433 434 435) () 0 () () () 0 0)
21 'get_vector_boundary_condition_by_number' 'boundary_conditions'
'get_vector_boundary_condition_by_number' 1 ((PROCEDURE UNKNOWN-INTENT
MODULE-PROC DECL UNKNOWN 0 0 SUBROUTINE ALWAYS_EXPLICIT) (UNKNOWN 0 0 0
UNKNOWN ()) 436 0 (437 438 439 440 441 442 443 444 445) () 0 () () () 0
0)
421 'field' '' 'field' 420 ((VARIABLE INOUT UNKNOWN-PROC UNKNOWN UNKNOWN
0 0 DUMMY) (DERIVED 194 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
422 'name' '' 'name' 420 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DUMMY) (CHARACTER 1 0 0 CHARACTER (())) 0 0 () () 0 () () () 0 0)
423 'type' '' 'type' 420 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DUMMY) (CHARACTER 1 0 0 CHARACTER (())) 0 0 () () 0 () () () 0 0)
424 'surface_element_list' '' 'surface_element_list' 420 ((VARIABLE IN
UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 DIMENSION DUMMY) (INTEGER 4 0 0 INTEGER
()) 0 0 () (1 0 ASSUMED_SHAPE (CONSTANT (INTEGER 4 0 0 INTEGER ()) 0 '1')
()) 0 () () () 0 0)
425 'applies' '' 'applies' 420 ((VARIABLE IN UNKNOWN-PROC UNKNOWN
UNKNOWN 0 0 DIMENSION OPTIONAL DUMMY) (LOGICAL 4 0 0 LOGICAL ()) 0 0 ()
(1 0 ASSUMED_SHAPE (CONSTANT (INTEGER 4 0 0 INTEGER ()) 0 '1') ()) 0 ()
() () 0 0)
426 'option_path' '' 'option_path' 420 ((VARIABLE IN UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0 OPTIONAL DUMMY) (CHARACTER 1 0 0 CHARACTER (())) 0 0
() () 0 () () () 0 0)
183 'surface_field' '' 'surface_field' 180 ((VARIABLE IN UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0 DUMMY) (DERIVED 188 0 0 DERIVED ()) 0 0 () () 0 () ()
() 0 0)
173 'field' '' 'field' 172 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0
0 DUMMY) (DERIVED 188 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
185 'field' '' 'field' 184 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0
0 DUMMY) (DERIVED 194 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
174 'name' '' 'name' 172 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DUMMY) (CHARACTER 1 0 0 CHARACTER (())) 0 0 () () 0 () () () 0 0)
175 'surface_field' '' 'surface_field' 172 ((VARIABLE IN UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0 DUMMY) (DERIVED 188 0 0 DERIVED ()) 0 0 () () 0 () ()
() 0 0)
177 'field' '' 'field' 176 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0
0 DUMMY) (DERIVED 194 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
178 'name' '' 'name' 176 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DUMMY) (CHARACTER 1 0 0 CHARACTER (())) 0 0 () () 0 () () () 0 0)
179 'surface_field' '' 'surface_field' 176 ((VARIABLE IN UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0 DUMMY) (DERIVED 194 0 0 DERIVED ()) 0 0 () () 0 () ()
() 0 0)
186 'name' '' 'name' 184 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DUMMY) (CHARACTER 1 0 0 CHARACTER (())) 0 0 () () 0 () () () 0 0)
187 'surface_field' '' 'surface_field' 184 ((VARIABLE IN UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0 DUMMY) (DERIVED 188 0 0 DERIVED ()) 0 0 () () 0 () ()
() 0 0)
190 'field' '' 'field' 189 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0
0 DUMMY) (DERIVED 188 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
191 'n' '' 'n' 189 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 DUMMY)
(INTEGER 4 0 0 INTEGER ()) 0 0 () () 0 () () () 0 0)
192 'name' '' 'name' 189 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DUMMY) (CHARACTER 1 0 0 CHARACTER (())) 0 0 () () 0 () () () 0 0)
193 'surface_field' '' 'surface_field' 189 ((VARIABLE UNKNOWN-INTENT
UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 POINTER RESULT ALWAYS_EXPLICIT) (
DERIVED 188 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
196 'field' '' 'field' 195 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0
0 DUMMY) (DERIVED 194 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
197 'n' '' 'n' 195 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 DUMMY)
(INTEGER 4 0 0 INTEGER ()) 0 0 () () 0 () () () 0 0)
198 'name' '' 'name' 195 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DUMMY) (CHARACTER 1 0 0 CHARACTER (())) 0 0 () () 0 () () () 0 0)
199 'surface_field' '' 'surface_field' 195 ((VARIABLE UNKNOWN-INTENT
UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 POINTER RESULT ALWAYS_EXPLICIT) (
DERIVED 194 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
211 'field' '' 'field' 210 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0
0 DUMMY) (DERIVED 194 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
212 'n' '' 'n' 210 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 DUMMY)
(INTEGER 4 0 0 INTEGER ()) 0 0 () () 0 () () () 0 0)
213 'name' '' 'name' 210 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DUMMY) (CHARACTER 1 0 0 CHARACTER (())) 0 0 () () 0 () () () 0 0)
214 'surface_field' '' 'surface_field' 210 ((VARIABLE UNKNOWN-INTENT
UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 POINTER RESULT ALWAYS_EXPLICIT) (
DERIVED 188 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
201 'field' '' 'field' 200 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0
0 DUMMY) (DERIVED 188 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
202 'bc_name' '' 'bc_name' 200 ((VARIABLE IN UNKNOWN-PROC UNKNOWN
UNKNOWN 0 0 DUMMY) (CHARACTER 1 0 0 CHARACTER (())) 0 0 () () 0 () () ()
0 0)
203 'name' '' 'name' 200 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DUMMY) (CHARACTER 1 0 0 CHARACTER (())) 0 0 () () 0 () () () 0 0)
204 'surface_field' '' 'surface_field' 200 ((VARIABLE UNKNOWN-INTENT
UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 POINTER RESULT ALWAYS_EXPLICIT) (
DERIVED 188 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
206 'field' '' 'field' 205 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0
0 DUMMY) (DERIVED 194 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
207 'bc_name' '' 'bc_name' 205 ((VARIABLE IN UNKNOWN-PROC UNKNOWN
UNKNOWN 0 0 DUMMY) (CHARACTER 1 0 0 CHARACTER (())) 0 0 () () 0 () () ()
0 0)
208 'name' '' 'name' 205 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DUMMY) (CHARACTER 1 0 0 CHARACTER (())) 0 0 () () 0 () () () 0 0)
209 'surface_field' '' 'surface_field' 205 ((VARIABLE UNKNOWN-INTENT
UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 POINTER RESULT ALWAYS_EXPLICIT) (
DERIVED 194 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
216 'field' '' 'field' 215 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0
0 DUMMY) (DERIVED 194 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
217 'bc_name' '' 'bc_name' 215 ((VARIABLE IN UNKNOWN-PROC UNKNOWN
UNKNOWN 0 0 DUMMY) (CHARACTER 1 0 0 CHARACTER (())) 0 0 () () 0 () () ()
0 0)
218 'name' '' 'name' 215 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DUMMY) (CHARACTER 1 0 0 CHARACTER (())) 0 0 () () 0 () () () 0 0)
219 'surface_field' '' 'surface_field' 215 ((VARIABLE UNKNOWN-INTENT
UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 POINTER RESULT ALWAYS_EXPLICIT) (
DERIVED 188 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
221 'field' '' 'field' 220 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0
0 DUMMY) (DERIVED 188 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
222 'n' '' 'n' 220 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 DUMMY)
(INTEGER 4 0 0 INTEGER ()) 0 0 () () 0 () () () 0 0)
223 'name' '' 'name' 220 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DUMMY) (CHARACTER 1 0 0 CHARACTER (())) 0 0 () () 0 () () () 0 0)
225 'field' '' 'field' 224 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0
0 DUMMY) (DERIVED 194 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
226 'n' '' 'n' 224 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 DUMMY)
(INTEGER 4 0 0 INTEGER ()) 0 0 () () 0 () () () 0 0)
227 'name' '' 'name' 224 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DUMMY) (CHARACTER 1 0 0 CHARACTER (())) 0 0 () () 0 () () () 0 0)
80 'field' '' 'field' 79 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DUMMY) (DERIVED 194 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
81 'n' '' 'n' 79 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 DUMMY) (
INTEGER 4 0 0 INTEGER ()) 0 0 () () 0 () () () 0 0)
82 'name' '' 'name' 79 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DUMMY) (CHARACTER 1 0 0 CHARACTER (())) 0 0 () () 0 () () () 0 0)
229 'field' '' 'field' 228 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0
0 DUMMY) (DERIVED 188 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
231 'field' '' 'field' 230 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0
0 DUMMY) (DERIVED 194 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
428 'field' '' 'field' 427 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0
0 DUMMY) (DERIVED 188 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
429 'n' '' 'n' 427 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 DUMMY)
(INTEGER 4 0 0 INTEGER ()) 0 0 () () 0 () () () 0 0)
430 'name' '' 'name' 427 ((VARIABLE OUT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
OPTIONAL DUMMY) (CHARACTER 1 0 0 CHARACTER (())) 0 0 () () 0 () () () 0
0)
431 'type' '' 'type' 427 ((VARIABLE OUT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
OPTIONAL DUMMY) (CHARACTER 1 0 0 CHARACTER (())) 0 0 () () 0 () () () 0
0)
432 'surface_element_list' '' 'surface_element_list' 427 ((VARIABLE
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 DIMENSION OPTIONAL
POINTER DUMMY) (INTEGER 4 0 0 INTEGER ()) 0 0 () (1 0 DEFERRED () ()) 0
() () () 0 0)
433 'surface_node_list' '' 'surface_node_list' 427 ((VARIABLE
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 DIMENSION OPTIONAL
POINTER DUMMY) (INTEGER 4 0 0 INTEGER ()) 0 0 () (1 0 DEFERRED () ()) 0
() () () 0 0)
434 'surface_mesh' '' 'surface_mesh' 427 ((VARIABLE UNKNOWN-INTENT
UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 OPTIONAL POINTER DUMMY) (DERIVED 71 0 0
DERIVED ()) 0 0 () () 0 () () () 0 0)
435 'option_path' '' 'option_path' 427 ((VARIABLE OUT UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0 OPTIONAL DUMMY) (CHARACTER 1 0 0 CHARACTER (())) 0 0
() () 0 () () () 0 0)
437 'field' '' 'field' 436 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0
0 DUMMY) (DERIVED 194 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
438 'n' '' 'n' 436 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 DUMMY)
(INTEGER 4 0 0 INTEGER ()) 0 0 () () 0 () () () 0 0)
439 'name' '' 'name' 436 ((VARIABLE OUT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
OPTIONAL DUMMY) (CHARACTER 1 0 0 CHARACTER (())) 0 0 () () 0 () () () 0
0)
440 'type' '' 'type' 436 ((VARIABLE OUT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
OPTIONAL DUMMY) (CHARACTER 1 0 0 CHARACTER (())) 0 0 () () 0 () () () 0
0)
441 'surface_element_list' '' 'surface_element_list' 436 ((VARIABLE
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 DIMENSION OPTIONAL
POINTER DUMMY) (INTEGER 4 0 0 INTEGER ()) 0 0 () (1 0 DEFERRED () ()) 0
() () () 0 0)
442 'surface_node_list' '' 'surface_node_list' 436 ((VARIABLE
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 DIMENSION OPTIONAL
POINTER DUMMY) (INTEGER 4 0 0 INTEGER ()) 0 0 () (1 0 DEFERRED () ()) 0
() () () 0 0)
443 'applies' '' 'applies' 436 ((VARIABLE OUT UNKNOWN-PROC UNKNOWN
UNKNOWN 0 0 DIMENSION OPTIONAL DUMMY) (LOGICAL 4 0 0 LOGICAL ()) 0 0 ()
(1 0 ASSUMED_SHAPE (CONSTANT (INTEGER 4 0 0 INTEGER ()) 0 '1') ()) 0 ()
() () 0 0)
444 'surface_mesh' '' 'surface_mesh' 436 ((VARIABLE UNKNOWN-INTENT
UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 OPTIONAL POINTER DUMMY) (DERIVED 71 0 0
DERIVED ()) 0 0 () () 0 () () () 0 0)
445 'option_path' '' 'option_path' 436 ((VARIABLE OUT UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0 OPTIONAL DUMMY) (CHARACTER 1 0 0 CHARACTER (())) 0 0
() () 0 () () () 0 0)
148 'field' '' 'field' 147 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0
0 DUMMY) (DERIVED 188 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
149 'name' '' 'name' 147 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DUMMY) (CHARACTER 1 0 0 CHARACTER (())) 0 0 () () 0 () () () 0 0)
150 'type' '' 'type' 147 ((VARIABLE OUT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
OPTIONAL DUMMY) (CHARACTER 1 0 0 CHARACTER (())) 0 0 () () 0 () () () 0
0)
151 'surface_node_list' '' 'surface_node_list' 147 ((VARIABLE
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 DIMENSION OPTIONAL
POINTER DUMMY) (INTEGER 4 0 0 INTEGER ()) 0 0 () (1 0 DEFERRED () ()) 0
() () () 0 0)
152 'surface_element_list' '' 'surface_element_list' 147 ((VARIABLE
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 DIMENSION OPTIONAL
POINTER DUMMY) (INTEGER 4 0 0 INTEGER ()) 0 0 () (1 0 DEFERRED () ()) 0
() () () 0 0)
153 'surface_mesh' '' 'surface_mesh' 147 ((VARIABLE UNKNOWN-INTENT
UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 OPTIONAL POINTER DUMMY) (DERIVED 71 0 0
DERIVED ()) 0 0 () () 0 () () () 0 0)
154 'option_path' '' 'option_path' 147 ((VARIABLE OUT UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0 OPTIONAL DUMMY) (CHARACTER 1 0 0 CHARACTER (())) 0 0
() () 0 () () () 0 0)
156 'field' '' 'field' 155 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0
0 DUMMY) (DERIVED 194 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
157 'name' '' 'name' 155 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DUMMY) (CHARACTER 1 0 0 CHARACTER (())) 0 0 () () 0 () () () 0 0)
158 'type' '' 'type' 155 ((VARIABLE OUT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
OPTIONAL DUMMY) (CHARACTER 1 0 0 CHARACTER (())) 0 0 () () 0 () () () 0
0)
159 'surface_element_list' '' 'surface_element_list' 155 ((VARIABLE
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 DIMENSION OPTIONAL
POINTER DUMMY) (INTEGER 4 0 0 INTEGER ()) 0 0 () (1 0 DEFERRED () ()) 0
() () () 0 0)
160 'surface_node_list' '' 'surface_node_list' 155 ((VARIABLE
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 DIMENSION OPTIONAL
POINTER DUMMY) (INTEGER 4 0 0 INTEGER ()) 0 0 () (1 0 DEFERRED () ()) 0
() () () 0 0)
161 'applies' '' 'applies' 155 ((VARIABLE OUT UNKNOWN-PROC UNKNOWN
UNKNOWN 0 0 DIMENSION OPTIONAL DUMMY) (LOGICAL 4 0 0 LOGICAL ()) 0 0 ()
(1 0 ASSUMED_SHAPE (CONSTANT (INTEGER 4 0 0 INTEGER ()) 0 '1') ()) 0 ()
() () 0 0)
162 'surface_mesh' '' 'surface_mesh' 155 ((VARIABLE UNKNOWN-INTENT
UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 OPTIONAL POINTER DUMMY) (DERIVED 71 0 0
DERIVED ()) 0 0 () () 0 () () () 0 0)
163 'option_path' '' 'option_path' 155 ((VARIABLE OUT UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0 OPTIONAL DUMMY) (CHARACTER 1 0 0 CHARACTER (())) 0 0
() () 0 () () () 0 0)
76 'mesh' '' 'mesh' 75 ((VARIABLE INOUT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DUMMY) (DERIVED 71 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
77 'periodic_bc_list' '' 'periodic_bc_list' 75 ((VARIABLE OUT
UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 DIMENSION DUMMY) (LOGICAL 4 0 0 LOGICAL
()) 0 0 () (1 0 ASSUMED_SHAPE (CONSTANT (INTEGER 4 0 0 INTEGER ()) 0 '1')
()) 0 () () () 0 0)
233 'field' '' 'field' 232 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0
0 TARGET DUMMY) (DERIVED 188 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
234 'types' '' 'types' 232 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0
0 DIMENSION DUMMY) (CHARACTER 1 0 0 CHARACTER (())) 0 0 () (1 0
ASSUMED_SHAPE (CONSTANT (INTEGER 4 0 0 INTEGER ()) 0 '1') ()) 0 () () ()
0 0)
235 'boundary_value' '' 'boundary_value' 232 ((VARIABLE OUT UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0 DUMMY) (DERIVED 188 0 0 DERIVED ()) 0 0 () () 0 () ()
() 0 0)
236 'bc_type_list' '' 'bc_type_list' 232 ((VARIABLE OUT UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0 DIMENSION DUMMY) (INTEGER 4 0 0 INTEGER ()) 0 0 () (
1 0 ASSUMED_SHAPE (CONSTANT (INTEGER 4 0 0 INTEGER ()) 0 '1') ()) 0 () ()
() 0 0)
237 'bc_number_list' '' 'bc_number_list' 232 ((VARIABLE OUT UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0 DIMENSION OPTIONAL DUMMY) (INTEGER 4 0 0 INTEGER ())
0 0 () (1 0 ASSUMED_SHAPE (CONSTANT (INTEGER 4 0 0 INTEGER ()) 0 '1') ())
0 () () () 0 0)
238 'boundary_second_value' '' 'boundary_second_value' 232 ((VARIABLE
OUT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 OPTIONAL DUMMY) (DERIVED 188 0 0
DERIVED ()) 0 0 () () 0 () () () 0 0)
247 'types' '' 'types' 245 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0
0 DIMENSION DUMMY) (CHARACTER 1 0 0 CHARACTER (())) 0 0 () (1 0
ASSUMED_SHAPE (CONSTANT (INTEGER 4 0 0 INTEGER ()) 0 '1') ()) 0 () () ()
0 0)
240 'field' '' 'field' 239 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0
0 TARGET DUMMY) (DERIVED 194 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
241 'types' '' 'types' 239 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0
0 DIMENSION DUMMY) (CHARACTER 1 0 0 CHARACTER (())) 0 0 () (1 0
ASSUMED_SHAPE (CONSTANT (INTEGER 4 0 0 INTEGER ()) 0 '1') ()) 0 () () ()
0 0)
242 'boundary_value' '' 'boundary_value' 239 ((VARIABLE OUT UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0 DUMMY) (DERIVED 194 0 0 DERIVED ()) 0 0 () () 0 () ()
() 0 0)
243 'bc_type_list' '' 'bc_type_list' 239 ((VARIABLE OUT UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0 DIMENSION DUMMY) (INTEGER 4 0 0 INTEGER ()) 0 0 () (
2 0 ASSUMED_SHAPE (CONSTANT (INTEGER 4 0 0 INTEGER ()) 0 '1') () (
CONSTANT (INTEGER 4 0 0 INTEGER ()) 0 '1') ()) 0 () () () 0 0)
244 'bc_number_list' '' 'bc_number_list' 239 ((VARIABLE OUT UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0 DIMENSION OPTIONAL DUMMY) (INTEGER 4 0 0 INTEGER ())
0 0 () (2 0 ASSUMED_SHAPE (CONSTANT (INTEGER 4 0 0 INTEGER ()) 0 '1') ()
(CONSTANT (INTEGER 4 0 0 INTEGER ()) 0 '1') ()) 0 () () () 0 0)
246 'field' '' 'field' 245 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0
0 DUMMY) (DERIVED 188 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
73 'mesh' '' 'mesh' 72 ((VARIABLE INOUT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DUMMY) (DERIVED 71 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
248 'bc_type_node_list' '' 'bc_type_node_list' 245 ((VARIABLE OUT
UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 DIMENSION DUMMY) (INTEGER 4 0 0 INTEGER
()) 0 0 () (1 0 ASSUMED_SHAPE (CONSTANT (INTEGER 4 0 0 INTEGER ()) 0 '1')
()) 0 () () () 0 0)
251 'types' '' 'types' 249 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0
0 DIMENSION DUMMY) (CHARACTER 1 0 0 CHARACTER (())) 0 0 () (1 0
ASSUMED_SHAPE (CONSTANT (INTEGER 4 0 0 INTEGER ()) 0 '1') ()) 0 () () ()
0 0)
252 'bc_type_node_list' '' 'bc_type_node_list' 249 ((VARIABLE OUT
UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 DIMENSION DUMMY) (INTEGER 4 0 0 INTEGER
()) 0 0 () (2 0 ASSUMED_SHAPE (CONSTANT (INTEGER 4 0 0 INTEGER ()) 0 '1')
() (CONSTANT (INTEGER 4 0 0 INTEGER ()) 0 '1') ()) 0 () () () 0 0)
250 'field' '' 'field' 249 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0
0 DUMMY) (DERIVED 194 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
255 'node' '' 'node' 253 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DUMMY) (INTEGER 4 0 0 INTEGER ()) 0 0 () () 0 () () () 0 0)
256 'rhs' '' 'rhs' 253 ((VARIABLE INOUT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
OPTIONAL DUMMY) (DERIVED 188 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
254 'matrix' '' 'matrix' 253 ((VARIABLE INOUT UNKNOWN-PROC UNKNOWN
UNKNOWN 0 0 DUMMY) (DERIVED 324 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
257 'reference_value' '' 'reference_value' 253 ((VARIABLE IN
UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 OPTIONAL DUMMY) (REAL 8 0 0 REAL ()) 0
0 () () 0 () () () 0 0)
258 'reference_node_owned' '' 'reference_node_owned' 253 ((VARIABLE IN
UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 OPTIONAL DUMMY) (LOGICAL 4 0 0 LOGICAL
()) 0 0 () () 0 () () () 0 0)
260 'matrix' '' 'matrix' 259 ((VARIABLE INOUT UNKNOWN-PROC UNKNOWN
UNKNOWN 0 0 DUMMY) (DERIVED 446 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
261 'node' '' 'node' 259 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DUMMY) (INTEGER 4 0 0 INTEGER ()) 0 0 () () 0 () () () 0 0)
262 'rhs' '' 'rhs' 259 ((VARIABLE INOUT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
OPTIONAL DUMMY) (DERIVED 194 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
263 'reference_value' '' 'reference_value' 259 ((VARIABLE IN
UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 DIMENSION OPTIONAL DUMMY) (REAL 8 0 0
REAL ()) 0 0 () (1 0 EXPLICIT (CONSTANT (INTEGER 4 0 0 INTEGER ()) 0 '1')
(FUNCTION (INTEGER 4 0 0 INTEGER ()) 0 447 (('' (VARIABLE (DERIVED 446 0
0 DERIVED ()) 0 260 ())) ('' (CONSTANT (INTEGER 4 0 0 INTEGER ()) 0 '1')))
'petsc_csr_blocks_withdim' 1 448)) 0 () () () 0 0)
264 'reference_node_owned' '' 'reference_node_owned' 259 ((VARIABLE IN
UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 OPTIONAL DUMMY) (LOGICAL 4 0 0 LOGICAL
()) 0 0 () () 0 () () () 0 0)
266 'field' '' 'field' 265 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0
0 DUMMY) (DERIVED 188 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
267 'type' '' 'type' 265 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DUMMY) (CHARACTER 1 0 0 CHARACTER (())) 0 0 () () 0 () () () 0 0)
269 'field' '' 'field' 268 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0
0 DUMMY) (DERIVED 194 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
270 'type' '' 'type' 268 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DUMMY) (CHARACTER 1 0 0 CHARACTER (())) 0 0 () () 0 () () () 0 0)
272 'field' '' 'field' 271 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0
0 DUMMY) (DERIVED 188 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
273 'name' '' 'name' 271 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DUMMY) (CHARACTER 1 0 0 CHARACTER (())) 0 0 () () 0 () () () 0 0)
275 'field' '' 'field' 274 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0
0 DUMMY) (DERIVED 194 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
276 'name' '' 'name' 274 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DUMMY) (CHARACTER 1 0 0 CHARACTER (())) 0 0 () () 0 () () () 0 0)
84 'states' '' 'states' 83 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0
0 DIMENSION DUMMY) (DERIVED 449 0 0 DERIVED ()) 0 0 () (1 0
ASSUMED_SHAPE (CONSTANT (INTEGER 4 0 0 INTEGER ()) 0 '1') ()) 0 () () ()
0 0)
278 'field' '' 'field' 277 ((VARIABLE INOUT UNKNOWN-PROC UNKNOWN UNKNOWN
0 0 DUMMY) (DERIVED 188 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
280 'field' '' 'field' 279 ((VARIABLE INOUT UNKNOWN-PROC UNKNOWN UNKNOWN
0 0 DUMMY) (DERIVED 194 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
282 'matrix' '' 'matrix' 281 ((VARIABLE INOUT UNKNOWN-PROC UNKNOWN
UNKNOWN 0 0 DUMMY) (DERIVED 324 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
283 'rhs' '' 'rhs' 281 ((VARIABLE INOUT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DUMMY) (DERIVED 188 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
284 'field' '' 'field' 281 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0
0 DUMMY) (DERIVED 188 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
285 'dt' '' 'dt' 281 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
OPTIONAL DUMMY) (REAL 8 0 0 REAL ()) 0 0 () () 0 () () () 0 0)
287 'lhs' '' 'lhs' 286 ((VARIABLE INOUT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DUMMY) (DERIVED 188 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
288 'rhs' '' 'rhs' 286 ((VARIABLE INOUT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DUMMY) (DERIVED 188 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
289 'field' '' 'field' 286 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0
0 DUMMY) (DERIVED 188 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
290 'dt' '' 'dt' 286 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DUMMY) (REAL 8 0 0 REAL ()) 0 0 () () 0 () () () 0 0)
292 'matrix' '' 'matrix' 291 ((VARIABLE INOUT UNKNOWN-PROC UNKNOWN
UNKNOWN 0 0 DUMMY) (DERIVED 450 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
293 'rhs' '' 'rhs' 291 ((VARIABLE INOUT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
OPTIONAL DUMMY) (DERIVED 194 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
294 'field' '' 'field' 291 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0
0 DUMMY) (DERIVED 194 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
295 'dt' '' 'dt' 291 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
OPTIONAL DUMMY) (REAL 8 0 0 REAL ()) 0 0 () () 0 () () () 0 0)
297 'matrix' '' 'matrix' 296 ((VARIABLE INOUT UNKNOWN-PROC UNKNOWN
UNKNOWN 0 0 DUMMY) (DERIVED 446 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
298 'rhs' '' 'rhs' 296 ((VARIABLE INOUT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
OPTIONAL DUMMY) (DERIVED 194 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
299 'field' '' 'field' 296 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0
0 DUMMY) (DERIVED 194 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
300 'dt' '' 'dt' 296 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
OPTIONAL DUMMY) (REAL 8 0 0 REAL ()) 0 0 () () 0 () () () 0 0)
302 'matrix' '' 'matrix' 301 ((VARIABLE INOUT UNKNOWN-PROC UNKNOWN
UNKNOWN 0 0 DUMMY) (DERIVED 324 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
303 'rhs' '' 'rhs' 301 ((VARIABLE INOUT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
OPTIONAL DUMMY) (DERIVED 194 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
304 'field' '' 'field' 301 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0
0 DUMMY) (DERIVED 194 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
305 'dt' '' 'dt' 301 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
OPTIONAL DUMMY) (REAL 8 0 0 REAL ()) 0 0 () () 0 () () () 0 0)
306 'dim' '' 'dim' 301 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DUMMY) (INTEGER 4 0 0 INTEGER ()) 0 0 () () 0 () () () 0 0)
310 'field' '' 'field' 307 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0
0 DUMMY) (DERIVED 194 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
309 'rhs' '' 'rhs' 307 ((VARIABLE INOUT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
OPTIONAL DUMMY) (DERIVED 194 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
308 'lhs' '' 'lhs' 307 ((VARIABLE INOUT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DUMMY) (DERIVED 194 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
311 'dt' '' 'dt' 307 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
OPTIONAL DUMMY) (REAL 8 0 0 REAL ()) 0 0 () () 0 () () () 0 0)
313 'lhs' '' 'lhs' 312 ((VARIABLE INOUT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DUMMY) (DERIVED 188 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
314 'rhs' '' 'rhs' 312 ((VARIABLE INOUT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
OPTIONAL DUMMY) (DERIVED 194 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
315 'field' '' 'field' 312 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0
0 DUMMY) (DERIVED 194 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
316 'dt' '' 'dt' 312 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
OPTIONAL DUMMY) (REAL 8 0 0 REAL ()) 0 0 () () 0 () () () 0 0)
317 'dim' '' 'dim' 312 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DUMMY) (INTEGER 4 0 0 INTEGER ()) 0 0 () () 0 () () () 0 0)
67 'input_states' '' 'input_states' 66 ((VARIABLE IN UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0 DIMENSION DUMMY) (DERIVED 449 0 0 DERIVED ()) 0 0 ()
(1 0 ASSUMED_SHAPE (CONSTANT (INTEGER 4 0 0 INTEGER ()) 0 '1') ()) 0 ()
() () 0 0)
68 'collapsed_states' '' 'collapsed_states' 66 ((VARIABLE INOUT
UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 DIMENSION DUMMY) (DERIVED 449 0 0
DERIVED ()) 0 0 () (1 0 EXPLICIT (CONSTANT (INTEGER 4 0 0 INTEGER ()) 0
'1') (FUNCTION (INTEGER 4 0 0 INTEGER ()) 0 451 (('' (VARIABLE (DERIVED
449 0 0 DERIVED ()) 1 67 ((ARRAY (FULL 0))))) ('' ()) ('' ())) '' 0 'size'))
0 () () () 0 0)
69 'bctype' '' 'bctype' 66 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0
0 OPTIONAL DUMMY) (CHARACTER 1 0 0 CHARACTER (())) 0 0 () () 0 () () ()
0 0)
446 'petsc_csr_matrix' 'sparse_tools_petsc' 'petsc_csr_matrix' 1 ((
DERIVED UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 POINTER_COMP) (
UNKNOWN 0 0 0 UNKNOWN ()) 0 0 () () 0 ((452 'm' (INTEGER 8 0 0 INTEGER ())
() (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0)
UNKNOWN-ACCESS ()) (453 'row_numbering' (DERIVED 454 0 0 DERIVED ()) ()
(UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0)
UNKNOWN-ACCESS (STRUCTURE (DERIVED 454 0 0 DERIVED ()) 0 (((NULL (
UNKNOWN 0 0 0 UNKNOWN ()) 0) ()) (() ()) (() ()) (() ()) (() ()) ((NULL
(UNKNOWN 0 0 0 UNKNOWN ()) 0) ()) ((NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0) ())
((NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0) ()) ((CONSTANT (CHARACTER 1 0 0
CHARACTER (())) 0 101
'                                                                                                     ')
())) ())) (455 'column_numbering' (DERIVED 454 0 0 DERIVED ()) () (
UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0)
UNKNOWN-ACCESS (STRUCTURE (DERIVED 454 0 0 DERIVED ()) 0 (((NULL (
UNKNOWN 0 0 0 UNKNOWN ()) 0) ()) (() ()) (() ()) (() ()) (() ()) ((NULL
(UNKNOWN 0 0 0 UNKNOWN ()) 0) ()) ((NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0) ())
((NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0) ()) ((CONSTANT (CHARACTER 1 0 0
CHARACTER (())) 0 101
'                                                                                                     ')
())) ())) (456 'row_halo' (DERIVED 107 0 0 DERIVED ()) () (UNKNOWN-FL
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 POINTER) UNKNOWN-ACCESS
(NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0)) (457 'column_halo' (DERIVED 107 0 0
DERIVED ()) () (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0
0 POINTER) UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0)) (458
'refcount' (DERIVED 96 0 0 DERIVED ()) () (UNKNOWN-FL UNKNOWN-INTENT
UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 POINTER) UNKNOWN-ACCESS (NULL (UNKNOWN
0 0 0 UNKNOWN ()) 0)) (459 'name' (CHARACTER 1 0 0 CHARACTER ((CONSTANT
(INTEGER 4 0 0 INTEGER ()) 0 '101'))) () (UNKNOWN-FL UNKNOWN-INTENT
UNKNOWN-PROC UNKNOWN UNKNOWN 0 0) UNKNOWN-ACCESS (CONSTANT (CHARACTER 1
0 0 CHARACTER (())) 0 101
'                                                                                                     '))
(460 'is_assembled' (LOGICAL 4 0 0 LOGICAL ()) () (UNKNOWN-FL
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0) UNKNOWN-ACCESS (
CONSTANT (LOGICAL 4 0 0 LOGICAL ()) 0 0))) PUBLIC (() () () ()) () 0 0
75554753)
340 'scalar_boundary_conditions_ptr' 'fields_data_types'
'scalar_boundary_conditions_ptr' 1 ((DERIVED UNKNOWN-INTENT UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0 POINTER_COMP) (UNKNOWN 0 0 0 UNKNOWN ()) 0 0 () () 0
((461 'boundary_condition' (DERIVED 462 0 0 DERIVED ()) (1 0 DEFERRED ()
()) (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DIMENSION POINTER) UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0)))
PUBLIC (() () () ()) () 0 0 47502942)
356 'vector_boundary_conditions_ptr' 'fields_data_types'
'vector_boundary_conditions_ptr' 1 ((DERIVED UNKNOWN-INTENT UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0 POINTER_COMP) (UNKNOWN 0 0 0 UNKNOWN ()) 0 0 () () 0
((463 'boundary_condition' (DERIVED 464 0 0 DERIVED ()) (1 0 DEFERRED ()
()) (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DIMENSION POINTER) UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0)))
PUBLIC (() () () ()) () 0 0 84476181)
96 'refcount_type' 'reference_counting' 'refcount_type' 1 ((DERIVED
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 POINTER_COMP) (UNKNOWN 0
0 0 UNKNOWN ()) 0 0 () () 0 ((465 'prev' (DERIVED 96 0 0 DERIVED ()) ()
(UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 POINTER)
UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0)) (466 'next' (
DERIVED 96 0 0 DERIVED ()) () (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0 POINTER) UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN
()) 0)) (467 'count' (INTEGER 4 0 0 INTEGER ()) () (UNKNOWN-FL
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0) UNKNOWN-ACCESS (
CONSTANT (INTEGER 4 0 0 INTEGER ()) 0 '0')) (468 'id' (INTEGER 4 0 0
INTEGER ()) () (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0
0) UNKNOWN-ACCESS ()) (469 'name' (CHARACTER 1 0 0 CHARACTER ((CONSTANT
(INTEGER 4 0 0 INTEGER ()) 0 '101'))) () (UNKNOWN-FL UNKNOWN-INTENT
UNKNOWN-PROC UNKNOWN UNKNOWN 0 0) UNKNOWN-ACCESS ()) (470 'type' (
CHARACTER 1 0 0 CHARACTER ((CONSTANT (INTEGER 4 0 0 INTEGER ()) 0 '101')))
() (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0)
UNKNOWN-ACCESS ()) (471 'tagged' (LOGICAL 4 0 0 LOGICAL ()) () (
UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0)
UNKNOWN-ACCESS (CONSTANT (LOGICAL 4 0 0 LOGICAL ()) 0 0))) PUBLIC (() ()
() ()) () 0 0 25948645)
138 'integer_vector' 'futils' 'integer_vector' 1 ((DERIVED
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 POINTER_COMP) (UNKNOWN 0
0 0 UNKNOWN ()) 0 0 () () 0 ((472 'ptr' (INTEGER 4 0 0 INTEGER ()) (1 0
DEFERRED () ()) (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN
0 0 DIMENSION POINTER) UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0)))
PUBLIC (() () () ()) () 0 0 9661976)
454 'petsc_numbering_type' 'petsc_tools' 'petsc_numbering_type' 1 ((
DERIVED UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 POINTER_COMP) (
UNKNOWN 0 0 0 UNKNOWN ()) 0 0 () () 0 ((473 'halo' (DERIVED 107 0 0
DERIVED ()) () (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0
0 POINTER) UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0)) (474
'nprivatenodes' (INTEGER 4 0 0 INTEGER ()) () (UNKNOWN-FL UNKNOWN-INTENT
UNKNOWN-PROC UNKNOWN UNKNOWN 0 0) UNKNOWN-ACCESS ()) (475
'universal_length' (INTEGER 4 0 0 INTEGER ()) () (UNKNOWN-FL
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0) UNKNOWN-ACCESS ()) (
476 'offset' (INTEGER 4 0 0 INTEGER ()) () (UNKNOWN-FL UNKNOWN-INTENT
UNKNOWN-PROC UNKNOWN UNKNOWN 0 0) UNKNOWN-ACCESS ()) (477 'gnn2unn' (
INTEGER 4 0 0 INTEGER ()) (2 0 DEFERRED () () () ()) (UNKNOWN-FL
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 DIMENSION POINTER)
UNKNOWN-ACCESS ()) (478 'ghost_nodes' (INTEGER 4 0 0 INTEGER ()) (1 0
DEFERRED () ()) (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN
0 0 DIMENSION POINTER) UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0))
(479 'ghost2unn' (INTEGER 4 0 0 INTEGER ()) (2 0 DEFERRED () () () ()) (
UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 DIMENSION
POINTER) UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0)) (480
'refcount' (DERIVED 96 0 0 DERIVED ()) () (UNKNOWN-FL UNKNOWN-INTENT
UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 POINTER) UNKNOWN-ACCESS (NULL (UNKNOWN
0 0 0 UNKNOWN ()) 0)) (481 'name' (CHARACTER 1 0 0 CHARACTER ((CONSTANT
(INTEGER 4 0 0 INTEGER ()) 0 '101'))) () (UNKNOWN-FL UNKNOWN-INTENT
UNKNOWN-PROC UNKNOWN UNKNOWN 0 0) UNKNOWN-ACCESS (CONSTANT (CHARACTER 1
0 0 CHARACTER (())) 0 101
'                                                                                                     ')))
PUBLIC (() () () ()) () 0 0 16529124)
449 'state_type' 'state_module' 'state_type' 1 ((DERIVED UNKNOWN-INTENT
UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 POINTER_COMP) (UNKNOWN 0 0 0 UNKNOWN ())
0 0 () () 0 ((482 'name' (CHARACTER 1 0 0 CHARACTER ((CONSTANT (INTEGER
4 0 0 INTEGER ()) 0 '101'))) () (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0) UNKNOWN-ACCESS (CONSTANT (CHARACTER 1 0 0 CHARACTER
(())) 0 101
'                                                                                                     '))
(483 'option_path' (CHARACTER 1 0 0 CHARACTER ((CONSTANT (INTEGER 4 0 0
INTEGER ()) 0 '8192'))) () (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0) UNKNOWN-ACCESS (CONSTANT (CHARACTER 1 0 0 CHARACTER
(())) 0 8192
'/uninitialised_path/                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            '))
(484 'vector_names' (CHARACTER 1 0 0 CHARACTER ((CONSTANT (INTEGER 4 0 0
INTEGER ()) 0 '101'))) (1 0 DEFERRED () ()) (UNKNOWN-FL UNKNOWN-INTENT
UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 DIMENSION POINTER) UNKNOWN-ACCESS (
NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0)) (485 'scalar_names' (CHARACTER 1 0 0
CHARACTER ((CONSTANT (INTEGER 4 0 0 INTEGER ()) 0 '101'))) (1 0 DEFERRED
() ()) (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DIMENSION POINTER) UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0)) (
486 'mesh_names' (CHARACTER 1 0 0 CHARACTER ((CONSTANT (INTEGER 4 0 0
INTEGER ()) 0 '101'))) (1 0 DEFERRED () ()) (UNKNOWN-FL UNKNOWN-INTENT
UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 DIMENSION POINTER) UNKNOWN-ACCESS (
NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0)) (487 'halo_names' (CHARACTER 1 0 0
CHARACTER ((CONSTANT (INTEGER 4 0 0 INTEGER ()) 0 '101'))) (1 0 DEFERRED
() ()) (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DIMENSION POINTER) UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0)) (
488 'tensor_names' (CHARACTER 1 0 0 CHARACTER ((CONSTANT (INTEGER 4 0 0
INTEGER ()) 0 '101'))) (1 0 DEFERRED () ()) (UNKNOWN-FL UNKNOWN-INTENT
UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 DIMENSION POINTER) UNKNOWN-ACCESS (
NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0)) (489 'csr_sparsity_names' (
CHARACTER 1 0 0 CHARACTER ((CONSTANT (INTEGER 4 0 0 INTEGER ()) 0 '101')))
(1 0 DEFERRED () ()) (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN
UNKNOWN 0 0 DIMENSION POINTER) UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0
UNKNOWN ()) 0)) (490 'csr_matrix_names' (CHARACTER 1 0 0 CHARACTER ((
CONSTANT (INTEGER 4 0 0 INTEGER ()) 0 '101'))) (1 0 DEFERRED () ()) (
UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 DIMENSION
POINTER) UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0)) (491
'block_csr_matrix_names' (CHARACTER 1 0 0 CHARACTER ((CONSTANT (INTEGER
4 0 0 INTEGER ()) 0 '101'))) (1 0 DEFERRED () ()) (UNKNOWN-FL
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 DIMENSION POINTER)
UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0)) (492
'petsc_csr_matrix_names' (CHARACTER 1 0 0 CHARACTER ((CONSTANT (INTEGER
4 0 0 INTEGER ()) 0 '101'))) (1 0 DEFERRED () ()) (UNKNOWN-FL
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 DIMENSION POINTER)
UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0)) (493 'vector_fields'
(DERIVED 494 0 0 DERIVED ()) (1 0 DEFERRED () ()) (UNKNOWN-FL
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 DIMENSION POINTER)
UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0)) (495 'tensor_fields'
(DERIVED 496 0 0 DERIVED ()) (1 0 DEFERRED () ()) (UNKNOWN-FL
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 DIMENSION POINTER)
UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0)) (497 'scalar_fields'
(DERIVED 498 0 0 DERIVED ()) (1 0 DEFERRED () ()) (UNKNOWN-FL
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 DIMENSION POINTER)
UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0)) (499 'meshes' (
DERIVED 500 0 0 DERIVED ()) (1 0 DEFERRED () ()) (UNKNOWN-FL
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 DIMENSION POINTER)
UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0)) (501 'halos' (
DERIVED 502 0 0 DERIVED ()) (1 0 DEFERRED () ()) (UNKNOWN-FL
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 DIMENSION POINTER)
UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0)) (503 'csr_sparsities'
(DERIVED 504 0 0 DERIVED ()) (1 0 DEFERRED () ()) (UNKNOWN-FL
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 DIMENSION POINTER)
UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0)) (505 'csr_matrices'
(DERIVED 506 0 0 DERIVED ()) (1 0 DEFERRED () ()) (UNKNOWN-FL
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 DIMENSION POINTER)
UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0)) (507
'block_csr_matrices' (DERIVED 508 0 0 DERIVED ()) (1 0 DEFERRED () ()) (
UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 DIMENSION
POINTER) UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0)) (509
'petsc_csr_matrices' (DERIVED 510 0 0 DERIVED ()) (1 0 DEFERRED () ()) (
UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 DIMENSION
POINTER) UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0))) PUBLIC (()
() () ()) () 0 0 7575341)
450 'block_csr_matrix' 'sparse_tools' 'block_csr_matrix' 1 ((DERIVED
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 POINTER_COMP) (UNKNOWN 0
0 0 UNKNOWN ()) 0 0 () () 0 ((511 'sparsity' (DERIVED 319 0 0 DERIVED ())
() (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0)
UNKNOWN-ACCESS (STRUCTURE (DERIVED 319 0 0 DERIVED ()) 0 (((NULL (
UNKNOWN 0 0 0 UNKNOWN ()) 0) ()) ((NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0) ())
((NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0) ()) (() ()) ((NULL (UNKNOWN 0 0 0
UNKNOWN ()) 0) ()) ((NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0) ()) ((NULL (
UNKNOWN 0 0 0 UNKNOWN ()) 0) ()) ((CONSTANT (CHARACTER 1 0 0 CHARACTER (
())) 0 101
'                                                                                                     ')
()) ((CONSTANT (LOGICAL 4 0 0 LOGICAL ()) 0 0) ()) ((CONSTANT (LOGICAL 4
0 0 LOGICAL ()) 0 0) ())) ())) (512 'val' (DERIVED 513 0 0 DERIVED ()) (
2 0 DEFERRED () () () ()) (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0 DIMENSION POINTER) UNKNOWN-ACCESS (NULL (UNKNOWN 0 0
0 UNKNOWN ()) 0)) (514 'contiguous_val' (REAL 8 0 0 REAL ()) (1 0
DEFERRED () ()) (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN
0 0 DIMENSION POINTER) UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0))
(515 'ival' (DERIVED 138 0 0 DERIVED ()) (2 0 DEFERRED () () () ()) (
UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 DIMENSION
POINTER) UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0)) (516 'blocks'
(INTEGER 4 0 0 INTEGER ()) (1 0 EXPLICIT (CONSTANT (INTEGER 4 0 0
INTEGER ()) 0 '1') (CONSTANT (INTEGER 4 0 0 INTEGER ()) 0 '2')) (
UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 DIMENSION)
UNKNOWN-ACCESS (ARRAY (INTEGER 4 0 0 INTEGER ()) 1 (((CONSTANT (INTEGER
4 0 0 INTEGER ()) 0 '0') ()) ((CONSTANT (INTEGER 4 0 0 INTEGER ()) 0 '0')
())) ('2'))) (517 'clone' (LOGICAL 4 0 0 LOGICAL ()) () (UNKNOWN-FL
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0) UNKNOWN-ACCESS (
CONSTANT (LOGICAL 4 0 0 LOGICAL ()) 0 0)) (518 'external_val' (LOGICAL 4
0 0 LOGICAL ()) () (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN
UNKNOWN 0 0) UNKNOWN-ACCESS (CONSTANT (LOGICAL 4 0 0 LOGICAL ()) 0 0)) (
519 'columns' (INTEGER 4 0 0 INTEGER ()) () (UNKNOWN-FL UNKNOWN-INTENT
UNKNOWN-PROC UNKNOWN UNKNOWN 0 0) UNKNOWN-ACCESS ()) (520 'refcount' (
DERIVED 96 0 0 DERIVED ()) () (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0 POINTER) UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN
()) 0)) (521 'name' (CHARACTER 1 0 0 CHARACTER ((CONSTANT (INTEGER 4 0 0
INTEGER ()) 0 '101'))) () (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0) UNKNOWN-ACCESS (CONSTANT (CHARACTER 1 0 0 CHARACTER
(())) 0 101
'                                                                                                     '))
(522 'diagonal' (LOGICAL 4 0 0 LOGICAL ()) () (UNKNOWN-FL UNKNOWN-INTENT
UNKNOWN-PROC UNKNOWN UNKNOWN 0 0) UNKNOWN-ACCESS (CONSTANT (LOGICAL 4 0
0 LOGICAL ()) 0 0)) (523 'equal_diagonal_blocks' (LOGICAL 4 0 0 LOGICAL
()) () (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0)
UNKNOWN-ACCESS (CONSTANT (LOGICAL 4 0 0 LOGICAL ()) 0 0)) (524 'ksp' (
INTEGER 8 0 0 INTEGER ()) () (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0 POINTER) UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN
()) 0))) PUBLIC (() () () ()) () 0 0 61921235)
508 'block_csr_matrix_pointer' 'sparse_tools' 'block_csr_matrix_pointer'
1 ((DERIVED UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 POINTER_COMP)
(UNKNOWN 0 0 0 UNKNOWN ()) 0 0 () () 0 ((525 'ptr' (DERIVED 450 0 0
DERIVED ()) () (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0
0 POINTER) UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0))) PUBLIC (
() () () ()) () 0 0 83606449)
506 'csr_matrix_pointer' 'sparse_tools' 'csr_matrix_pointer' 1 ((
DERIVED UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 POINTER_COMP) (
UNKNOWN 0 0 0 UNKNOWN ()) 0 0 () () 0 ((526 'ptr' (DERIVED 324 0 0
DERIVED ()) () (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0
0 POINTER) UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0))) PUBLIC (
() () () ()) () 0 0 82770239)
504 'csr_sparsity_pointer' 'sparse_tools' 'csr_sparsity_pointer' 1 ((
DERIVED UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 POINTER_COMP) (
UNKNOWN 0 0 0 UNKNOWN ()) 0 0 () () 0 ((527 'ptr' (DERIVED 319 0 0
DERIVED ()) () (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0
0 POINTER) UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0))) PUBLIC (
() () () ()) () 0 0 9687815)
513 'real_vector' 'futils' 'real_vector' 1 ((DERIVED UNKNOWN-INTENT
UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 POINTER_COMP) (UNKNOWN 0 0 0 UNKNOWN ())
0 0 () () 0 ((528 'ptr' (REAL 8 0 0 REAL ()) (1 0 DEFERRED () ()) (
UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 DIMENSION
POINTER) UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0))) PUBLIC (()
() () ()) () 0 0 72870256)
447 'blocks' '' 'blocks' 259 ((PROCEDURE UNKNOWN-INTENT UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0 FUNCTION) (UNKNOWN 0 0 0 UNKNOWN ()) 0 0 () () 447 ()
() () 0 0)
451 'size' '(intrinsic)' 'size' 66 ((PROCEDURE UNKNOWN-INTENT
UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 FUNCTION) (UNKNOWN 0 0 0 UNKNOWN ()) 0
0 () () 451 () () () 0 0)
510 'petsc_csr_matrix_pointer' 'sparse_tools_petsc'
'petsc_csr_matrix_pointer' 1 ((DERIVED UNKNOWN-INTENT UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0 POINTER_COMP) (UNKNOWN 0 0 0 UNKNOWN ()) 0 0 () () 0
((529 'ptr' (DERIVED 446 0 0 DERIVED ()) () (UNKNOWN-FL UNKNOWN-INTENT
UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 POINTER) UNKNOWN-ACCESS (NULL (UNKNOWN
0 0 0 UNKNOWN ()) 0))) PUBLIC (() () () ()) () 0 0 85445535)
500 'mesh_pointer' 'fields_data_types' 'mesh_pointer' 1 ((DERIVED
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 POINTER_COMP) (UNKNOWN 0
0 0 UNKNOWN ()) 0 0 () () 0 ((530 'ptr' (DERIVED 71 0 0 DERIVED ()) () (
UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 POINTER)
UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0))) PUBLIC (() () () ())
() 0 0 65432384)
462 'scalar_boundary_condition' 'fields_data_types'
'scalar_boundary_condition' 1 ((DERIVED UNKNOWN-INTENT UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0 POINTER_COMP) (UNKNOWN 0 0 0 UNKNOWN ()) 0 0 () () 0
((531 'name' (CHARACTER 1 0 0 CHARACTER ((CONSTANT (INTEGER 4 0 0
INTEGER ()) 0 '101'))) () (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0) UNKNOWN-ACCESS ()) (532 'type' (CHARACTER 1 0 0
CHARACTER ((CONSTANT (INTEGER 4 0 0 INTEGER ()) 0 '101'))) () (
UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0)
UNKNOWN-ACCESS (CONSTANT (CHARACTER 1 0 0 CHARACTER (())) 0 101
'                                                                                                     '))
(533 'surface_element_list' (INTEGER 4 0 0 INTEGER ()) (1 0 DEFERRED ()
()) (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DIMENSION POINTER) UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0)) (
534 'surface_node_list' (INTEGER 4 0 0 INTEGER ()) (1 0 DEFERRED () ())
(UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 DIMENSION
POINTER) UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0)) (535
'surface_mesh' (DERIVED 71 0 0 DERIVED ()) () (UNKNOWN-FL UNKNOWN-INTENT
UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 POINTER) UNKNOWN-ACCESS ()) (536
'surface_fields' (DERIVED 188 0 0 DERIVED ()) (1 0 DEFERRED () ()) (
UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 DIMENSION
POINTER) UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0)) (537
'option_path' (CHARACTER 1 0 0 CHARACTER ((CONSTANT (INTEGER 4 0 0
INTEGER ()) 0 '8192'))) () (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0) UNKNOWN-ACCESS (CONSTANT (CHARACTER 1 0 0 CHARACTER
(())) 0 8192
'/uninitialised_path/                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            ')))
PUBLIC (() () () ()) () 0 0 50740068)
498 'scalar_field_pointer' 'fields_data_types' 'scalar_field_pointer' 1
((DERIVED UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 POINTER_COMP)
(UNKNOWN 0 0 0 UNKNOWN ()) 0 0 () () 0 ((538 'ptr' (DERIVED 188 0 0
DERIVED ()) () (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0
0 POINTER) UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0))) PUBLIC (
() () () ()) () 0 0 54828570)
496 'tensor_field_pointer' 'fields_data_types' 'tensor_field_pointer' 1
((DERIVED UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 POINTER_COMP)
(UNKNOWN 0 0 0 UNKNOWN ()) 0 0 () () 0 ((539 'ptr' (DERIVED 540 0 0
DERIVED ()) () (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0
0 POINTER) UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0))) PUBLIC (
() () () ()) () 0 0 16722823)
464 'vector_boundary_condition' 'fields_data_types'
'vector_boundary_condition' 1 ((DERIVED UNKNOWN-INTENT UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0 POINTER_COMP) (UNKNOWN 0 0 0 UNKNOWN ()) 0 0 () () 0
((541 'name' (CHARACTER 1 0 0 CHARACTER ((CONSTANT (INTEGER 4 0 0
INTEGER ()) 0 '101'))) () (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0) UNKNOWN-ACCESS ()) (542 'type' (CHARACTER 1 0 0
CHARACTER ((CONSTANT (INTEGER 4 0 0 INTEGER ()) 0 '101'))) () (
UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0)
UNKNOWN-ACCESS (CONSTANT (CHARACTER 1 0 0 CHARACTER (())) 0 101
'                                                                                                     '))
(543 'applies' (LOGICAL 4 0 0 LOGICAL ()) (1 0 EXPLICIT (CONSTANT (
INTEGER 4 0 0 INTEGER ()) 0 '1') (CONSTANT (INTEGER 4 0 0 INTEGER ()) 0
'3')) (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DIMENSION) UNKNOWN-ACCESS ()) (544 'surface_element_list' (INTEGER 4 0 0
INTEGER ()) (1 0 DEFERRED () ()) (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0 DIMENSION POINTER) UNKNOWN-ACCESS (NULL (UNKNOWN 0 0
0 UNKNOWN ()) 0)) (545 'surface_node_list' (INTEGER 4 0 0 INTEGER ()) (
1 0 DEFERRED () ()) (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN
UNKNOWN 0 0 DIMENSION POINTER) UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0
UNKNOWN ()) 0)) (546 'surface_mesh' (DERIVED 71 0 0 DERIVED ()) () (
UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 POINTER)
UNKNOWN-ACCESS ()) (547 'surface_fields' (DERIVED 194 0 0 DERIVED ()) (
1 0 DEFERRED () ()) (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN
UNKNOWN 0 0 DIMENSION POINTER) UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0
UNKNOWN ()) 0)) (548 'scalar_surface_fields' (DERIVED 188 0 0 DERIVED ())
(1 0 DEFERRED () ()) (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN
UNKNOWN 0 0 DIMENSION POINTER) UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0
UNKNOWN ()) 0)) (549 'option_path' (CHARACTER 1 0 0 CHARACTER ((
CONSTANT (INTEGER 4 0 0 INTEGER ()) 0 '8192'))) () (UNKNOWN-FL
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0) UNKNOWN-ACCESS (
CONSTANT (CHARACTER 1 0 0 CHARACTER (())) 0 8192
'/uninitialised_path/                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            ')))
PUBLIC (() () () ()) () 0 0 42205165)
494 'vector_field_pointer' 'fields_data_types' 'vector_field_pointer' 1
((DERIVED UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 POINTER_COMP)
(UNKNOWN 0 0 0 UNKNOWN ()) 0 0 () () 0 ((550 'ptr' (DERIVED 194 0 0
DERIVED ()) () (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0
0 POINTER) UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0))) PUBLIC (
() () () ()) () 0 0 74448721)
540 'tensor_field' 'fields_data_types' 'tensor_field' 1 ((DERIVED
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 POINTER_COMP) (UNKNOWN 0
0 0 UNKNOWN ()) 0 0 () () 0 ((551 'val' (REAL 8 0 0 REAL ()) (3 0
DEFERRED () () () () () ()) (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0 DIMENSION POINTER) UNKNOWN-ACCESS ()) (552 'wrapped'
(LOGICAL 4 0 0 LOGICAL ()) () (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0) UNKNOWN-ACCESS (CONSTANT (LOGICAL 4 0 0 LOGICAL ())
0 1)) (553 'field_type' (INTEGER 4 0 0 INTEGER ()) () (UNKNOWN-FL
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0) UNKNOWN-ACCESS (
CONSTANT (INTEGER 4 0 0 INTEGER ()) 0 '0')) (554 'name' (CHARACTER 1 0 0
CHARACTER ((CONSTANT (INTEGER 4 0 0 INTEGER ()) 0 '101'))) () (
UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0)
UNKNOWN-ACCESS ()) (555 'dim' (INTEGER 4 0 0 INTEGER ()) (1 0 EXPLICIT (
CONSTANT (INTEGER 4 0 0 INTEGER ()) 0 '1') (CONSTANT (INTEGER 4 0 0
INTEGER ()) 0 '2')) (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN
UNKNOWN 0 0 DIMENSION) UNKNOWN-ACCESS ()) (556 'option_path' (CHARACTER
1 0 0 CHARACTER ((CONSTANT (INTEGER 4 0 0 INTEGER ()) 0 '8192'))) () (
UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0)
UNKNOWN-ACCESS (CONSTANT (CHARACTER 1 0 0 CHARACTER (())) 0 8192
'/uninitialised_path/                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            '))
(557 'mesh' (DERIVED 71 0 0 DERIVED ()) () (UNKNOWN-FL UNKNOWN-INTENT
UNKNOWN-PROC UNKNOWN UNKNOWN 0 0) UNKNOWN-ACCESS (STRUCTURE (DERIVED 71
0 0 DERIVED ()) 0 ((() ()) ((CONSTANT (LOGICAL 4 0 0 LOGICAL ()) 0 1) ())
((STRUCTURE (DERIVED 88 0 0 DERIVED ()) 0 ((() ()) (() ()) (() ()) (() ())
((NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0) ()) ((NULL (UNKNOWN 0 0 0 UNKNOWN ())
0) ()) ((NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0) ()) ((NULL (UNKNOWN 0 0 0
UNKNOWN ()) 0) ()) ((NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0) ()) ((NULL (
UNKNOWN 0 0 0 UNKNOWN ()) 0) ()) ((NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0) ())
((STRUCTURE (DERIVED 89 0 0 DERIVED ()) 0 ((() ()) (() ()) (() ()) (() ())
((NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0) ()) ((NULL (UNKNOWN 0 0 0 UNKNOWN ())
0) ()) (() ()) ((NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0) ()) (() ())) ()) ())
((NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0) ()) ((NULL (UNKNOWN 0 0 0 UNKNOWN ())
0) ()) ((NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0) ()) ((NULL (UNKNOWN 0 0 0
UNKNOWN ()) 0) ()) (() ())) ()) ()) (() ()) (() ()) (() ()) ((CONSTANT (
CHARACTER 1 0 0 CHARACTER (())) 0 8192
'/uninitialised_path/                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            ')
()) ((CONSTANT (INTEGER 4 0 0 INTEGER ()) 0 '0') ()) ((NULL (UNKNOWN 0 0
0 UNKNOWN ()) 0) ()) ((NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0) ()) ((NULL (
UNKNOWN 0 0 0 UNKNOWN ()) 0) ()) ((NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0) ())
((NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0) ()) ((NULL (UNKNOWN 0 0 0 UNKNOWN ())
0) ()) ((NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0) ()) ((NULL (UNKNOWN 0 0 0
UNKNOWN ()) 0) ()) ((NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0) ()) ((CONSTANT (
LOGICAL 4 0 0 LOGICAL ()) 0 0) ())) ())) (558 'refcount' (DERIVED 96 0 0
DERIVED ()) () (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0
0 POINTER) UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0)) (559
'aliased' (LOGICAL 4 0 0 LOGICAL ()) () (UNKNOWN-FL UNKNOWN-INTENT
UNKNOWN-PROC UNKNOWN UNKNOWN 0 0) UNKNOWN-ACCESS (CONSTANT (LOGICAL 4 0
0 LOGICAL ()) 0 0))) PUBLIC (() () () ()) () 0 0 25110185)
364 'picker_ptr' 'picker_data_types' 'picker_ptr' 1 ((DERIVED
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 POINTER_COMP) (UNKNOWN 0
0 0 UNKNOWN ()) 0 0 () () 0 ((560 'ptr' (DERIVED 561 0 0 DERIVED ()) ()
(UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 POINTER)
UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0))) PUBLIC (() () () ())
() 0 0 71120007)
561 'picker_type' 'picker_data_types' 'picker_type' 1 ((DERIVED
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 POINTER_COMP) (UNKNOWN 0
0 0 UNKNOWN ()) 0 0 () () 0 ((562 'name' (CHARACTER 1 0 0 CHARACTER ((
CONSTANT (INTEGER 4 0 0 INTEGER ()) 0 '101'))) () (UNKNOWN-FL
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0) UNKNOWN-ACCESS ()) (
563 'refcount' (DERIVED 96 0 0 DERIVED ()) () (UNKNOWN-FL UNKNOWN-INTENT
UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 POINTER) UNKNOWN-ACCESS (NULL (UNKNOWN
0 0 0 UNKNOWN ()) 0)) (564 'picker_id' (INTEGER 4 0 0 INTEGER ()) () (
UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0)
UNKNOWN-ACCESS (CONSTANT (INTEGER 4 0 0 INTEGER ()) 0 '0')) (565
'last_mesh_movement' (INTEGER 4 0 0 INTEGER ()) () (UNKNOWN-FL
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0) UNKNOWN-ACCESS (
CONSTANT (INTEGER 4 0 0 INTEGER ()) 0 '0'))) PUBLIC (() () () ()) () 0 0
8821665)
448 'petsc_csr_blocks_withdim' 'sparse_tools_petsc'
'petsc_csr_blocks_withdim' 1 ((PROCEDURE UNKNOWN-INTENT MODULE-PROC DECL
UNKNOWN 0 0 FUNCTION PURE ALWAYS_EXPLICIT) (INTEGER 4 0 0 INTEGER ())
566 0 (567 568) () 569 () () () 0 0)
568 'dim' '' 'dim' 566 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
OPTIONAL DUMMY) (INTEGER 4 0 0 INTEGER ()) 0 0 () () 0 () () () 0 0)
569 'blocks' '' 'blocks' 566 ((VARIABLE UNKNOWN-INTENT UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0 RESULT ALWAYS_EXPLICIT) (INTEGER 4 0 0 INTEGER ()) 0
0 () () 0 () () () 0 0)
567 'matrix' '' 'matrix' 566 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN
0 0 DUMMY) (DERIVED 446 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
502 'halo_pointer' 'halo_data_types' 'halo_pointer' 1 ((DERIVED
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 POINTER_COMP) (UNKNOWN 0
0 0 UNKNOWN ()) 0 0 () () 0 ((570 'ptr' (DERIVED 107 0 0 DERIVED ()) ()
(UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 POINTER)
UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0))) PUBLIC (() () () ())
() 0 0 72023282)
5 'add_scalar_boundary_condition' 'boundary_conditions'
'add_scalar_boundary_condition' 1 ((PROCEDURE UNKNOWN-INTENT MODULE-PROC
DECL UNKNOWN 0 0 SUBROUTINE ALWAYS_EXPLICIT) (UNKNOWN 0 0 0 UNKNOWN ())
571 0 (572 573 574 575 576) () 0 () () () 0 0)
4 'add_vector_boundary_condition' 'boundary_conditions'
'add_vector_boundary_condition' 1 ((PROCEDURE UNKNOWN-INTENT MODULE-PROC
DECL UNKNOWN 0 0 SUBROUTINE ALWAYS_EXPLICIT) (UNKNOWN 0 0 0 UNKNOWN ())
577 0 (578 579 580 581 582 583) () 0 () () () 0 0)
42 'remove_scalar_boundary_condition' 'boundary_conditions'
'remove_scalar_boundary_condition' 1 ((PROCEDURE UNKNOWN-INTENT
MODULE-PROC DECL UNKNOWN 0 0 SUBROUTINE) (UNKNOWN 0 0 0 UNKNOWN ()) 584
0 (585 586) () 0 () () () 0 0)
41 'remove_vector_boundary_condition' 'boundary_conditions'
'remove_vector_boundary_condition' 1 ((PROCEDURE UNKNOWN-INTENT
MODULE-PROC DECL UNKNOWN 0 0 SUBROUTINE) (UNKNOWN 0 0 0 UNKNOWN ()) 587
0 (588 589) () 0 () () () 0 0)
3 'add_scalar_boundary_condition_surface_elements' 'boundary_conditions'
'add_scalar_boundary_condition_surface_elements' 1 ((PROCEDURE
UNKNOWN-INTENT MODULE-PROC DECL UNKNOWN 0 0 SUBROUTINE ALWAYS_EXPLICIT)
(UNKNOWN 0 0 0 UNKNOWN ()) 590 0 (591 592 593 594 595) () 0 () () () 0 0)
572 'field' '' 'field' 571 ((VARIABLE INOUT UNKNOWN-PROC UNKNOWN UNKNOWN
0 0 DUMMY) (DERIVED 188 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
575 'boundary_ids' '' 'boundary_ids' 571 ((VARIABLE IN UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0 DIMENSION DUMMY) (INTEGER 4 0 0 INTEGER ()) 0 0 () (
1 0 ASSUMED_SHAPE (CONSTANT (INTEGER 4 0 0 INTEGER ()) 0 '1') ()) 0 () ()
() 0 0)
576 'option_path' '' 'option_path' 571 ((VARIABLE IN UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0 OPTIONAL DUMMY) (CHARACTER 1 0 0 CHARACTER (())) 0 0
() () 0 () () () 0 0)
573 'name' '' 'name' 571 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DUMMY) (CHARACTER 1 0 0 CHARACTER (())) 0 0 () () 0 () () () 0 0)
574 'type' '' 'type' 571 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DUMMY) (CHARACTER 1 0 0 CHARACTER (())) 0 0 () () 0 () () () 0 0)
591 'field' '' 'field' 590 ((VARIABLE INOUT UNKNOWN-PROC UNKNOWN UNKNOWN
0 0 DUMMY) (DERIVED 188 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
592 'name' '' 'name' 590 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DUMMY) (CHARACTER 1 0 0 CHARACTER (())) 0 0 () () 0 () () () 0 0)
593 'type' '' 'type' 590 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DUMMY) (CHARACTER 1 0 0 CHARACTER (())) 0 0 () () 0 () () () 0 0)
594 'surface_element_list' '' 'surface_element_list' 590 ((VARIABLE IN
UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 DIMENSION DUMMY) (INTEGER 4 0 0 INTEGER
()) 0 0 () (1 0 ASSUMED_SHAPE (CONSTANT (INTEGER 4 0 0 INTEGER ()) 0 '1')
()) 0 () () () 0 0)
595 'option_path' '' 'option_path' 590 ((VARIABLE IN UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0 OPTIONAL DUMMY) (CHARACTER 1 0 0 CHARACTER (())) 0 0
() () 0 () () () 0 0)
580 'type' '' 'type' 577 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DUMMY) (CHARACTER 1 0 0 CHARACTER (())) 0 0 () () 0 () () () 0 0)
579 'name' '' 'name' 577 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DUMMY) (CHARACTER 1 0 0 CHARACTER (())) 0 0 () () 0 () () () 0 0)
578 'field' '' 'field' 577 ((VARIABLE INOUT UNKNOWN-PROC UNKNOWN UNKNOWN
0 0 DUMMY) (DERIVED 194 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
581 'boundary_ids' '' 'boundary_ids' 577 ((VARIABLE IN UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0 DIMENSION DUMMY) (INTEGER 4 0 0 INTEGER ()) 0 0 () (
1 0 ASSUMED_SHAPE (CONSTANT (INTEGER 4 0 0 INTEGER ()) 0 '1') ()) 0 () ()
() 0 0)
582 'applies' '' 'applies' 577 ((VARIABLE IN UNKNOWN-PROC UNKNOWN
UNKNOWN 0 0 DIMENSION OPTIONAL DUMMY) (LOGICAL 4 0 0 LOGICAL ()) 0 0 ()
(1 0 ASSUMED_SHAPE (CONSTANT (INTEGER 4 0 0 INTEGER ()) 0 '1') ()) 0 ()
() () 0 0)
583 'option_path' '' 'option_path' 577 ((VARIABLE IN UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0 OPTIONAL DUMMY) (CHARACTER 1 0 0 CHARACTER (())) 0 0
() () 0 () () () 0 0)
586 'name' '' 'name' 584 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DUMMY) (CHARACTER 1 0 0 CHARACTER (())) 0 0 () () 0 () () () 0 0)
585 'field' '' 'field' 584 ((VARIABLE INOUT UNKNOWN-PROC UNKNOWN UNKNOWN
0 0 DUMMY) (DERIVED 188 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
588 'field' '' 'field' 587 ((VARIABLE INOUT UNKNOWN-PROC UNKNOWN UNKNOWN
0 0 DUMMY) (DERIVED 194 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
589 'name' '' 'name' 587 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DUMMY) (CHARACTER 1 0 0 CHARACTER (())) 0 0 () () 0 () () () 0 0)
)

('derive_collapsed_bcs' 0 65 'get_dg_surface_mesh' 0 70
'get_periodic_boundary_condition' 0 74 'has_scalar_surface_field' 0 78
'set_dirichlet_consistent' 0 45)
