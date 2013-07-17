GFORTRAN module version '6' created from Halos.F90 on Mon Jun 10 19:47:32 2013
MD5:a436f2314aa40077870adc9e3557298c -- If you edit this, you'll get what you deserve.

(() () () () () () () () () () () () () () () () () () () () () () ()
() () () ())

()

(('deallocate' 'halos_allocates' 2 3) ('derive_l1_from_l2_halo'
'halos_derivation' 4) ('get_universal_numbering' 'halos_numbering' 5 6)
('allocate' 'halos_allocates' 7 8) ('halo_communicator' 'halos_base' 9)
('halo_max' 'halos_communications' 10 11 12) ('halo_update'
'halos_communications' 13 14 15 16 17 18 19 20 21 22 23 24 25) (
'halo_verifies' 'halos_communications' 26 27 28 29 30) ('incref'
'halos_allocates' 31) ('node_count' 'halos_base' 32) (
'pending_communication' 'halos_debug' 33 34) ('nullify' 'halos_allocates'
35) ('reallocate' 'halos_allocates' 36) ('read_halos' 'halos_registration'
37 38) ('nodes_owned' 'halos_ownership' 39) ('reorder_halo' 'halos_repair'
40 41) ('node_owned' 'halos_ownership' 42) ('has_references'
'halos_allocates' 43) ('halo_universal_number' 'halos_numbering' 44 45)
('serial_storage_halo' 'halos_base' 46 47) ('zero' 'halos_base' 48))

(('mpi_fortran_argv_null' 49 0 0 'mpi_fortran_argv_null') (
'mpi_fortran_in_place' 50 0 0 'mpi_fortran_in_place') (
'mpi_fortran_status_ignore' 51 0 0 'mpi_fortran_status_ignore') (
'mpi_fortran_statuses_ignore' 52 0 0 'mpi_fortran_statuses_ignore') (
'petscfortran1' 53 0 0 'petscfortran1') ('petscfortran10' 54 0 0
'petscfortran10') ('petscfortran2' 55 0 0 'petscfortran2') (
'petscfortran3' 56 0 0 'petscfortran3') ('petscfortran4' 57 0 0
'petscfortran4') ('petscfortran5' 58 0 0 'petscfortran5') (
'petscfortran6' 59 0 0 'petscfortran6') ('mpi_fortran_argvs_null' 60 0 0
'mpi_fortran_argvs_null') ('petscfortran7' 61 0 0 'petscfortran7') (
'petscfortran8' 62 0 0 'petscfortran8') ('petscfortran9' 63 0 0
'petscfortran9') ('mpi_fortran_bottom' 64 0 0 'mpi_fortran_bottom') (
'mpi_fortran_errcodes_ignore' 65 0 0 'mpi_fortran_errcodes_ignore'))

()

()

(60 'mpi_argvs_null' 'mpi_interfaces' 'mpi_argvs_null' 1 ((VARIABLE
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 IN_COMMON) (REAL 8 0 0
REAL ()) 0 0 () () 0 () () () 0 0)
49 'mpi_argv_null' 'mpi_interfaces' 'mpi_argv_null' 1 ((VARIABLE
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 DIMENSION IN_COMMON) (
CHARACTER 1 0 0 CHARACTER ((CONSTANT (INTEGER 4 0 0 INTEGER ()) 0 '1')))
0 0 () (1 0 EXPLICIT (CONSTANT (INTEGER 4 0 0 INTEGER ()) 0 '1') (
CONSTANT (INTEGER 4 0 0 INTEGER ()) 0 '1')) 0 () () () 0 0)
2 'deallocate_halo_vector' 'halos_allocates' 'deallocate_halo_vector' 1
((PROCEDURE UNKNOWN-INTENT MODULE-PROC DECL UNKNOWN 0 0 SUBROUTINE
ALWAYS_EXPLICIT) (UNKNOWN 0 0 0 UNKNOWN ()) 66 0 (67) () 0 () () () 0 0)
43 'has_references_halo_type' 'halos_allocates' 'has_references_halo_type'
1 ((PROCEDURE UNKNOWN-INTENT MODULE-PROC DECL UNKNOWN 0 0 FUNCTION PURE)
(LOGICAL 4 0 0 LOGICAL ()) 68 0 (69) () 70 () () () 0 0)
36 'reallocate_halo' 'halos_allocates' 'reallocate_halo' 1 ((PROCEDURE
UNKNOWN-INTENT MODULE-PROC DECL UNKNOWN 0 0 SUBROUTINE ALWAYS_EXPLICIT)
(UNKNOWN 0 0 0 UNKNOWN ()) 71 0 (72 73 74) () 0 () () () 0 0)
35 'nullify_halo' 'halos_allocates' 'nullify_halo' 1 ((PROCEDURE
UNKNOWN-INTENT MODULE-PROC DECL UNKNOWN 0 0 SUBROUTINE) (UNKNOWN 0 0 0
UNKNOWN ()) 75 0 (76) () 0 () () () 0 0)
77 'refcount_type' 'reference_counting' 'refcount_type' 1 ((DERIVED
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 POINTER_COMP) (UNKNOWN 0
0 0 UNKNOWN ()) 0 0 () () 0 ((78 'prev' (DERIVED 77 0 0 DERIVED ()) () (
UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 POINTER)
UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0)) (79 'next' (DERIVED
77 0 0 DERIVED ()) () (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN
UNKNOWN 0 0 POINTER) UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0))
(80 'count' (INTEGER 4 0 0 INTEGER ()) () (UNKNOWN-FL UNKNOWN-INTENT
UNKNOWN-PROC UNKNOWN UNKNOWN 0 0) UNKNOWN-ACCESS (CONSTANT (INTEGER 4 0
0 INTEGER ()) 0 '0')) (81 'id' (INTEGER 4 0 0 INTEGER ()) () (
UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0)
UNKNOWN-ACCESS ()) (82 'name' (CHARACTER 1 0 0 CHARACTER ((CONSTANT (
INTEGER 4 0 0 INTEGER ()) 0 '101'))) () (UNKNOWN-FL UNKNOWN-INTENT
UNKNOWN-PROC UNKNOWN UNKNOWN 0 0) UNKNOWN-ACCESS ()) (83 'type' (
CHARACTER 1 0 0 CHARACTER ((CONSTANT (INTEGER 4 0 0 INTEGER ()) 0 '101')))
() (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0)
UNKNOWN-ACCESS ()) (84 'tagged' (LOGICAL 4 0 0 LOGICAL ()) () (
UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0)
UNKNOWN-ACCESS (CONSTANT (LOGICAL 4 0 0 LOGICAL ()) 0 0))) PUBLIC (() ()
() ()) () 0 0 25948645)
85 'integer_vector' 'futils' 'integer_vector' 1 ((DERIVED UNKNOWN-INTENT
UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 POINTER_COMP) (UNKNOWN 0 0 0 UNKNOWN ())
0 0 () () 0 ((86 'ptr' (INTEGER 4 0 0 INTEGER ()) (1 0 DEFERRED () ()) (
UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 DIMENSION
POINTER) UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0))) PUBLIC (()
() () ()) () 0 0 9661976)
31 'incref_halo_type' 'halos_allocates' 'incref_halo_type' 1 ((
PROCEDURE UNKNOWN-INTENT MODULE-PROC DECL UNKNOWN 0 0 SUBROUTINE
ALWAYS_EXPLICIT) (UNKNOWN 0 0 0 UNKNOWN ()) 87 0 (88) () 0 () () () 0 0)
51 'mpi_status_ignore' 'mpi_interfaces' 'mpi_status_ignore' 1 ((
VARIABLE UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 DIMENSION
IN_COMMON) (INTEGER 4 0 0 INTEGER ()) 0 0 () (1 0 EXPLICIT (CONSTANT (
INTEGER 4 0 0 INTEGER ()) 0 '1') (CONSTANT (INTEGER 4 0 0 INTEGER ()) 0
'5')) 0 () () () 0 0)
52 'mpi_statuses_ignore' 'mpi_interfaces' 'mpi_statuses_ignore' 1 ((
VARIABLE UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 IN_COMMON) (
REAL 8 0 0 REAL ()) 0 0 () () 0 () () () 0 0)
64 'mpi_bottom' 'mpi_interfaces' 'mpi_bottom' 1 ((VARIABLE
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 IN_COMMON) (INTEGER 4 0
0 INTEGER ()) 0 0 () () 0 () () () 0 0)
50 'mpi_in_place' 'mpi_interfaces' 'mpi_in_place' 1 ((VARIABLE
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 IN_COMMON) (INTEGER 4 0
0 INTEGER ()) 0 0 () () 0 () () () 0 0)
65 'mpi_errcodes_ignore' 'mpi_interfaces' 'mpi_errcodes_ignore' 1 ((
VARIABLE UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 DIMENSION
IN_COMMON) (INTEGER 4 0 0 INTEGER ()) 0 0 () (1 0 EXPLICIT (CONSTANT (
INTEGER 4 0 0 INTEGER ()) 0 '1') (CONSTANT (INTEGER 4 0 0 INTEGER ()) 0
'1')) 0 () () () 0 0)
34 'pending_communication_communicator' 'parallel_tools'
'pending_communication_communicator' 1 ((PROCEDURE UNKNOWN-INTENT
MODULE-PROC DECL UNKNOWN 0 0 FUNCTION ALWAYS_EXPLICIT) (LOGICAL 4 0 0
LOGICAL ()) 89 0 (90) () 91 () () () 0 0)
33 'pending_communication_halo' 'halos_debug' 'pending_communication_halo'
1 ((PROCEDURE UNKNOWN-INTENT MODULE-PROC DECL UNKNOWN 0 0 FUNCTION) (
LOGICAL 4 0 0 LOGICAL ()) 92 0 (93) () 94 () () () 0 0)
58 'petsc_null_double' 'petscsysdef' 'petsc_null_double' 1 ((VARIABLE
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 IN_COMMON) (REAL 8 0 0
REAL ()) 0 0 () () 0 () () () 0 0)
57 'petsc_null_scalar' 'petscsysdef' 'petsc_null_scalar' 1 ((VARIABLE
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 IN_COMMON) (REAL 8 0 0
REAL ()) 0 0 () () 0 () () () 0 0)
61 'petsc_null_truth' 'petscsysdef' 'petsc_null_truth' 1 ((VARIABLE
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 IN_COMMON) (LOGICAL 4 0
0 LOGICAL ()) 0 0 () () 0 () () () 0 0)
63 'petsc_comm_world' 'petscsysdef' 'petsc_comm_world' 1 ((VARIABLE
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 IN_COMMON) (INTEGER 4 0
0 INTEGER ()) 0 0 () () 0 () () () 0 0)
62 'petsc_null_object' 'petscsysdef' 'petsc_null_object' 1 ((VARIABLE
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 IN_COMMON) (INTEGER 8 0
0 INTEGER ()) 0 0 () () 0 () () () 0 0)
59 'petsc_null_real' 'petscsysdef' 'petsc_null_real' 1 ((VARIABLE
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 IN_COMMON) (REAL 8 0 0
REAL ()) 0 0 () () 0 () () () 0 0)
56 'petsc_null' 'petscsysdef' 'petsc_null' 1 ((VARIABLE UNKNOWN-INTENT
UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 IN_COMMON) (INTEGER 4 0 0 INTEGER ()) 0
0 () () 0 () () () 0 0)
54 'petsc_comm_self' 'petscsysdef' 'petsc_comm_self' 1 ((VARIABLE
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 IN_COMMON) (INTEGER 4 0
0 INTEGER ()) 0 0 () () 0 () () () 0 0)
53 'petsc_null_character' 'petscsysdef' 'petsc_null_character' 1 ((
VARIABLE UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN IMPLICIT-SAVE 0 0 IN_COMMON)
(CHARACTER 1 0 0 CHARACTER ((CONSTANT (INTEGER 4 0 0 INTEGER ()) 0 '80')))
0 0 () () 0 () () () 0 0)
55 'petsc_null_integer' 'petscsysdef' 'petsc_null_integer' 1 ((VARIABLE
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 IN_COMMON) (INTEGER 4 0
0 INTEGER ()) 0 0 () () 0 () () () 0 0)
95 'vector_field' 'fields_data_types' 'vector_field' 1 ((DERIVED
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 POINTER_COMP) (UNKNOWN 0
0 0 UNKNOWN ()) 0 0 () () 0 ((96 'val' (REAL 8 0 0 REAL ()) (2 0
DEFERRED () () () ()) (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN
UNKNOWN 0 0 DIMENSION POINTER) UNKNOWN-ACCESS ()) (97 'wrapped' (
LOGICAL 4 0 0 LOGICAL ()) () (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0) UNKNOWN-ACCESS (CONSTANT (LOGICAL 4 0 0 LOGICAL ())
0 1)) (98 'field_type' (INTEGER 4 0 0 INTEGER ()) () (UNKNOWN-FL
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0) UNKNOWN-ACCESS (
CONSTANT (INTEGER 4 0 0 INTEGER ()) 0 '0')) (99 'bc' (DERIVED 100 0 0
DERIVED ()) () (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0
0 POINTER) UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0)) (101 'name'
(CHARACTER 1 0 0 CHARACTER ((CONSTANT (INTEGER 4 0 0 INTEGER ()) 0 '101')))
() (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0)
UNKNOWN-ACCESS ()) (102 'dim' (INTEGER 4 0 0 INTEGER ()) () (UNKNOWN-FL
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0) UNKNOWN-ACCESS ()) (
103 'option_path' (CHARACTER 1 0 0 CHARACTER ((CONSTANT (INTEGER 4 0 0
INTEGER ()) 0 '8192'))) () (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0) UNKNOWN-ACCESS (CONSTANT (CHARACTER 1 0 0 CHARACTER
(())) 0 8192
'/uninitialised_path/                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            '))
(104 'mesh' (DERIVED 105 0 0 DERIVED ()) () (UNKNOWN-FL UNKNOWN-INTENT
UNKNOWN-PROC UNKNOWN UNKNOWN 0 0) UNKNOWN-ACCESS (STRUCTURE (DERIVED 105
0 0 DERIVED ()) 0 ((() ()) ((CONSTANT (LOGICAL 4 0 0 LOGICAL ()) 0 1) ())
((STRUCTURE (DERIVED 106 0 0 DERIVED ()) 0 ((() ()) (() ()) (() ()) (()
()) ((NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0) ()) ((NULL (UNKNOWN 0 0 0
UNKNOWN ()) 0) ()) ((NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0) ()) ((NULL (
UNKNOWN 0 0 0 UNKNOWN ()) 0) ()) ((NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0) ())
((NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0) ()) ((NULL (UNKNOWN 0 0 0 UNKNOWN ())
0) ()) ((STRUCTURE (DERIVED 107 0 0 DERIVED ()) 0 ((() ()) (() ()) (() ())
(() ()) ((NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0) ()) ((NULL (UNKNOWN 0 0 0
UNKNOWN ()) 0) ()) (() ()) ((NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0) ()) (()
())) ()) ()) ((NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0) ()) ((NULL (UNKNOWN 0
0 0 UNKNOWN ()) 0) ()) ((NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0) ()) ((NULL (
UNKNOWN 0 0 0 UNKNOWN ()) 0) ()) (() ())) ()) ()) (() ()) (() ()) (() ())
((CONSTANT (CHARACTER 1 0 0 CHARACTER (())) 0 8192
'/uninitialised_path/                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            ')
()) ((CONSTANT (INTEGER 4 0 0 INTEGER ()) 0 '0') ()) ((NULL (UNKNOWN 0 0
0 UNKNOWN ()) 0) ()) ((NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0) ()) ((NULL (
UNKNOWN 0 0 0 UNKNOWN ()) 0) ()) ((NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0) ())
((NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0) ()) ((NULL (UNKNOWN 0 0 0 UNKNOWN ())
0) ()) ((NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0) ()) ((NULL (UNKNOWN 0 0 0
UNKNOWN ()) 0) ()) ((NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0) ()) ((NULL (
UNKNOWN 0 0 0 UNKNOWN ()) 0) ()) ((CONSTANT (LOGICAL 4 0 0 LOGICAL ()) 0
0) ())) ())) (108 'refcount' (DERIVED 77 0 0 DERIVED ()) () (UNKNOWN-FL
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 POINTER) UNKNOWN-ACCESS
(NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0)) (109 'aliased' (LOGICAL 4 0 0
LOGICAL ()) () (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0
0) UNKNOWN-ACCESS (CONSTANT (LOGICAL 4 0 0 LOGICAL ()) 0 0)) (110 'picker'
(DERIVED 111 0 0 DERIVED ()) () (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0 POINTER) UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN
()) 0))) PUBLIC (() () () ()) () 0 0 43598963)
105 'mesh_type' 'fields_data_types' 'mesh_type' 1 ((DERIVED
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 POINTER_COMP) (UNKNOWN 0
0 0 UNKNOWN ()) 0 0 () () 0 ((112 'ndglno' (INTEGER 4 0 0 INTEGER ()) (
1 0 DEFERRED () ()) (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN
UNKNOWN 0 0 DIMENSION POINTER) UNKNOWN-ACCESS ()) (113 'wrapped' (
LOGICAL 4 0 0 LOGICAL ()) () (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0) UNKNOWN-ACCESS (CONSTANT (LOGICAL 4 0 0 LOGICAL ())
0 1)) (114 'shape' (DERIVED 106 0 0 DERIVED ()) () (UNKNOWN-FL
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0) UNKNOWN-ACCESS (
STRUCTURE (DERIVED 106 0 0 DERIVED ()) 0 ((() ()) (() ()) (() ()) (() ())
((NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0) ()) ((NULL (UNKNOWN 0 0 0 UNKNOWN ())
0) ()) ((NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0) ()) ((NULL (UNKNOWN 0 0 0
UNKNOWN ()) 0) ()) ((NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0) ()) ((NULL (
UNKNOWN 0 0 0 UNKNOWN ()) 0) ()) ((NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0) ())
((STRUCTURE (DERIVED 107 0 0 DERIVED ()) 0 ((() ()) (() ()) (() ()) (()
()) ((NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0) ()) ((NULL (UNKNOWN 0 0 0
UNKNOWN ()) 0) ()) (() ()) ((NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0) ()) (()
())) ()) ()) ((NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0) ()) ((NULL (UNKNOWN 0
0 0 UNKNOWN ()) 0) ()) ((NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0) ()) ((NULL (
UNKNOWN 0 0 0 UNKNOWN ()) 0) ()) (() ())) ())) (115 'elements' (INTEGER
4 0 0 INTEGER ()) () (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN
UNKNOWN 0 0) UNKNOWN-ACCESS ()) (116 'nodes' (INTEGER 4 0 0 INTEGER ())
() (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0)
UNKNOWN-ACCESS ()) (117 'name' (CHARACTER 1 0 0 CHARACTER ((CONSTANT (
INTEGER 4 0 0 INTEGER ()) 0 '101'))) () (UNKNOWN-FL UNKNOWN-INTENT
UNKNOWN-PROC UNKNOWN UNKNOWN 0 0) UNKNOWN-ACCESS ()) (118 'option_path'
(CHARACTER 1 0 0 CHARACTER ((CONSTANT (INTEGER 4 0 0 INTEGER ()) 0 '8192')))
() (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0)
UNKNOWN-ACCESS (CONSTANT (CHARACTER 1 0 0 CHARACTER (())) 0 8192
'/uninitialised_path/                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            '))
(119 'continuity' (INTEGER 4 0 0 INTEGER ()) () (UNKNOWN-FL
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0) UNKNOWN-ACCESS (
CONSTANT (INTEGER 4 0 0 INTEGER ()) 0 '0')) (120 'refcount' (DERIVED 77
0 0 DERIVED ()) () (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN
UNKNOWN 0 0 POINTER) UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0))
(121 'faces' (DERIVED 122 0 0 DERIVED ()) () (UNKNOWN-FL UNKNOWN-INTENT
UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 POINTER) UNKNOWN-ACCESS (NULL (UNKNOWN
0 0 0 UNKNOWN ()) 0)) (123 'subdomain_mesh' (DERIVED 124 0 0 DERIVED ())
() (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 POINTER)
UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0)) (125 'adj_lists' (
DERIVED 126 0 0 DERIVED ()) () (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0 POINTER) UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN
()) 0)) (127 'columns' (INTEGER 4 0 0 INTEGER ()) (1 0 DEFERRED () ()) (
UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 DIMENSION
POINTER) UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0)) (128
'element_columns' (INTEGER 4 0 0 INTEGER ()) (1 0 DEFERRED () ()) (
UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 DIMENSION
POINTER) UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0)) (129
'region_ids' (INTEGER 4 0 0 INTEGER ()) (1 0 DEFERRED () ()) (
UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 DIMENSION
POINTER) UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0)) (130 'halos'
(DERIVED 131 0 0 DERIVED ()) (1 0 DEFERRED () ()) (UNKNOWN-FL
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 DIMENSION POINTER)
UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0)) (132 'element_halos'
(DERIVED 131 0 0 DERIVED ()) (1 0 DEFERRED () ()) (UNKNOWN-FL
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 DIMENSION POINTER)
UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0)) (133 'colourings' (
DERIVED 134 0 0 DERIVED ()) (1 0 DEFERRED () ()) (UNKNOWN-FL
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 DIMENSION POINTER)
UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0)) (135 'periodic' (
LOGICAL 4 0 0 LOGICAL ()) () (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0) UNKNOWN-ACCESS (CONSTANT (LOGICAL 4 0 0 LOGICAL ())
0 0))) PUBLIC (() () () ()) () 0 0 2572791)
100 'vector_boundary_conditions_ptr' 'fields_data_types'
'vector_boundary_conditions_ptr' 1 ((DERIVED UNKNOWN-INTENT UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0 POINTER_COMP) (UNKNOWN 0 0 0 UNKNOWN ()) 0 0 () () 0
((136 'boundary_condition' (DERIVED 137 0 0 DERIVED ()) (1 0 DEFERRED ()
()) (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DIMENSION POINTER) UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0)))
PUBLIC (() () () ()) () 0 0 84476181)
106 'element_type' 'elements' 'element_type' 1 ((DERIVED UNKNOWN-INTENT
UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 POINTER_COMP) (UNKNOWN 0 0 0 UNKNOWN ())
0 0 () () 0 ((138 'dim' (INTEGER 4 0 0 INTEGER ()) () (UNKNOWN-FL
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0) UNKNOWN-ACCESS ()) (
139 'loc' (INTEGER 4 0 0 INTEGER ()) () (UNKNOWN-FL UNKNOWN-INTENT
UNKNOWN-PROC UNKNOWN UNKNOWN 0 0) UNKNOWN-ACCESS ()) (140 'ngi' (
INTEGER 4 0 0 INTEGER ()) () (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0) UNKNOWN-ACCESS ()) (141 'degree' (INTEGER 4 0 0
INTEGER ()) () (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0
0) UNKNOWN-ACCESS ()) (142 'n' (REAL 8 0 0 REAL ()) (2 0 DEFERRED () ()
() ()) (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DIMENSION POINTER) UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0)) (
143 'dn' (REAL 8 0 0 REAL ()) (3 0 DEFERRED () () () () () ()) (
UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 DIMENSION
POINTER) UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0)) (144 'n_s'
(REAL 8 0 0 REAL ()) (3 0 DEFERRED () () () () () ()) (UNKNOWN-FL
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 DIMENSION POINTER)
UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0)) (145 'dn_s' (REAL 8
0 0 REAL ()) (4 0 DEFERRED () () () () () () () ()) (UNKNOWN-FL
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 DIMENSION POINTER)
UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0)) (146 'spoly' (
DERIVED 147 0 0 DERIVED ()) (2 0 DEFERRED () () () ()) (UNKNOWN-FL
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 DIMENSION POINTER)
UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0)) (148 'dspoly' (
DERIVED 147 0 0 DERIVED ()) (2 0 DEFERRED () () () ()) (UNKNOWN-FL
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 DIMENSION POINTER)
UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0)) (149 'numbering' (
DERIVED 150 0 0 DERIVED ()) () (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0 POINTER) UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN
()) 0)) (151 'quadrature' (DERIVED 107 0 0 DERIVED ()) () (UNKNOWN-FL
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0) UNKNOWN-ACCESS (
STRUCTURE (DERIVED 107 0 0 DERIVED ()) 0 ((() ()) (() ()) (() ()) (() ())
((NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0) ()) ((NULL (UNKNOWN 0 0 0 UNKNOWN ())
0) ()) (() ()) ((NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0) ()) (() ())) ())) (
152 'surface_quadrature' (DERIVED 107 0 0 DERIVED ()) () (UNKNOWN-FL
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 POINTER) UNKNOWN-ACCESS
(NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0)) (153 'superconvergence' (DERIVED
154 0 0 DERIVED ()) () (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN
UNKNOWN 0 0 POINTER) UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0))
(155 'constraints' (DERIVED 156 0 0 DERIVED ()) () (UNKNOWN-FL
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 POINTER) UNKNOWN-ACCESS
(NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0)) (157 'refcount' (DERIVED 77 0 0
DERIVED ()) () (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0
0 POINTER) UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0)) (158 'name'
(CHARACTER 1 0 0 CHARACTER ((CONSTANT (INTEGER 4 0 0 INTEGER ()) 0 '0')))
() (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0)
UNKNOWN-ACCESS ())) PUBLIC (() () () ()) () 0 0 79461029)
107 'quadrature_type' 'quadrature' 'quadrature_type' 1 ((DERIVED
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 POINTER_COMP) (UNKNOWN 0
0 0 UNKNOWN ()) 0 0 () () 0 ((159 'dim' (INTEGER 4 0 0 INTEGER ()) () (
UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0)
UNKNOWN-ACCESS ()) (160 'degree' (INTEGER 4 0 0 INTEGER ()) () (
UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0)
UNKNOWN-ACCESS ()) (161 'vertices' (INTEGER 4 0 0 INTEGER ()) () (
UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0)
UNKNOWN-ACCESS ()) (162 'ngi' (INTEGER 4 0 0 INTEGER ()) () (UNKNOWN-FL
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0) UNKNOWN-ACCESS ()) (
163 'weight' (REAL 8 0 0 REAL ()) (1 0 DEFERRED () ()) (UNKNOWN-FL
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 DIMENSION POINTER)
UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0)) (164 'l' (REAL 8 0 0
REAL ()) (2 0 DEFERRED () () () ()) (UNKNOWN-FL UNKNOWN-INTENT
UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 DIMENSION POINTER) UNKNOWN-ACCESS (
NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0)) (165 'name' (CHARACTER 1 0 0
CHARACTER ((CONSTANT (INTEGER 4 0 0 INTEGER ()) 0 '0'))) () (UNKNOWN-FL
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0) UNKNOWN-ACCESS ()) (
166 'refcount' (DERIVED 77 0 0 DERIVED ()) () (UNKNOWN-FL UNKNOWN-INTENT
UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 POINTER) UNKNOWN-ACCESS (NULL (UNKNOWN
0 0 0 UNKNOWN ()) 0)) (167 'family' (INTEGER 4 0 0 INTEGER ()) () (
UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0)
UNKNOWN-ACCESS ())) PUBLIC (() () () ()) () 0 0 59837722)
111 'picker_ptr' 'picker_data_types' 'picker_ptr' 1 ((DERIVED
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 POINTER_COMP) (UNKNOWN 0
0 0 UNKNOWN ()) 0 0 () () 0 ((168 'ptr' (DERIVED 169 0 0 DERIVED ()) ()
(UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 POINTER)
UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0))) PUBLIC (() () () ())
() 0 0 71120007)
122 'mesh_faces' 'fields_data_types' 'mesh_faces' 1 ((DERIVED
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 POINTER_COMP) (UNKNOWN 0
0 0 UNKNOWN ()) 0 0 () () 0 ((170 'shape' (DERIVED 106 0 0 DERIVED ()) ()
(UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 POINTER)
UNKNOWN-ACCESS ()) (171 'face_list' (DERIVED 172 0 0 DERIVED ()) () (
UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0)
UNKNOWN-ACCESS (STRUCTURE (DERIVED 172 0 0 DERIVED ()) 0 (((STRUCTURE (
DERIVED 173 0 0 DERIVED ()) 0 (((NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0) ())
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
())) ())) (174 'face_lno' (INTEGER 4 0 0 INTEGER ()) (1 0 DEFERRED () ())
(UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 DIMENSION
POINTER) UNKNOWN-ACCESS ()) (175 'surface_mesh' (DERIVED 105 0 0 DERIVED
()) () (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0)
UNKNOWN-ACCESS (STRUCTURE (DERIVED 105 0 0 DERIVED ()) 0 ((() ()) ((
CONSTANT (LOGICAL 4 0 0 LOGICAL ()) 0 1) ()) ((STRUCTURE (DERIVED 106 0
0 DERIVED ()) 0 ((() ()) (() ()) (() ()) (() ()) ((NULL (UNKNOWN 0 0 0
UNKNOWN ()) 0) ()) ((NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0) ()) ((NULL (
UNKNOWN 0 0 0 UNKNOWN ()) 0) ()) ((NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0) ())
((NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0) ()) ((NULL (UNKNOWN 0 0 0 UNKNOWN ())
0) ()) ((NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0) ()) ((STRUCTURE (DERIVED 107
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
UNKNOWN ()) 0) ()) ((NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0) ()) ((NULL (
UNKNOWN 0 0 0 UNKNOWN ()) 0) ()) ((CONSTANT (LOGICAL 4 0 0 LOGICAL ()) 0
0) ())) ())) (176 'surface_node_list' (INTEGER 4 0 0 INTEGER ()) (1 0
DEFERRED () ()) (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN
0 0 DIMENSION POINTER) UNKNOWN-ACCESS ()) (177 'face_element_list' (
INTEGER 4 0 0 INTEGER ()) (1 0 DEFERRED () ()) (UNKNOWN-FL
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 DIMENSION POINTER)
UNKNOWN-ACCESS ()) (178 'boundary_ids' (INTEGER 4 0 0 INTEGER ()) (1 0
DEFERRED () ()) (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN
0 0 DIMENSION POINTER) UNKNOWN-ACCESS ()) (179 'coplanar_ids' (INTEGER 4
0 0 INTEGER ()) (1 0 DEFERRED () ()) (UNKNOWN-FL UNKNOWN-INTENT
UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 DIMENSION POINTER) UNKNOWN-ACCESS (
NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0)) (180 'dg_surface_mesh' (DERIVED 105
0 0 DERIVED ()) () (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN
UNKNOWN 0 0 POINTER) UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0))
(181 'has_internal_boundaries' (LOGICAL 4 0 0 LOGICAL ()) () (
UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0)
UNKNOWN-ACCESS (CONSTANT (LOGICAL 4 0 0 LOGICAL ()) 0 0))) PUBLIC (() ()
() ()) () 0 0 11936185)
126 'adjacency_cache' 'fields_data_types' 'adjacency_cache' 1 ((DERIVED
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 POINTER_COMP) (UNKNOWN 0
0 0 UNKNOWN ()) 0 0 () () 0 ((182 'nnlist' (DERIVED 173 0 0 DERIVED ())
() (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 POINTER)
UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0)) (183 'nelist' (
DERIVED 173 0 0 DERIVED ()) () (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0 POINTER) UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN
()) 0)) (184 'eelist' (DERIVED 173 0 0 DERIVED ()) () (UNKNOWN-FL
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 POINTER) UNKNOWN-ACCESS
(NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0))) PUBLIC (() () () ()) () 0 0
37419158)
124 'mesh_subdomain_mesh' 'fields_data_types' 'mesh_subdomain_mesh' 1 (
(DERIVED UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 POINTER_COMP) (
UNKNOWN 0 0 0 UNKNOWN ()) 0 0 () () 0 ((185 'element_list' (INTEGER 4 0
0 INTEGER ()) (1 0 DEFERRED () ()) (UNKNOWN-FL UNKNOWN-INTENT
UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 DIMENSION POINTER) UNKNOWN-ACCESS ()) (
186 'node_list' (INTEGER 4 0 0 INTEGER ()) (1 0 DEFERRED () ()) (
UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 DIMENSION
POINTER) UNKNOWN-ACCESS ())) PUBLIC (() () () ()) () 0 0 85302181)
134 'integer_set_vector' 'integer_set_module' 'integer_set_vector' 1 ((
DERIVED UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 POINTER_COMP) (
UNKNOWN 0 0 0 UNKNOWN ()) 0 0 () () 0 ((187 'sets' (DERIVED 188 0 0
DERIVED ()) (1 0 DEFERRED () ()) (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0 DIMENSION POINTER) UNKNOWN-ACCESS ())) PUBLIC (() ()
() ()) () 0 0 40688950)
3 'deallocate_halo' 'halos_allocates' 'deallocate_halo' 1 ((PROCEDURE
UNKNOWN-INTENT MODULE-PROC DECL UNKNOWN 0 0 SUBROUTINE) (UNKNOWN 0 0 0
UNKNOWN ()) 189 0 (190) () 0 () () () 0 0)
169 'picker_type' 'picker_data_types' 'picker_type' 1 ((DERIVED
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 POINTER_COMP) (UNKNOWN 0
0 0 UNKNOWN ()) 0 0 () () 0 ((191 'name' (CHARACTER 1 0 0 CHARACTER ((
CONSTANT (INTEGER 4 0 0 INTEGER ()) 0 '101'))) () (UNKNOWN-FL
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0) UNKNOWN-ACCESS ()) (
192 'refcount' (DERIVED 77 0 0 DERIVED ()) () (UNKNOWN-FL UNKNOWN-INTENT
UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 POINTER) UNKNOWN-ACCESS (NULL (UNKNOWN
0 0 0 UNKNOWN ()) 0)) (193 'picker_id' (INTEGER 4 0 0 INTEGER ()) () (
UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0)
UNKNOWN-ACCESS (CONSTANT (INTEGER 4 0 0 INTEGER ()) 0 '0')) (194
'last_mesh_movement' (INTEGER 4 0 0 INTEGER ()) () (UNKNOWN-FL
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0) UNKNOWN-ACCESS (
CONSTANT (INTEGER 4 0 0 INTEGER ()) 0 '0'))) PUBLIC (() () () ()) () 0 0
8821665)
173 'csr_sparsity' 'sparse_tools' 'csr_sparsity' 1 ((DERIVED
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 POINTER_COMP) (UNKNOWN 0
0 0 UNKNOWN ()) 0 0 () () 0 ((195 'findrm' (INTEGER 4 0 0 INTEGER ()) (
1 0 DEFERRED () ()) (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN
UNKNOWN 0 0 DIMENSION POINTER) UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0
UNKNOWN ()) 0)) (196 'centrm' (INTEGER 4 0 0 INTEGER ()) (1 0 DEFERRED ()
()) (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DIMENSION POINTER) UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0)) (
197 'colm' (INTEGER 4 0 0 INTEGER ()) (1 0 DEFERRED () ()) (UNKNOWN-FL
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 DIMENSION POINTER)
UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0)) (198 'columns' (
INTEGER 4 0 0 INTEGER ()) () (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0) UNKNOWN-ACCESS ()) (199 'row_halo' (DERIVED 131 0 0
DERIVED ()) () (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0
0 POINTER) UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0)) (200
'column_halo' (DERIVED 131 0 0 DERIVED ()) () (UNKNOWN-FL UNKNOWN-INTENT
UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 POINTER) UNKNOWN-ACCESS (NULL (UNKNOWN
0 0 0 UNKNOWN ()) 0)) (201 'refcount' (DERIVED 77 0 0 DERIVED ()) () (
UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 POINTER)
UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0)) (202 'name' (
CHARACTER 1 0 0 CHARACTER ((CONSTANT (INTEGER 4 0 0 INTEGER ()) 0 '101')))
() (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0)
UNKNOWN-ACCESS (CONSTANT (CHARACTER 1 0 0 CHARACTER (())) 0 101
'                                                                                                     '))
(203 'wrapped' (LOGICAL 4 0 0 LOGICAL ()) () (UNKNOWN-FL UNKNOWN-INTENT
UNKNOWN-PROC UNKNOWN UNKNOWN 0 0) UNKNOWN-ACCESS (CONSTANT (LOGICAL 4 0
0 LOGICAL ()) 0 0)) (204 'sorted_rows' (LOGICAL 4 0 0 LOGICAL ()) () (
UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0)
UNKNOWN-ACCESS (CONSTANT (LOGICAL 4 0 0 LOGICAL ()) 0 0))) PUBLIC (() ()
() ()) () 0 0 78882345)
172 'csr_matrix' 'sparse_tools' 'csr_matrix' 1 ((DERIVED UNKNOWN-INTENT
UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 POINTER_COMP) (UNKNOWN 0 0 0 UNKNOWN ())
0 0 () () 0 ((205 'sparsity' (DERIVED 173 0 0 DERIVED ()) () (
UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0)
UNKNOWN-ACCESS (STRUCTURE (DERIVED 173 0 0 DERIVED ()) 0 (((NULL (
UNKNOWN 0 0 0 UNKNOWN ()) 0) ()) ((NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0) ())
((NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0) ()) (() ()) ((NULL (UNKNOWN 0 0 0
UNKNOWN ()) 0) ()) ((NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0) ()) ((NULL (
UNKNOWN 0 0 0 UNKNOWN ()) 0) ()) ((CONSTANT (CHARACTER 1 0 0 CHARACTER (
())) 0 101
'                                                                                                     ')
()) ((CONSTANT (LOGICAL 4 0 0 LOGICAL ()) 0 0) ()) ((CONSTANT (LOGICAL 4
0 0 LOGICAL ()) 0 0) ())) ())) (206 'val' (REAL 8 0 0 REAL ()) (1 0
DEFERRED () ()) (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN
0 0 DIMENSION POINTER) UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0))
(207 'ival' (INTEGER 4 0 0 INTEGER ()) (1 0 DEFERRED () ()) (UNKNOWN-FL
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 DIMENSION POINTER)
UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0)) (208 'clone' (
LOGICAL 4 0 0 LOGICAL ()) () (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0) UNKNOWN-ACCESS (CONSTANT (LOGICAL 4 0 0 LOGICAL ())
0 0)) (209 'external_val' (LOGICAL 4 0 0 LOGICAL ()) () (UNKNOWN-FL
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0) UNKNOWN-ACCESS (
CONSTANT (LOGICAL 4 0 0 LOGICAL ()) 0 0)) (210 'inactive' (DERIVED 211 0
0 DERIVED ()) () (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN
0 0 POINTER) UNKNOWN-ACCESS ()) (212 'ksp' (INTEGER 8 0 0 INTEGER ()) ()
(UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 POINTER)
UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0)) (213 'refcount' (
DERIVED 77 0 0 DERIVED ()) () (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0 POINTER) UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN
()) 0)) (214 'name' (CHARACTER 1 0 0 CHARACTER ((CONSTANT (INTEGER 4 0 0
INTEGER ()) 0 '101'))) () (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0) UNKNOWN-ACCESS (CONSTANT (CHARACTER 1 0 0 CHARACTER
(())) 0 101
'                                                                                                     ')))
PUBLIC (() () () ()) () 0 0 66371681)
211 'logical_array_ptr' 'sparse_tools' 'logical_array_ptr' 1 ((DERIVED
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 POINTER_COMP) (UNKNOWN 0
0 0 UNKNOWN ()) 0 0 () () 0 ((215 'ptr' (LOGICAL 4 0 0 LOGICAL ()) (1 0
DEFERRED () ()) (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN
0 0 DIMENSION POINTER) UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0)))
PUBLIC (() () () ()) () 0 0 30974511)
188 'integer_set' 'integer_set_module' 'integer_set' 1 ((DERIVED
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0) (UNKNOWN 0 0 0 UNKNOWN
()) 0 0 () () 0 ((216 'address' (DERIVED 217 0 0 DERIVED ()) () (
UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0)
UNKNOWN-ACCESS ())) PUBLIC (() () () ()) () 0 0 67217964)
217 'c_ptr' '__iso_c_binding' '' 1 ((DERIVED UNKNOWN-INTENT UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0 IS_BIND_C IS_C_INTEROP IS_ISO_C) (DERIVED 217 1 1
UNKNOWN ()) 0 0 () () 0 ((218 '__c_ptr_c_address' (INTEGER 8 1 0 INTEGER
()) () (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0)
UNKNOWN-ACCESS ())) UNKNOWN-ACCESS () () 2 39 0)
137 'vector_boundary_condition' 'fields_data_types'
'vector_boundary_condition' 1 ((DERIVED UNKNOWN-INTENT UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0 POINTER_COMP) (UNKNOWN 0 0 0 UNKNOWN ()) 0 0 () () 0
((219 'name' (CHARACTER 1 0 0 CHARACTER ((CONSTANT (INTEGER 4 0 0
INTEGER ()) 0 '101'))) () (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0) UNKNOWN-ACCESS ()) (220 'type' (CHARACTER 1 0 0
CHARACTER ((CONSTANT (INTEGER 4 0 0 INTEGER ()) 0 '101'))) () (
UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0)
UNKNOWN-ACCESS (CONSTANT (CHARACTER 1 0 0 CHARACTER (())) 0 101
'                                                                                                     '))
(221 'applies' (LOGICAL 4 0 0 LOGICAL ()) (1 0 EXPLICIT (CONSTANT (
INTEGER 4 0 0 INTEGER ()) 0 '1') (CONSTANT (INTEGER 4 0 0 INTEGER ()) 0
'3')) (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DIMENSION) UNKNOWN-ACCESS ()) (222 'surface_element_list' (INTEGER 4 0 0
INTEGER ()) (1 0 DEFERRED () ()) (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0 DIMENSION POINTER) UNKNOWN-ACCESS (NULL (UNKNOWN 0 0
0 UNKNOWN ()) 0)) (223 'surface_node_list' (INTEGER 4 0 0 INTEGER ()) (
1 0 DEFERRED () ()) (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN
UNKNOWN 0 0 DIMENSION POINTER) UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0
UNKNOWN ()) 0)) (224 'surface_mesh' (DERIVED 105 0 0 DERIVED ()) () (
UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 POINTER)
UNKNOWN-ACCESS ()) (225 'surface_fields' (DERIVED 95 0 0 DERIVED ()) (1
0 DEFERRED () ()) (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN
UNKNOWN 0 0 DIMENSION POINTER) UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0
UNKNOWN ()) 0)) (226 'scalar_surface_fields' (DERIVED 227 0 0 DERIVED ())
(1 0 DEFERRED () ()) (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN
UNKNOWN 0 0 DIMENSION POINTER) UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0
UNKNOWN ()) 0)) (228 'option_path' (CHARACTER 1 0 0 CHARACTER ((
CONSTANT (INTEGER 4 0 0 INTEGER ()) 0 '8192'))) () (UNKNOWN-FL
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0) UNKNOWN-ACCESS (
CONSTANT (CHARACTER 1 0 0 CHARACTER (())) 0 8192
'/uninitialised_path/                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            ')))
PUBLIC (() () () ()) () 0 0 42205165)
229 'scalar_boundary_conditions_ptr' 'fields_data_types'
'scalar_boundary_conditions_ptr' 1 ((DERIVED UNKNOWN-INTENT UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0 POINTER_COMP) (UNKNOWN 0 0 0 UNKNOWN ()) 0 0 () () 0
((230 'boundary_condition' (DERIVED 231 0 0 DERIVED ()) (1 0 DEFERRED ()
()) (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DIMENSION POINTER) UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0)))
PUBLIC (() () () ()) () 0 0 47502942)
231 'scalar_boundary_condition' 'fields_data_types'
'scalar_boundary_condition' 1 ((DERIVED UNKNOWN-INTENT UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0 POINTER_COMP) (UNKNOWN 0 0 0 UNKNOWN ()) 0 0 () () 0
((232 'name' (CHARACTER 1 0 0 CHARACTER ((CONSTANT (INTEGER 4 0 0
INTEGER ()) 0 '101'))) () (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0) UNKNOWN-ACCESS ()) (233 'type' (CHARACTER 1 0 0
CHARACTER ((CONSTANT (INTEGER 4 0 0 INTEGER ()) 0 '101'))) () (
UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0)
UNKNOWN-ACCESS (CONSTANT (CHARACTER 1 0 0 CHARACTER (())) 0 101
'                                                                                                     '))
(234 'surface_element_list' (INTEGER 4 0 0 INTEGER ()) (1 0 DEFERRED ()
()) (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DIMENSION POINTER) UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0)) (
235 'surface_node_list' (INTEGER 4 0 0 INTEGER ()) (1 0 DEFERRED () ())
(UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 DIMENSION
POINTER) UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0)) (236
'surface_mesh' (DERIVED 105 0 0 DERIVED ()) () (UNKNOWN-FL
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 POINTER) UNKNOWN-ACCESS
()) (237 'surface_fields' (DERIVED 227 0 0 DERIVED ()) (1 0 DEFERRED ()
()) (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DIMENSION POINTER) UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0)) (
238 'option_path' (CHARACTER 1 0 0 CHARACTER ((CONSTANT (INTEGER 4 0 0
INTEGER ()) 0 '8192'))) () (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0) UNKNOWN-ACCESS (CONSTANT (CHARACTER 1 0 0 CHARACTER
(())) 0 8192
'/uninitialised_path/                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            ')))
PUBLIC (() () () ()) () 0 0 50740068)
8 'allocate_halo' 'halos_allocates' 'allocate_halo' 1 ((PROCEDURE
UNKNOWN-INTENT MODULE-PROC DECL UNKNOWN 0 0 SUBROUTINE ALWAYS_EXPLICIT)
(UNKNOWN 0 0 0 UNKNOWN ()) 239 0 (240 241 242 243 244 245 246 247 248) ()
0 () () () 0 0)
147 'polynomial' 'polynomials' 'polynomial' 1 ((DERIVED UNKNOWN-INTENT
UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 POINTER_COMP) (UNKNOWN 0 0 0 UNKNOWN ())
0 0 () () 0 ((249 'coefs' (REAL 8 0 0 REAL ()) (1 0 DEFERRED () ()) (
UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 DIMENSION
POINTER) UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0)) (250 'degree'
(INTEGER 4 0 0 INTEGER ()) () (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0) UNKNOWN-ACCESS (CONSTANT (INTEGER 4 0 0 INTEGER ())
0 '-1'))) PUBLIC (() () () ()) () 0 0 87989236)
150 'ele_numbering_type' 'element_numbering' 'ele_numbering_type' 1 ((
DERIVED UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 POINTER_COMP) (
UNKNOWN 0 0 0 UNKNOWN ()) 0 0 () () 0 ((251 'faces' (INTEGER 4 0 0
INTEGER ()) () (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0
0) UNKNOWN-ACCESS ()) (252 'vertices' (INTEGER 4 0 0 INTEGER ()) () (
UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0)
UNKNOWN-ACCESS ()) (253 'edges' (INTEGER 4 0 0 INTEGER ()) () (
UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0)
UNKNOWN-ACCESS ()) (254 'boundaries' (INTEGER 4 0 0 INTEGER ()) () (
UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0)
UNKNOWN-ACCESS ()) (255 'degree' (INTEGER 4 0 0 INTEGER ()) () (
UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0)
UNKNOWN-ACCESS ()) (256 'dimension' (INTEGER 4 0 0 INTEGER ()) () (
UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0)
UNKNOWN-ACCESS ()) (257 'nodes' (INTEGER 4 0 0 INTEGER ()) () (
UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0)
UNKNOWN-ACCESS ()) (258 'type' (INTEGER 4 0 0 INTEGER ()) () (
UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0)
UNKNOWN-ACCESS (CONSTANT (INTEGER 4 0 0 INTEGER ()) 0 '1')) (259 'family'
(INTEGER 4 0 0 INTEGER ()) () (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0) UNKNOWN-ACCESS ()) (260 'count2number' (INTEGER 4 0
0 INTEGER ()) (3 0 DEFERRED () () () () () ()) (UNKNOWN-FL
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 DIMENSION POINTER)
UNKNOWN-ACCESS ()) (261 'number2count' (INTEGER 4 0 0 INTEGER ()) (2 0
DEFERRED () () () ()) (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN
UNKNOWN 0 0 DIMENSION POINTER) UNKNOWN-ACCESS ()) (262 'boundary_coord'
(INTEGER 4 0 0 INTEGER ()) (1 0 DEFERRED () ()) (UNKNOWN-FL
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 DIMENSION POINTER)
UNKNOWN-ACCESS ()) (263 'boundary_val' (INTEGER 4 0 0 INTEGER ()) (1 0
DEFERRED () ()) (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN
0 0 DIMENSION POINTER) UNKNOWN-ACCESS ())) PUBLIC (() () () ()) () 0 0
96431082)
156 'constraints_type' 'elements' 'constraints_type' 1 ((DERIVED
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 POINTER_COMP) (UNKNOWN 0
0 0 UNKNOWN ()) 0 0 () () 0 ((264 'type' (INTEGER 4 0 0 INTEGER ()) () (
UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0)
UNKNOWN-ACCESS ()) (265 'dim' (INTEGER 4 0 0 INTEGER ()) () (UNKNOWN-FL
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0) UNKNOWN-ACCESS ()) (
266 'degree' (INTEGER 4 0 0 INTEGER ()) () (UNKNOWN-FL UNKNOWN-INTENT
UNKNOWN-PROC UNKNOWN UNKNOWN 0 0) UNKNOWN-ACCESS ()) (267 'loc' (
INTEGER 4 0 0 INTEGER ()) () (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0) UNKNOWN-ACCESS ()) (268 'n_constraints' (INTEGER 4
0 0 INTEGER ()) () (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN
UNKNOWN 0 0) UNKNOWN-ACCESS ()) (269 'orthogonal' (REAL 8 0 0 REAL ()) (
3 0 DEFERRED () () () () () ()) (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0 DIMENSION POINTER) UNKNOWN-ACCESS (NULL (UNKNOWN 0 0
0 UNKNOWN ()) 0))) PUBLIC (() () () ()) () 0 0 57548971)
154 'superconvergence_type' 'elements' 'superconvergence_type' 1 ((
DERIVED UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 POINTER_COMP) (
UNKNOWN 0 0 0 UNKNOWN ()) 0 0 () () 0 ((270 'nsp' (INTEGER 4 0 0 INTEGER
()) () (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0)
UNKNOWN-ACCESS ()) (271 'l' (REAL 8 0 0 REAL ()) (2 0 DEFERRED () () ()
()) (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DIMENSION POINTER) UNKNOWN-ACCESS ()) (272 'n' (REAL 8 0 0 REAL ()) (2 0
DEFERRED () () () ()) (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN
UNKNOWN 0 0 DIMENSION POINTER) UNKNOWN-ACCESS ()) (273 'dn' (REAL 8 0 0
REAL ()) (3 0 DEFERRED () () () () () ()) (UNKNOWN-FL UNKNOWN-INTENT
UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 DIMENSION POINTER) UNKNOWN-ACCESS ()))
PUBLIC (() () () ()) () 0 0 18282395)
227 'scalar_field' 'fields_data_types' 'scalar_field' 1 ((DERIVED
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 POINTER_COMP) (UNKNOWN 0
0 0 UNKNOWN ()) 0 0 () () 0 ((274 'val' (REAL 8 0 0 REAL ()) (1 0
DEFERRED () ()) (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN
0 0 DIMENSION POINTER) UNKNOWN-ACCESS ()) (275 'val_stride' (INTEGER 4 0
0 INTEGER ()) () (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN
0 0) UNKNOWN-ACCESS (CONSTANT (INTEGER 4 0 0 INTEGER ()) 0 '1')) (276
'wrapped' (LOGICAL 4 0 0 LOGICAL ()) () (UNKNOWN-FL UNKNOWN-INTENT
UNKNOWN-PROC UNKNOWN UNKNOWN 0 0) UNKNOWN-ACCESS (CONSTANT (LOGICAL 4 0
0 LOGICAL ()) 0 1)) (277 'field_type' (INTEGER 4 0 0 INTEGER ()) () (
UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0)
UNKNOWN-ACCESS (CONSTANT (INTEGER 4 0 0 INTEGER ()) 0 '0')) (278 'bc' (
DERIVED 229 0 0 DERIVED ()) () (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0 POINTER) UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN
()) 0)) (279 'name' (CHARACTER 1 0 0 CHARACTER ((CONSTANT (INTEGER 4 0 0
INTEGER ()) 0 '101'))) () (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0) UNKNOWN-ACCESS ()) (280 'option_path' (CHARACTER 1
0 0 CHARACTER ((CONSTANT (INTEGER 4 0 0 INTEGER ()) 0 '8192'))) () (
UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0)
UNKNOWN-ACCESS (CONSTANT (CHARACTER 1 0 0 CHARACTER (())) 0 8192
'/uninitialised_path/                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            '))
(281 'mesh' (DERIVED 105 0 0 DERIVED ()) () (UNKNOWN-FL UNKNOWN-INTENT
UNKNOWN-PROC UNKNOWN UNKNOWN 0 0) UNKNOWN-ACCESS (STRUCTURE (DERIVED 105
0 0 DERIVED ()) 0 ((() ()) ((CONSTANT (LOGICAL 4 0 0 LOGICAL ()) 0 1) ())
((STRUCTURE (DERIVED 106 0 0 DERIVED ()) 0 ((() ()) (() ()) (() ()) (()
()) ((NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0) ()) ((NULL (UNKNOWN 0 0 0
UNKNOWN ()) 0) ()) ((NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0) ()) ((NULL (
UNKNOWN 0 0 0 UNKNOWN ()) 0) ()) ((NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0) ())
((NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0) ()) ((NULL (UNKNOWN 0 0 0 UNKNOWN ())
0) ()) ((STRUCTURE (DERIVED 107 0 0 DERIVED ()) 0 ((() ()) (() ()) (() ())
(() ()) ((NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0) ()) ((NULL (UNKNOWN 0 0 0
UNKNOWN ()) 0) ()) (() ()) ((NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0) ()) (()
())) ()) ()) ((NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0) ()) ((NULL (UNKNOWN 0
0 0 UNKNOWN ()) 0) ()) ((NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0) ()) ((NULL (
UNKNOWN 0 0 0 UNKNOWN ()) 0) ()) (() ())) ()) ()) (() ()) (() ()) (() ())
((CONSTANT (CHARACTER 1 0 0 CHARACTER (())) 0 8192
'/uninitialised_path/                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            ')
()) ((CONSTANT (INTEGER 4 0 0 INTEGER ()) 0 '0') ()) ((NULL (UNKNOWN 0 0
0 UNKNOWN ()) 0) ()) ((NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0) ()) ((NULL (
UNKNOWN 0 0 0 UNKNOWN ()) 0) ()) ((NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0) ())
((NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0) ()) ((NULL (UNKNOWN 0 0 0 UNKNOWN ())
0) ()) ((NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0) ()) ((NULL (UNKNOWN 0 0 0
UNKNOWN ()) 0) ()) ((NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0) ()) ((NULL (
UNKNOWN 0 0 0 UNKNOWN ()) 0) ()) ((CONSTANT (LOGICAL 4 0 0 LOGICAL ()) 0
0) ())) ())) (282 'refcount' (DERIVED 77 0 0 DERIVED ()) () (UNKNOWN-FL
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 POINTER) UNKNOWN-ACCESS
(NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0)) (283 'aliased' (LOGICAL 4 0 0
LOGICAL ()) () (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0
0) UNKNOWN-ACCESS (CONSTANT (LOGICAL 4 0 0 LOGICAL ()) 0 0)) (284
'py_locweight' (REAL 8 0 0 REAL ()) (2 0 DEFERRED () () () ()) (
UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 DIMENSION
POINTER) UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0)) (285
'py_func' (CHARACTER 1 0 0 CHARACTER ((CONSTANT (INTEGER 4 0 0 INTEGER ())
0 '8192'))) () (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0
0) UNKNOWN-ACCESS ()) (286 'py_positions' (DERIVED 95 0 0 DERIVED ()) ()
(UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 POINTER)
UNKNOWN-ACCESS ()) (287 'py_positions_same_mesh' (LOGICAL 4 0 0 LOGICAL
()) () (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0)
UNKNOWN-ACCESS ()) (288 'py_dim' (INTEGER 4 0 0 INTEGER ()) () (
UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0)
UNKNOWN-ACCESS ()) (289 'py_positions_shape' (DERIVED 106 0 0 DERIVED ())
() (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 POINTER)
UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0))) PUBLIC (() () () ())
() 0 0 64912956)
4 'derive_l1_from_l2_halo_mesh' 'halos_derivation'
'derive_l1_from_l2_halo_mesh' 1 ((PROCEDURE UNKNOWN-INTENT MODULE-PROC
DECL UNKNOWN 0 0 SUBROUTINE ALWAYS_EXPLICIT) (UNKNOWN 0 0 0 UNKNOWN ())
290 0 (291 292 293) () 0 () () () 0 0)
294 'integer_hash_table' 'integer_hash_table_module' 'integer_hash_table'
1 ((DERIVED UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0) (UNKNOWN 0
0 0 UNKNOWN ()) 0 0 () () 0 ((295 'address' (DERIVED 217 0 0 DERIVED ())
() (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0)
UNKNOWN-ACCESS ())) PUBLIC (() () () ()) () 0 0 67448272)
15 'halo_update_scalar' 'halos_communications' 'halo_update_scalar' 1 (
(PROCEDURE UNKNOWN-INTENT MODULE-PROC DECL UNKNOWN 0 0 SUBROUTINE
ALWAYS_EXPLICIT) (UNKNOWN 0 0 0 UNKNOWN ()) 296 0 (297 298) () 0 () () ()
0 0)
14 'halo_update_vector' 'halos_communications' 'halo_update_vector' 1 (
(PROCEDURE UNKNOWN-INTENT MODULE-PROC DECL UNKNOWN 0 0 SUBROUTINE
ALWAYS_EXPLICIT) (UNKNOWN 0 0 0 UNKNOWN ()) 299 0 (300 301) () 0 () () ()
0 0)
16 'halo_update_tensor_on_halo' 'halos_communications'
'halo_update_tensor_on_halo' 1 ((PROCEDURE UNKNOWN-INTENT MODULE-PROC
DECL UNKNOWN 0 0 SUBROUTINE) (UNKNOWN 0 0 0 UNKNOWN ()) 302 0 (303 304)
() 0 () () () 0 0)
13 'halo_update_tensor' 'halos_communications' 'halo_update_tensor' 1 (
(PROCEDURE UNKNOWN-INTENT MODULE-PROC DECL UNKNOWN 0 0 SUBROUTINE
ALWAYS_EXPLICIT) (UNKNOWN 0 0 0 UNKNOWN ()) 305 0 (306 307) () 0 () () ()
0 0)
18 'halo_update_scalar_on_halo' 'halos_communications'
'halo_update_scalar_on_halo' 1 ((PROCEDURE UNKNOWN-INTENT MODULE-PROC
DECL UNKNOWN 0 0 SUBROUTINE) (UNKNOWN 0 0 0 UNKNOWN ()) 308 0 (309 310)
() 0 () () () 0 0)
20 'halo_update_array_real_block' 'halos_communications'
'halo_update_array_real_block' 1 ((PROCEDURE UNKNOWN-INTENT MODULE-PROC
DECL UNKNOWN 0 0 SUBROUTINE ALWAYS_EXPLICIT) (UNKNOWN 0 0 0 UNKNOWN ())
311 0 (312 313) () 0 () () () 0 0)
19 'halo_update_array_real_block2' 'halos_communications'
'halo_update_array_real_block2' 1 ((PROCEDURE UNKNOWN-INTENT MODULE-PROC
DECL UNKNOWN 0 0 SUBROUTINE ALWAYS_EXPLICIT) (UNKNOWN 0 0 0 UNKNOWN ())
314 0 (315 316) () 0 () () () 0 0)
17 'halo_update_vector_on_halo' 'halos_communications'
'halo_update_vector_on_halo' 1 ((PROCEDURE UNKNOWN-INTENT MODULE-PROC
DECL UNKNOWN 0 0 SUBROUTINE) (UNKNOWN 0 0 0 UNKNOWN ()) 317 0 (318 319)
() 0 () () () 0 0)
22 'halo_update_array_integer_star' 'halos_communications'
'halo_update_array_integer_star' 1 ((PROCEDURE UNKNOWN-INTENT
MODULE-PROC DECL UNKNOWN 0 0 SUBROUTINE) (UNKNOWN 0 0 0 UNKNOWN ()) 320
0 (321 322 323) () 0 () () () 0 0)
24 'halo_update_array_integer_block' 'halos_communications'
'halo_update_array_integer_block' 1 ((PROCEDURE UNKNOWN-INTENT
MODULE-PROC DECL UNKNOWN 0 0 SUBROUTINE ALWAYS_EXPLICIT) (UNKNOWN 0 0 0
UNKNOWN ()) 324 0 (325 326) () 0 () () () 0 0)
25 'halo_update_array_integer' 'halos_communications'
'halo_update_array_integer' 1 ((PROCEDURE UNKNOWN-INTENT MODULE-PROC
DECL UNKNOWN 0 0 SUBROUTINE ALWAYS_EXPLICIT) (UNKNOWN 0 0 0 UNKNOWN ())
327 0 (328 329) () 0 () () () 0 0)
23 'halo_update_array_integer_block2' 'halos_communications'
'halo_update_array_integer_block2' 1 ((PROCEDURE UNKNOWN-INTENT
MODULE-PROC DECL UNKNOWN 0 0 SUBROUTINE ALWAYS_EXPLICIT) (UNKNOWN 0 0 0
UNKNOWN ()) 330 0 (331 332) () 0 () () () 0 0)
21 'halo_update_array_real' 'halos_communications'
'halo_update_array_real' 1 ((PROCEDURE UNKNOWN-INTENT MODULE-PROC DECL
UNKNOWN 0 0 SUBROUTINE ALWAYS_EXPLICIT) (UNKNOWN 0 0 0 UNKNOWN ()) 333 0
(334 335) () 0 () () () 0 0)
28 'halo_verifies_scalar' 'halos_communications' 'halo_verifies_scalar'
1 ((PROCEDURE UNKNOWN-INTENT MODULE-PROC DECL UNKNOWN 0 0 FUNCTION) (
LOGICAL 4 0 0 LOGICAL ()) 336 0 (337 338) () 339 () () () 0 0)
30 'halo_verifies_array_integer' 'halos_communications'
'halo_verifies_array_integer' 1 ((PROCEDURE UNKNOWN-INTENT MODULE-PROC
DECL UNKNOWN 0 0 FUNCTION ALWAYS_EXPLICIT) (LOGICAL 4 0 0 LOGICAL ())
340 0 (341 342) () 343 () () () 0 0)
10 'halo_max_scalar' 'halos_communications' 'halo_max_scalar' 1 ((
PROCEDURE UNKNOWN-INTENT MODULE-PROC DECL UNKNOWN 0 0 SUBROUTINE
ALWAYS_EXPLICIT) (UNKNOWN 0 0 0 UNKNOWN ()) 344 0 (345 346) () 0 () () ()
0 0)
29 'halo_verifies_array_real' 'halos_communications'
'halo_verifies_array_real' 1 ((PROCEDURE UNKNOWN-INTENT MODULE-PROC DECL
UNKNOWN 0 0 FUNCTION ALWAYS_EXPLICIT) (LOGICAL 4 0 0 LOGICAL ()) 347 0 (
348 349) () 350 () () () 0 0)
12 'halo_max_array_real' 'halos_communications' 'halo_max_array_real' 1
((PROCEDURE UNKNOWN-INTENT MODULE-PROC DECL UNKNOWN 0 0 SUBROUTINE
ALWAYS_EXPLICIT) (UNKNOWN 0 0 0 UNKNOWN ()) 351 0 (352 353) () 0 () () ()
0 0)
11 'halo_max_scalar_on_halo' 'halos_communications'
'halo_max_scalar_on_halo' 1 ((PROCEDURE UNKNOWN-INTENT MODULE-PROC DECL
UNKNOWN 0 0 SUBROUTINE) (UNKNOWN 0 0 0 UNKNOWN ()) 354 0 (355 356) () 0
() () () 0 0)
27 'halo_verifies_vector_dim' 'halos_communications'
'halo_verifies_vector_dim' 1 ((PROCEDURE UNKNOWN-INTENT MODULE-PROC DECL
UNKNOWN 0 0 FUNCTION) (LOGICAL 4 0 0 LOGICAL ()) 357 0 (358 359 360) ()
361 () () () 0 0)
26 'halo_verifies_vector' 'halos_communications' 'halo_verifies_vector'
1 ((PROCEDURE UNKNOWN-INTENT MODULE-PROC DECL UNKNOWN 0 0 FUNCTION) (
LOGICAL 4 0 0 LOGICAL ()) 362 0 (363 364) () 365 () () () 0 0)
366 'tensor_field' 'fields_data_types' 'tensor_field' 1 ((DERIVED
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 POINTER_COMP) (UNKNOWN 0
0 0 UNKNOWN ()) 0 0 () () 0 ((367 'val' (REAL 8 0 0 REAL ()) (3 0
DEFERRED () () () () () ()) (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0 DIMENSION POINTER) UNKNOWN-ACCESS ()) (368 'wrapped'
(LOGICAL 4 0 0 LOGICAL ()) () (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0) UNKNOWN-ACCESS (CONSTANT (LOGICAL 4 0 0 LOGICAL ())
0 1)) (369 'field_type' (INTEGER 4 0 0 INTEGER ()) () (UNKNOWN-FL
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0) UNKNOWN-ACCESS (
CONSTANT (INTEGER 4 0 0 INTEGER ()) 0 '0')) (370 'name' (CHARACTER 1 0 0
CHARACTER ((CONSTANT (INTEGER 4 0 0 INTEGER ()) 0 '101'))) () (
UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0)
UNKNOWN-ACCESS ()) (371 'dim' (INTEGER 4 0 0 INTEGER ()) (1 0 EXPLICIT (
CONSTANT (INTEGER 4 0 0 INTEGER ()) 0 '1') (CONSTANT (INTEGER 4 0 0
INTEGER ()) 0 '2')) (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN
UNKNOWN 0 0 DIMENSION) UNKNOWN-ACCESS ()) (372 'option_path' (CHARACTER
1 0 0 CHARACTER ((CONSTANT (INTEGER 4 0 0 INTEGER ()) 0 '8192'))) () (
UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0)
UNKNOWN-ACCESS (CONSTANT (CHARACTER 1 0 0 CHARACTER (())) 0 8192
'/uninitialised_path/                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            '))
(373 'mesh' (DERIVED 105 0 0 DERIVED ()) () (UNKNOWN-FL UNKNOWN-INTENT
UNKNOWN-PROC UNKNOWN UNKNOWN 0 0) UNKNOWN-ACCESS (STRUCTURE (DERIVED 105
0 0 DERIVED ()) 0 ((() ()) ((CONSTANT (LOGICAL 4 0 0 LOGICAL ()) 0 1) ())
((STRUCTURE (DERIVED 106 0 0 DERIVED ()) 0 ((() ()) (() ()) (() ()) (()
()) ((NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0) ()) ((NULL (UNKNOWN 0 0 0
UNKNOWN ()) 0) ()) ((NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0) ()) ((NULL (
UNKNOWN 0 0 0 UNKNOWN ()) 0) ()) ((NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0) ())
((NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0) ()) ((NULL (UNKNOWN 0 0 0 UNKNOWN ())
0) ()) ((STRUCTURE (DERIVED 107 0 0 DERIVED ()) 0 ((() ()) (() ()) (() ())
(() ()) ((NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0) ()) ((NULL (UNKNOWN 0 0 0
UNKNOWN ()) 0) ()) (() ()) ((NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0) ()) (()
())) ()) ()) ((NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0) ()) ((NULL (UNKNOWN 0
0 0 UNKNOWN ()) 0) ()) ((NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0) ()) ((NULL (
UNKNOWN 0 0 0 UNKNOWN ()) 0) ()) (() ())) ()) ()) (() ()) (() ()) (() ())
((CONSTANT (CHARACTER 1 0 0 CHARACTER (())) 0 8192
'/uninitialised_path/                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            ')
()) ((CONSTANT (INTEGER 4 0 0 INTEGER ()) 0 '0') ()) ((NULL (UNKNOWN 0 0
0 UNKNOWN ()) 0) ()) ((NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0) ()) ((NULL (
UNKNOWN 0 0 0 UNKNOWN ()) 0) ()) ((NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0) ())
((NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0) ()) ((NULL (UNKNOWN 0 0 0 UNKNOWN ())
0) ()) ((NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0) ()) ((NULL (UNKNOWN 0 0 0
UNKNOWN ()) 0) ()) ((NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0) ()) ((NULL (
UNKNOWN 0 0 0 UNKNOWN ()) 0) ()) ((CONSTANT (LOGICAL 4 0 0 LOGICAL ()) 0
0) ())) ())) (374 'refcount' (DERIVED 77 0 0 DERIVED ()) () (UNKNOWN-FL
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 POINTER) UNKNOWN-ACCESS
(NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0)) (375 'aliased' (LOGICAL 4 0 0
LOGICAL ()) () (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0
0) UNKNOWN-ACCESS (CONSTANT (LOGICAL 4 0 0 LOGICAL ()) 0 0))) PUBLIC (()
() () ()) () 0 0 25110185)
5 'get_universal_numbering_multiple_components' 'halos_numbering'
'get_universal_numbering_multiple_components' 1 ((PROCEDURE
UNKNOWN-INTENT MODULE-PROC DECL UNKNOWN 0 0 SUBROUTINE ALWAYS_EXPLICIT)
(UNKNOWN 0 0 0 UNKNOWN ()) 376 0 (377 378) () 0 () () () 0 0)
44 'halo_universal_number_vector' 'halos_numbering'
'halo_universal_number_vector' 1 ((PROCEDURE UNKNOWN-INTENT MODULE-PROC
DECL UNKNOWN 0 0 DIMENSION FUNCTION ALWAYS_EXPLICIT) (INTEGER 4 0 0
INTEGER ()) 379 0 (380 381) (1 0 EXPLICIT (CONSTANT (INTEGER 4 0 0
INTEGER ()) 0 '1') (FUNCTION (INTEGER 4 0 0 INTEGER ()) 0 382 (('' (
VARIABLE (INTEGER 4 0 0 INTEGER ()) 1 381 ((ARRAY (FULL 0))))) ('' ()) (
'' ())) '' 0 'size')) 383 () () () 0 0)
384 'halo_proc_count' 'halos_base' 'halo_proc_count' 1 ((PROCEDURE
UNKNOWN-INTENT MODULE-PROC DECL UNKNOWN 0 0 FUNCTION PURE) (INTEGER 4 0
0 INTEGER ()) 385 0 (386) () 384 () () () 0 0)
32 'node_count_halo' 'halos_base' 'node_count_halo' 1 ((PROCEDURE
UNKNOWN-INTENT MODULE-PROC DECL UNKNOWN 0 0 FUNCTION PURE) (INTEGER 4 0
0 INTEGER ()) 387 0 (388) () 389 () () () 0 0)
39 'nodes_owned_halo' 'halos_ownership' 'nodes_owned_halo' 1 ((
PROCEDURE UNKNOWN-INTENT MODULE-PROC DECL UNKNOWN 0 0 DIMENSION FUNCTION
ALWAYS_EXPLICIT) (LOGICAL 4 0 0 LOGICAL ()) 390 0 (391 392) (1 0
EXPLICIT (CONSTANT (INTEGER 4 0 0 INTEGER ()) 0 '1') (FUNCTION (INTEGER
4 0 0 INTEGER ()) 0 393 (('' (VARIABLE (INTEGER 4 0 0 INTEGER ()) 1 392
((ARRAY (FULL 0))))) ('' ()) ('' ())) '' 0 'size')) 394 () () () 0 0)
42 'node_owned_halo' 'halos_ownership' 'node_owned_halo' 1 ((PROCEDURE
UNKNOWN-INTENT MODULE-PROC DECL UNKNOWN 0 0 FUNCTION) (LOGICAL 4 0 0
LOGICAL ()) 395 0 (396 397) () 42 () () () 0 0)
38 'read_halos_mesh' 'halos_registration' 'read_halos_mesh' 1 ((
PROCEDURE UNKNOWN-INTENT MODULE-PROC DECL UNKNOWN 0 0 SUBROUTINE
ALWAYS_EXPLICIT) (UNKNOWN 0 0 0 UNKNOWN ()) 398 0 (399 400 401) () 0 ()
() () 0 0)
46 'serial_storage_halo_multiple' 'halos_base'
'serial_storage_halo_multiple' 1 ((PROCEDURE UNKNOWN-INTENT MODULE-PROC
DECL UNKNOWN 0 0 DIMENSION FUNCTION ALWAYS_EXPLICIT) (LOGICAL 4 0 0
LOGICAL ()) 402 0 (403) (1 0 EXPLICIT (CONSTANT (INTEGER 4 0 0 INTEGER ())
0 '1') (FUNCTION (INTEGER 4 0 0 INTEGER ()) 0 404 (('' (VARIABLE (
DERIVED 131 0 0 DERIVED ()) 1 403 ((ARRAY (FULL 0))))) ('' ()) ('' ())) ''
0 'size')) 405 () () () 0 0)
37 'read_halos_positions' 'halos_registration' 'read_halos_positions' 1
((PROCEDURE UNKNOWN-INTENT MODULE-PROC DECL UNKNOWN 0 0 SUBROUTINE
ALWAYS_EXPLICIT) (UNKNOWN 0 0 0 UNKNOWN ()) 406 0 (407 408 409) () 0 ()
() () 0 0)
9 'halo_communicator_halo' 'halos_base' 'halo_communicator_halo' 1 ((
PROCEDURE UNKNOWN-INTENT MODULE-PROC DECL UNKNOWN 0 0 FUNCTION PURE) (
INTEGER 4 0 0 INTEGER ()) 410 0 (411) () 412 () () () 0 0)
47 'serial_storage_halo_single' 'halos_base' 'serial_storage_halo_single'
1 ((PROCEDURE UNKNOWN-INTENT MODULE-PROC DECL UNKNOWN 0 0 FUNCTION) (
LOGICAL 4 0 0 LOGICAL ()) 413 0 (414) () 415 () () () 0 0)
40 'reorder_halo_halo' 'halos_repair' 'reorder_halo_halo' 1 ((PROCEDURE
UNKNOWN-INTENT MODULE-PROC DECL UNKNOWN 0 0 SUBROUTINE) (UNKNOWN 0 0 0
UNKNOWN ()) 416 0 (417 418) () 0 () () () 0 0)
41 'reorder_halo_vector' 'halos_repair' 'reorder_halo_vector' 1 ((
PROCEDURE UNKNOWN-INTENT MODULE-PROC DECL UNKNOWN 0 0 SUBROUTINE) (
UNKNOWN 0 0 0 UNKNOWN ()) 419 0 (420 421) () 0 () () () 0 0)
48 'zero_halo' 'halos_base' 'zero_halo' 1 ((PROCEDURE UNKNOWN-INTENT
MODULE-PROC DECL UNKNOWN 0 0 SUBROUTINE) (UNKNOWN 0 0 0 UNKNOWN ()) 422
0 (423) () 0 () () () 0 0)
7 'allocate_halo_halo' 'halos_allocates' 'allocate_halo_halo' 1 ((
PROCEDURE UNKNOWN-INTENT MODULE-PROC DECL UNKNOWN 0 0 SUBROUTINE) (
UNKNOWN 0 0 0 UNKNOWN ()) 424 0 (425 426) () 0 () () () 0 0)
427 'create_global_to_universal_numbering' 'halos_numbering'
'create_global_to_universal_numbering' 1 ((PROCEDURE UNKNOWN-INTENT
MODULE-PROC DECL UNKNOWN 0 0 SUBROUTINE ALWAYS_EXPLICIT) (UNKNOWN 0 0 0
UNKNOWN ()) 428 0 (429 430) () 0 () () () 0 0)
431 'create_ownership' 'halos_ownership' 'create_ownership' 1 ((
PROCEDURE UNKNOWN-INTENT MODULE-PROC DECL UNKNOWN 0 0 SUBROUTINE) (
UNKNOWN 0 0 0 UNKNOWN ()) 432 0 (433) () 0 () () () 0 0)
434 'deallocate_ownership_cache' 'halos_allocates'
'deallocate_ownership_cache' 1 ((PROCEDURE UNKNOWN-INTENT MODULE-PROC
DECL UNKNOWN 0 0 SUBROUTINE IMPLICIT_PURE) (UNKNOWN 0 0 0 UNKNOWN ())
435 0 (436) () 0 () () () 0 0)
437 'deallocate_universal_numbering_cache' 'halos_allocates'
'deallocate_universal_numbering_cache' 1 ((PROCEDURE UNKNOWN-INTENT
MODULE-PROC DECL UNKNOWN 0 0 SUBROUTINE IMPLICIT_PURE) (UNKNOWN 0 0 0
UNKNOWN ()) 438 0 (439) () 0 () () () 0 0)
440 'derive_element_halo_from_node_halo' 'halos_derivation'
'derive_element_halo_from_node_halo' 1 ((PROCEDURE UNKNOWN-INTENT
MODULE-PROC DECL UNKNOWN 0 0 SUBROUTINE ALWAYS_EXPLICIT) (UNKNOWN 0 0 0
UNKNOWN ()) 441 0 (442 443 444) () 0 () () () 0 0)
445 'derive_maximal_surface_element_halo' 'halos_derivation'
'derive_maximal_surface_element_halo' 1 ((PROCEDURE UNKNOWN-INTENT
MODULE-PROC DECL UNKNOWN 0 0 FUNCTION ALWAYS_EXPLICIT) (DERIVED 131 0 0
DERIVED ()) 446 0 (447 448 449 450) () 451 () () () 0 0)
452 'derive_nonperiodic_halos_from_periodic_halos' 'halos_derivation'
'derive_nonperiodic_halos_from_periodic_halos' 1 ((PROCEDURE
UNKNOWN-INTENT MODULE-PROC DECL UNKNOWN 0 0 SUBROUTINE) (UNKNOWN 0 0 0
UNKNOWN ()) 453 0 (454 455 456) () 0 () () () 0 0)
457 'derive_sub_halo' 'halos_derivation' 'derive_sub_halo' 1 ((
PROCEDURE UNKNOWN-INTENT MODULE-PROC DECL UNKNOWN 0 0 FUNCTION
ALWAYS_EXPLICIT) (DERIVED 131 0 0 DERIVED ()) 458 0 (459 460) () 461 ()
() () 0 0)
462 'ele_owner' 'halos_derivation' 'ele_owner' 1 ((PROCEDURE
UNKNOWN-INTENT MODULE-PROC DECL UNKNOWN 0 0 FUNCTION) (INTEGER 4 0 0
INTEGER ()) 463 0 (464 465 466) () 462 () () () 0 0)
467 'ewrite_universal_numbers' 'halos_numbering'
'ewrite_universal_numbers' 1 ((PROCEDURE UNKNOWN-INTENT MODULE-PROC DECL
UNKNOWN 0 0 SUBROUTINE) (UNKNOWN 0 0 0 UNKNOWN ()) 468 0 (469 470) () 0
() () () 0 0)
471 'extract_all_halo_receives' 'halos_base' 'extract_all_halo_receives'
1 ((PROCEDURE UNKNOWN-INTENT MODULE-PROC DECL UNKNOWN 0 0 SUBROUTINE
ALWAYS_EXPLICIT) (UNKNOWN 0 0 0 UNKNOWN ()) 472 0 (473 474 475 476) () 0
() () () 0 0)
477 'extract_all_halo_sends' 'halos_base' 'extract_all_halo_sends' 1 ((
PROCEDURE UNKNOWN-INTENT MODULE-PROC DECL UNKNOWN 0 0 SUBROUTINE
ALWAYS_EXPLICIT) (UNKNOWN 0 0 0 UNKNOWN ()) 478 0 (479 480 481 482) () 0
() () () 0 0)
483 'extract_raw_halo_data' 'halos_registration' 'extract_raw_halo_data'
1 ((PROCEDURE UNKNOWN-INTENT MODULE-PROC DECL UNKNOWN 0 0 SUBROUTINE
ALWAYS_EXPLICIT) (UNKNOWN 0 0 0 UNKNOWN ()) 484 0 (485 486 487 488 489
490) () 0 () () () 0 0)
491 'form_halo_from_raw_data' 'halos_registration'
'form_halo_from_raw_data' 1 ((PROCEDURE UNKNOWN-INTENT MODULE-PROC DECL
UNKNOWN 0 0 SUBROUTINE ALWAYS_EXPLICIT) (UNKNOWN 0 0 0 UNKNOWN ()) 492 0
(493 494 495 496 497 498 499 500 501) () 0 () () () 0 0)
502 'generate_substate_halos' 'halos_registration'
'generate_substate_halos' 1 ((PROCEDURE UNKNOWN-INTENT MODULE-PROC DECL
UNKNOWN 0 0 SUBROUTINE ALWAYS_EXPLICIT) (UNKNOWN 0 0 0 UNKNOWN ()) 503 0
(504 505 506 507) () 0 () () () 0 0)
508 'get_node_owners' 'halos_ownership' 'get_node_owners' 1 ((PROCEDURE
UNKNOWN-INTENT MODULE-PROC DECL UNKNOWN 0 0 SUBROUTINE ALWAYS_EXPLICIT)
(UNKNOWN 0 0 0 UNKNOWN ()) 509 0 (510 511) () 0 () () () 0 0)
512 'get_owned_nodes' 'halos_ownership' 'get_owned_nodes' 1 ((PROCEDURE
UNKNOWN-INTENT MODULE-PROC DECL UNKNOWN 0 0 SUBROUTINE ALWAYS_EXPLICIT)
(UNKNOWN 0 0 0 UNKNOWN ()) 513 0 (514 515) () 0 () () () 0 0)
6 'get_universal_numbering' 'halos_numbering' 'get_universal_numbering'
1 ((PROCEDURE UNKNOWN-INTENT MODULE-PROC DECL UNKNOWN 0 0 SUBROUTINE
GENERIC) (UNKNOWN 0 0 0 UNKNOWN ()) 516 0 (517 518) () 0 () () () 0 0)
519 'get_universal_numbering_inverse' 'halos_numbering'
'get_universal_numbering_inverse' 1 ((PROCEDURE UNKNOWN-INTENT
MODULE-PROC DECL UNKNOWN 0 0 SUBROUTINE) (UNKNOWN 0 0 0 UNKNOWN ()) 520
0 (521 522) () 0 () () () 0 0)
523 'halo_all_receives_count' 'halos_base' 'halo_all_receives_count' 1 (
(PROCEDURE UNKNOWN-INTENT MODULE-PROC DECL UNKNOWN 0 0 FUNCTION PURE) (
INTEGER 4 0 0 INTEGER ()) 524 0 (525) () 523 () () () 0 0)
526 'halo_all_sends_count' 'halos_base' 'halo_all_sends_count' 1 ((
PROCEDURE UNKNOWN-INTENT MODULE-PROC DECL UNKNOWN 0 0 FUNCTION PURE) (
INTEGER 4 0 0 INTEGER ()) 527 0 (528) () 526 () () () 0 0)
529 'halo_all_unique_receives_count' 'halos_base'
'halo_all_unique_receives_count' 1 ((PROCEDURE UNKNOWN-INTENT
MODULE-PROC DECL UNKNOWN 0 0 FUNCTION) (INTEGER 4 0 0 INTEGER ()) 530 0
(531) () 529 () () () 0 0)
532 'halo_data_type' 'halos_base' 'halo_data_type' 1 ((PROCEDURE
UNKNOWN-INTENT MODULE-PROC DECL UNKNOWN 0 0 FUNCTION PURE) (INTEGER 4 0
0 INTEGER ()) 533 0 (534) () 532 () () () 0 0)
535 'halo_name' 'halos_base' 'halo_name' 1 ((PROCEDURE UNKNOWN-INTENT
MODULE-PROC DECL UNKNOWN 0 0 FUNCTION PURE) (CHARACTER 1 0 0 CHARACTER (
(FUNCTION (INTEGER 4 0 0 INTEGER ()) 0 536 (('' (VARIABLE (CHARACTER 1 0
0 CHARACTER ((CONSTANT (INTEGER 4 0 0 INTEGER ()) 0 '101'))) 0 537 ((
COMPONENT 131 538 'name')))) ('' ())) '__len_trim1' 0 'lnblnk'))) 539 0
(537) () 535 () () () 0 0)
540 'halo_node_owner' 'halos_ownership' 'halo_node_owner' 1 ((PROCEDURE
UNKNOWN-INTENT MODULE-PROC DECL UNKNOWN 0 0 FUNCTION ALWAYS_EXPLICIT) (
INTEGER 4 0 0 INTEGER ()) 541 0 (542 543 544) () 545 () () () 0 0)
546 'halo_node_owners' 'halos_ownership' 'halo_node_owners' 1 ((
PROCEDURE UNKNOWN-INTENT MODULE-PROC DECL UNKNOWN 0 0 DIMENSION FUNCTION
ALWAYS_EXPLICIT) (INTEGER 4 0 0 INTEGER ()) 547 0 (548 549 550) (1 0
EXPLICIT (CONSTANT (INTEGER 4 0 0 INTEGER ()) 0 '1') (FUNCTION (INTEGER
4 0 0 INTEGER ()) 0 551 (('' (VARIABLE (INTEGER 4 0 0 INTEGER ()) 1 549
((ARRAY (FULL 0))))) ('' ()) ('' ())) '' 0 'size')) 552 () () () 0 0)
553 'halo_nowned_nodes' 'halos_base' 'halo_nowned_nodes' 1 ((PROCEDURE
UNKNOWN-INTENT MODULE-PROC DECL UNKNOWN 0 0 FUNCTION PURE) (INTEGER 4 0
0 INTEGER ()) 554 0 (555) () 553 () () () 0 0)
556 'halo_order_general' 'halo_data_types' 'halo_order_general' 1 ((
PARAMETER UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN IMPLICIT-SAVE 0 0) (
INTEGER 4 0 0 INTEGER ()) 0 0 () (CONSTANT (INTEGER 4 0 0 INTEGER ()) 0
'1') () 0 () () () 0 0)
557 'halo_order_trailing_receives' 'halo_data_types'
'halo_order_trailing_receives' 1 ((PARAMETER UNKNOWN-INTENT UNKNOWN-PROC
UNKNOWN IMPLICIT-SAVE 0 0) (INTEGER 4 0 0 INTEGER ()) 0 0 () (CONSTANT (
INTEGER 4 0 0 INTEGER ()) 0 '2') () 0 () () () 0 0)
558 'halo_ordering_scheme' 'halos_base' 'halo_ordering_scheme' 1 ((
PROCEDURE UNKNOWN-INTENT MODULE-PROC DECL UNKNOWN 0 0 FUNCTION PURE) (
INTEGER 4 0 0 INTEGER ()) 559 0 (560) () 558 () () () 0 0)
561 'halo_pointer' 'halo_data_types' 'halo_pointer' 1 ((DERIVED
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 POINTER_COMP) (UNKNOWN 0
0 0 UNKNOWN ()) 0 0 () () 0 ((562 'ptr' (DERIVED 131 0 0 DERIVED ()) ()
(UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 POINTER)
UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0))) PUBLIC (() () () ())
() 0 0 72023282)
563 'halo_receive' 'halos_base' 'halo_receive' 1 ((PROCEDURE
UNKNOWN-INTENT MODULE-PROC DECL UNKNOWN 0 0 FUNCTION) (INTEGER 4 0 0
INTEGER ()) 564 0 (565 566 567) () 563 () () () 0 0)
568 'halo_receive_count' 'halos_base' 'halo_receive_count' 1 ((
PROCEDURE UNKNOWN-INTENT MODULE-PROC DECL UNKNOWN 0 0 FUNCTION PURE) (
INTEGER 4 0 0 INTEGER ()) 569 0 (570 571) () 568 () () () 0 0)
572 'halo_receive_counts' 'halos_base' 'halo_receive_counts' 1 ((
PROCEDURE UNKNOWN-INTENT MODULE-PROC DECL UNKNOWN 0 0 SUBROUTINE
IMPLICIT_PURE) (UNKNOWN 0 0 0 UNKNOWN ()) 573 0 (574 575) () 0 () () ()
0 0)
576 'halo_receives' 'halos_base' 'halo_receives' 1 ((PROCEDURE
UNKNOWN-INTENT MODULE-PROC DECL UNKNOWN 0 0 DIMENSION POINTER FUNCTION
ALWAYS_EXPLICIT) (INTEGER 4 0 0 INTEGER ()) 577 0 (578 579) (1 0
DEFERRED () ()) 576 () () () 0 0)
580 'halo_send' 'halos_base' 'halo_send' 1 ((PROCEDURE UNKNOWN-INTENT
MODULE-PROC DECL UNKNOWN 0 0 FUNCTION) (INTEGER 4 0 0 INTEGER ()) 581 0
(582 583 584) () 580 () () () 0 0)
585 'halo_send_count' 'halos_base' 'halo_send_count' 1 ((PROCEDURE
UNKNOWN-INTENT MODULE-PROC DECL UNKNOWN 0 0 FUNCTION PURE) (INTEGER 4 0
0 INTEGER ()) 586 0 (587 588) () 585 () () () 0 0)
589 'halo_send_counts' 'halos_base' 'halo_send_counts' 1 ((PROCEDURE
UNKNOWN-INTENT MODULE-PROC DECL UNKNOWN 0 0 SUBROUTINE IMPLICIT_PURE) (
UNKNOWN 0 0 0 UNKNOWN ()) 590 0 (591 592) () 0 () () () 0 0)
593 'halo_sends' 'halos_base' 'halo_sends' 1 ((PROCEDURE UNKNOWN-INTENT
MODULE-PROC DECL UNKNOWN 0 0 DIMENSION POINTER FUNCTION ALWAYS_EXPLICIT)
(INTEGER 4 0 0 INTEGER ()) 594 0 (595 596) (1 0 DEFERRED () ()) 593 () ()
() 0 0)
131 'halo_type' 'halo_data_types' 'halo_type' 1 ((DERIVED UNKNOWN-INTENT
UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 POINTER_COMP) (UNKNOWN 0 0 0 UNKNOWN ())
0 0 () () 0 ((538 'name' (CHARACTER 1 0 0 CHARACTER ((CONSTANT (INTEGER
4 0 0 INTEGER ()) 0 '101'))) () (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0) UNKNOWN-ACCESS ()) (597 'refcount' (DERIVED 77 0 0
DERIVED ()) () (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0
0 POINTER) UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0)) (598
'data_type' (INTEGER 4 0 0 INTEGER ()) () (UNKNOWN-FL UNKNOWN-INTENT
UNKNOWN-PROC UNKNOWN UNKNOWN 0 0) UNKNOWN-ACCESS (CONSTANT (INTEGER 4 0
0 INTEGER ()) 0 '0')) (599 'ordering_scheme' (INTEGER 4 0 0 INTEGER ())
() (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0)
UNKNOWN-ACCESS (CONSTANT (INTEGER 4 0 0 INTEGER ()) 0 '0')) (600
'communicator' (INTEGER 4 0 0 INTEGER ()) () (UNKNOWN-FL UNKNOWN-INTENT
UNKNOWN-PROC UNKNOWN UNKNOWN 0 0) UNKNOWN-ACCESS ()) (601 'nprocs' (
INTEGER 4 0 0 INTEGER ()) () (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0) UNKNOWN-ACCESS (CONSTANT (INTEGER 4 0 0 INTEGER ())
0 '0')) (602 'sends' (DERIVED 85 0 0 DERIVED ()) (1 0 DEFERRED () ()) (
UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 DIMENSION
POINTER) UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0)) (603
'receives' (DERIVED 85 0 0 DERIVED ()) (1 0 DEFERRED () ()) (UNKNOWN-FL
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 DIMENSION POINTER)
UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0)) (604 'nowned_nodes'
(INTEGER 4 0 0 INTEGER ()) () (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0) UNKNOWN-ACCESS (CONSTANT (INTEGER 4 0 0 INTEGER ())
0 '-1')) (605 'owners' (INTEGER 4 0 0 INTEGER ()) (1 0 DEFERRED () ()) (
UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 DIMENSION
POINTER) UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0)) (606
'unn_count' (INTEGER 4 0 0 INTEGER ()) () (UNKNOWN-FL UNKNOWN-INTENT
UNKNOWN-PROC UNKNOWN UNKNOWN 0 0) UNKNOWN-ACCESS (CONSTANT (INTEGER 4 0
0 INTEGER ()) 0 '-1')) (607 'owned_nodes_unn_base' (INTEGER 4 0 0
INTEGER ()) (1 0 DEFERRED () ()) (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0 DIMENSION POINTER) UNKNOWN-ACCESS (NULL (UNKNOWN 0 0
0 UNKNOWN ()) 0)) (608 'my_owned_nodes_unn_base' (INTEGER 4 0 0 INTEGER
()) () (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0)
UNKNOWN-ACCESS (CONSTANT (INTEGER 4 0 0 INTEGER ()) 0 '-1')) (609
'receives_gnn_to_unn' (INTEGER 4 0 0 INTEGER ()) (1 0 DEFERRED () ()) (
UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 DIMENSION
POINTER) UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0)) (610
'gnn_to_unn' (INTEGER 4 0 0 INTEGER ()) (1 0 DEFERRED () ()) (
UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 DIMENSION
POINTER) UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0))) PUBLIC (()
() () ()) () 0 0 63994053)
611 'halo_type_cg_node' 'halo_data_types' 'halo_type_cg_node' 1 ((
PARAMETER UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN IMPLICIT-SAVE 0 0) (
INTEGER 4 0 0 INTEGER ()) 0 0 () (CONSTANT (INTEGER 4 0 0 INTEGER ()) 0
'1') () 0 () () () 0 0)
612 'halo_type_dg_node' 'halo_data_types' 'halo_type_dg_node' 1 ((
PARAMETER UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN IMPLICIT-SAVE 0 0) (
INTEGER 4 0 0 INTEGER ()) 0 0 () (CONSTANT (INTEGER 4 0 0 INTEGER ()) 0
'2') () 0 () () () 0 0)
613 'halo_type_element' 'halo_data_types' 'halo_type_element' 1 ((
PARAMETER UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN IMPLICIT-SAVE 0 0) (
INTEGER 4 0 0 INTEGER ()) 0 0 () (CONSTANT (INTEGER 4 0 0 INTEGER ()) 0
'3') () 0 () () () 0 0)
614 'halo_unique_receive_count' 'halos_base' 'halo_unique_receive_count'
1 ((PROCEDURE UNKNOWN-INTENT MODULE-PROC DECL UNKNOWN 0 0 FUNCTION) (
INTEGER 4 0 0 INTEGER ()) 615 0 (616 617) () 614 () () () 0 0)
618 'halo_universal_node_owners' 'halos_ownership'
'halo_universal_node_owners' 1 ((PROCEDURE UNKNOWN-INTENT MODULE-PROC
DECL UNKNOWN 0 0 DIMENSION FUNCTION ALWAYS_EXPLICIT) (INTEGER 4 0 0
INTEGER ()) 619 0 (620 621) (1 0 EXPLICIT (CONSTANT (INTEGER 4 0 0
INTEGER ()) 0 '1') (FUNCTION (INTEGER 4 0 0 INTEGER ()) 0 622 (('' (
VARIABLE (INTEGER 4 0 0 INTEGER ()) 1 621 ((ARRAY (FULL 0))))) ('' ()) (
'' ())) '' 0 'size')) 623 () () () 0 0)
45 'halo_universal_number' 'halos_numbering' 'halo_universal_number' 1 (
(PROCEDURE UNKNOWN-INTENT MODULE-PROC DECL UNKNOWN 0 0 FUNCTION GENERIC)
(INTEGER 4 0 0 INTEGER ()) 624 0 (625 626) () 627 () () () 0 0)
628 'halo_universal_numbers' 'halos_numbering' 'halo_universal_numbers'
1 ((PROCEDURE UNKNOWN-INTENT MODULE-PROC DECL UNKNOWN 0 0 DIMENSION
FUNCTION ALWAYS_EXPLICIT) (INTEGER 4 0 0 INTEGER ()) 629 0 (630 631) (1
0 EXPLICIT (CONSTANT (INTEGER 4 0 0 INTEGER ()) 0 '1') (FUNCTION (
INTEGER 4 0 0 INTEGER ()) 0 632 (('' (VARIABLE (INTEGER 4 0 0 INTEGER ())
1 631 ((ARRAY (FULL 0))))) ('' ()) ('' ())) '' 0 'size')) 633 () () () 0
0)
634 'halo_valid_for_communication' 'halos_debug'
'halo_valid_for_communication' 1 ((PROCEDURE UNKNOWN-INTENT MODULE-PROC
DECL UNKNOWN 0 0 FUNCTION) (LOGICAL 4 0 0 LOGICAL ()) 635 0 (636) () 637
() () () 0 0)
638 'halos' 'halos' 'halos' 1 ((MODULE UNKNOWN-INTENT UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0) (UNKNOWN 0 0 0 UNKNOWN ()) 0 0 () () 0 () () () 0 0)
639 'has_global_to_universal_numbering' 'halos_numbering'
'has_global_to_universal_numbering' 1 ((PROCEDURE UNKNOWN-INTENT
MODULE-PROC DECL UNKNOWN 0 0 FUNCTION) (LOGICAL 4 0 0 LOGICAL ()) 640 0
(641) () 642 () () () 0 0)
643 'has_nowned_nodes' 'halos_base' 'has_nowned_nodes' 1 ((PROCEDURE
UNKNOWN-INTENT MODULE-PROC DECL UNKNOWN 0 0 FUNCTION PURE) (LOGICAL 4 0
0 LOGICAL ()) 644 0 (645) () 643 () () () 0 0)
646 'has_ownership' 'halos_ownership' 'has_ownership' 1 ((PROCEDURE
UNKNOWN-INTENT MODULE-PROC DECL UNKNOWN 0 0 FUNCTION PURE) (LOGICAL 4 0
0 LOGICAL ()) 647 0 (648) () 646 () () () 0 0)
649 'invert_comms_sizes' 'halos_derivation' 'invert_comms_sizes' 1 ((
PROCEDURE UNKNOWN-INTENT MODULE-PROC DECL UNKNOWN 0 0 DIMENSION FUNCTION
ALWAYS_EXPLICIT) (INTEGER 4 0 0 INTEGER ()) 650 0 (651 652) (1 0
EXPLICIT (CONSTANT (INTEGER 4 0 0 INTEGER ()) 0 '1') (FUNCTION (INTEGER
4 0 0 INTEGER ()) 0 653 (('' (VARIABLE (INTEGER 4 0 0 INTEGER ()) 1 651
((ARRAY (FULL 0))))) ('' ()) ('' ())) '' 0 'size')) 654 () () () 0 0)
655 'max_halo_node' 'halos_base' 'max_halo_node' 1 ((PROCEDURE
UNKNOWN-INTENT MODULE-PROC DECL UNKNOWN 0 0 FUNCTION) (INTEGER 4 0 0
INTEGER ()) 656 0 (657) () 658 () () () 0 0)
659 'max_halo_receive_node' 'halos_base' 'max_halo_receive_node' 1 ((
PROCEDURE UNKNOWN-INTENT MODULE-PROC DECL UNKNOWN 0 0 FUNCTION) (
INTEGER 4 0 0 INTEGER ()) 660 0 (661) () 662 () () () 0 0)
663 'max_halo_send_node' 'halos_base' 'max_halo_send_node' 1 ((
PROCEDURE UNKNOWN-INTENT MODULE-PROC DECL UNKNOWN 0 0 FUNCTION) (
INTEGER 4 0 0 INTEGER ()) 664 0 (665) () 666 () () () 0 0)
667 'min_halo_node' 'halos_base' 'min_halo_node' 1 ((PROCEDURE
UNKNOWN-INTENT MODULE-PROC DECL UNKNOWN 0 0 FUNCTION) (INTEGER 4 0 0
INTEGER ()) 668 0 (669) () 670 () () () 0 0)
671 'min_halo_receive_node' 'halos_base' 'min_halo_receive_node' 1 ((
PROCEDURE UNKNOWN-INTENT MODULE-PROC DECL UNKNOWN 0 0 FUNCTION) (
INTEGER 4 0 0 INTEGER ()) 672 0 (673) () 674 () () () 0 0)
675 'min_halo_send_node' 'halos_base' 'min_halo_send_node' 1 ((
PROCEDURE UNKNOWN-INTENT MODULE-PROC DECL UNKNOWN 0 0 FUNCTION) (
INTEGER 4 0 0 INTEGER ()) 676 0 (677) () 678 () () () 0 0)
679 'print_halo' 'halos_debug' 'print_halo' 1 ((PROCEDURE UNKNOWN-INTENT
MODULE-PROC DECL UNKNOWN 0 0 SUBROUTINE) (UNKNOWN 0 0 0 UNKNOWN ()) 680
0 (681 682) () 0 () () () 0 0)
683 'reorder_element_halo' 'halos_repair' 'reorder_element_halo' 1 ((
PROCEDURE UNKNOWN-INTENT MODULE-PROC DECL UNKNOWN 0 0 SUBROUTINE) (
UNKNOWN 0 0 0 UNKNOWN ()) 684 0 (685 686 687) () 0 () () () 0 0)
688 'reorder_halo_from_element_halo' 'halos_repair'
'reorder_halo_from_element_halo' 1 ((PROCEDURE UNKNOWN-INTENT
MODULE-PROC DECL UNKNOWN 0 0 SUBROUTINE) (UNKNOWN 0 0 0 UNKNOWN ()) 689
0 (690 691 692) () 0 () () () 0 0)
693 'reorder_halo_receives' 'halos_repair' 'reorder_halo_receives' 1 ((
PROCEDURE UNKNOWN-INTENT MODULE-PROC DECL UNKNOWN 0 0 SUBROUTINE) (
UNKNOWN 0 0 0 UNKNOWN ()) 694 0 (695 696) () 0 () () () 0 0)
697 'reorder_l1_from_l2_halo' 'halos_repair' 'reorder_l1_from_l2_halo' 1
((PROCEDURE UNKNOWN-INTENT MODULE-PROC DECL UNKNOWN 0 0 SUBROUTINE
ALWAYS_EXPLICIT) (UNKNOWN 0 0 0 UNKNOWN ()) 698 0 (699 700 701) () 0 ()
() () 0 0)
702 'set_all_halo_receives' 'halos_base' 'set_all_halo_receives' 1 ((
PROCEDURE UNKNOWN-INTENT MODULE-PROC DECL UNKNOWN 0 0 SUBROUTINE
ALWAYS_EXPLICIT) (UNKNOWN 0 0 0 UNKNOWN ()) 703 0 (704 705) () 0 () () ()
0 0)
706 'set_all_halo_sends' 'halos_base' 'set_all_halo_sends' 1 ((
PROCEDURE UNKNOWN-INTENT MODULE-PROC DECL UNKNOWN 0 0 SUBROUTINE
ALWAYS_EXPLICIT) (UNKNOWN 0 0 0 UNKNOWN ()) 707 0 (708 709) () 0 () () ()
0 0)
710 'set_halo_communicator' 'halos_base' 'set_halo_communicator' 1 ((
PROCEDURE UNKNOWN-INTENT MODULE-PROC DECL UNKNOWN 0 0 SUBROUTINE) (
UNKNOWN 0 0 0 UNKNOWN ()) 711 0 (712 713) () 0 () () () 0 0)
714 'set_halo_data_type' 'halos_base' 'set_halo_data_type' 1 ((
PROCEDURE UNKNOWN-INTENT MODULE-PROC DECL UNKNOWN 0 0 SUBROUTINE) (
UNKNOWN 0 0 0 UNKNOWN ()) 715 0 (716 717) () 0 () () () 0 0)
718 'set_halo_name' 'halos_base' 'set_halo_name' 1 ((PROCEDURE
UNKNOWN-INTENT MODULE-PROC DECL UNKNOWN 0 0 SUBROUTINE IMPLICIT_PURE) (
UNKNOWN 0 0 0 UNKNOWN ()) 719 0 (720 721) () 0 () () () 0 0)
722 'set_halo_nowned_nodes' 'halos_base' 'set_halo_nowned_nodes' 1 ((
PROCEDURE UNKNOWN-INTENT MODULE-PROC DECL UNKNOWN 0 0 SUBROUTINE) (
UNKNOWN 0 0 0 UNKNOWN ()) 723 0 (724 725) () 0 () () () 0 0)
726 'set_halo_ordering_scheme' 'halos_base' 'set_halo_ordering_scheme' 1
((PROCEDURE UNKNOWN-INTENT MODULE-PROC DECL UNKNOWN 0 0 SUBROUTINE) (
UNKNOWN 0 0 0 UNKNOWN ()) 727 0 (728 729) () 0 () () () 0 0)
730 'set_halo_receive' 'halos_base' 'set_halo_receive' 1 ((PROCEDURE
UNKNOWN-INTENT MODULE-PROC DECL UNKNOWN 0 0 SUBROUTINE) (UNKNOWN 0 0 0
UNKNOWN ()) 731 0 (732 733 734 735) () 0 () () () 0 0)
736 'set_halo_receives' 'halos_base' 'set_halo_receives' 1 ((PROCEDURE
UNKNOWN-INTENT MODULE-PROC DECL UNKNOWN 0 0 SUBROUTINE) (UNKNOWN 0 0 0
UNKNOWN ()) 737 0 (738 739 740) () 0 () () () 0 0)
741 'set_halo_send' 'halos_base' 'set_halo_send' 1 ((PROCEDURE
UNKNOWN-INTENT MODULE-PROC DECL UNKNOWN 0 0 SUBROUTINE) (UNKNOWN 0 0 0
UNKNOWN ()) 742 0 (743 744 745 746) () 0 () () () 0 0)
747 'set_halo_sends' 'halos_base' 'set_halo_sends' 1 ((PROCEDURE
UNKNOWN-INTENT MODULE-PROC DECL UNKNOWN 0 0 SUBROUTINE) (UNKNOWN 0 0 0
UNKNOWN ()) 748 0 (749 750 751) () 0 () () () 0 0)
752 'set_halo_universal_number' 'halos_numbering'
'set_halo_universal_number' 1 ((PROCEDURE UNKNOWN-INTENT MODULE-PROC
DECL UNKNOWN 0 0 SUBROUTINE ALWAYS_EXPLICIT) (UNKNOWN 0 0 0 UNKNOWN ())
753 0 (754 755 756 757) () 0 () () () 0 0)
758 'trailing_receives_consistent' 'halos_debug'
'trailing_receives_consistent' 1 ((PROCEDURE UNKNOWN-INTENT MODULE-PROC
DECL UNKNOWN 0 0 FUNCTION) (LOGICAL 4 0 0 LOGICAL ()) 759 0 (760) () 761
() () () 0 0)
762 'universal_numbering_count' 'halos_numbering'
'universal_numbering_count' 1 ((PROCEDURE UNKNOWN-INTENT MODULE-PROC
DECL UNKNOWN 0 0 FUNCTION PURE) (INTEGER 4 0 0 INTEGER ()) 763 0 (764) ()
765 () () () 0 0)
766 'valid_global_to_universal_numbering' 'halos_numbering'
'valid_global_to_universal_numbering' 1 ((PROCEDURE UNKNOWN-INTENT
MODULE-PROC DECL UNKNOWN 0 0 FUNCTION) (LOGICAL 4 0 0 LOGICAL ()) 767 0
(768) () 769 () () () 0 0)
770 'valid_halo_communicator' 'halos_debug' 'valid_halo_communicator' 1
((PROCEDURE UNKNOWN-INTENT MODULE-PROC DECL UNKNOWN 0 0 FUNCTION) (
LOGICAL 4 0 0 LOGICAL ()) 771 0 (772) () 773 () () () 0 0)
774 'valid_halo_node_counts' 'halos_debug' 'valid_halo_node_counts' 1 (
(PROCEDURE UNKNOWN-INTENT MODULE-PROC DECL UNKNOWN 0 0 FUNCTION) (
LOGICAL 4 0 0 LOGICAL ()) 775 0 (776) () 777 () () () 0 0)
778 'valid_serial_halo' 'halos_debug' 'valid_serial_halo' 1 ((PROCEDURE
UNKNOWN-INTENT MODULE-PROC DECL UNKNOWN 0 0 FUNCTION IMPLICIT_PURE) (
LOGICAL 4 0 0 LOGICAL ()) 779 0 (780) () 781 () () () 0 0)
782 'verify_halos' 'halos_registration' 'verify_halos' 1 ((PROCEDURE
UNKNOWN-INTENT MODULE-PROC DECL UNKNOWN 0 0 SUBROUTINE) (UNKNOWN 0 0 0
UNKNOWN ()) 783 0 (784) () 0 () () () 0 0)
785 'write_halos' 'halos_registration' 'write_halos' 1 ((PROCEDURE
UNKNOWN-INTENT MODULE-PROC DECL UNKNOWN 0 0 SUBROUTINE) (UNKNOWN 0 0 0
UNKNOWN ()) 786 0 (787 788) () 0 () () () 0 0)
789 'write_universal_numbering' 'halos_diagnostics'
'write_universal_numbering' 1 ((PROCEDURE UNKNOWN-INTENT MODULE-PROC
DECL UNKNOWN 0 0 SUBROUTINE) (UNKNOWN 0 0 0 UNKNOWN ()) 790 0 (791 792
793 794) () 0 () () () 0 0)
682 'priority' '' 'priority' 680 ((VARIABLE IN UNKNOWN-PROC UNKNOWN
UNKNOWN 0 0 DUMMY) (INTEGER 4 0 0 INTEGER ()) 0 0 () () 0 () () () 0 0)
436 'halo' '' 'halo' 435 ((VARIABLE INOUT UNKNOWN-PROC UNKNOWN UNKNOWN 0
0 DUMMY) (DERIVED 131 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
241 'nsends' '' 'nsends' 239 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN
0 0 DIMENSION DUMMY) (INTEGER 4 0 0 INTEGER ()) 0 0 () (1 0
ASSUMED_SHAPE (CONSTANT (INTEGER 4 0 0 INTEGER ()) 0 '1') ()) 0 () () ()
0 0)
240 'halo' '' 'halo' 239 ((VARIABLE OUT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DUMMY) (DERIVED 131 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
244 'communicator' '' 'communicator' 239 ((VARIABLE IN UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0 OPTIONAL DUMMY) (INTEGER 4 0 0 INTEGER ()) 0 0 () ()
0 () () () 0 0)
243 'name' '' 'name' 239 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
OPTIONAL DUMMY) (CHARACTER 1 0 0 CHARACTER (())) 0 0 () () 0 () () () 0
0)
242 'nreceives' '' 'nreceives' 239 ((VARIABLE IN UNKNOWN-PROC UNKNOWN
UNKNOWN 0 0 DIMENSION DUMMY) (INTEGER 4 0 0 INTEGER ()) 0 0 () (1 0
ASSUMED_SHAPE (CONSTANT (INTEGER 4 0 0 INTEGER ()) 0 '1') ()) 0 () () ()
0 0)
67 'halos' '' 'halos' 66 ((VARIABLE INOUT UNKNOWN-PROC UNKNOWN UNKNOWN 0
0 DIMENSION DUMMY) (DERIVED 131 0 0 DERIVED ()) 0 0 () (1 0
ASSUMED_SHAPE (CONSTANT (INTEGER 4 0 0 INTEGER ()) 0 '1') ()) 0 () () ()
0 0)
248 'ordering_scheme' '' 'ordering_scheme' 239 ((VARIABLE IN
UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 OPTIONAL DUMMY) (INTEGER 4 0 0 INTEGER
()) 0 0 () () 0 () () () 0 0)
88 'object' '' 'object' 87 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0
0 TARGET DUMMY) (DERIVED 131 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
69 'object' '' 'object' 68 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0
0 DUMMY) (DERIVED 131 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
70 'has_references' '' 'has_references' 68 ((VARIABLE UNKNOWN-INTENT
UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 RESULT) (LOGICAL 4 0 0 LOGICAL ()) 0 0
() () 0 () () () 0 0)
76 'halo' '' 'halo' 75 ((VARIABLE INOUT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DUMMY) (DERIVED 131 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
247 'data_type' '' 'data_type' 239 ((VARIABLE IN UNKNOWN-PROC UNKNOWN
UNKNOWN 0 0 OPTIONAL DUMMY) (INTEGER 4 0 0 INTEGER ()) 0 0 () () 0 () ()
() 0 0)
246 'nowned_nodes' '' 'nowned_nodes' 239 ((VARIABLE IN UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0 OPTIONAL DUMMY) (INTEGER 4 0 0 INTEGER ()) 0 0 () ()
0 () () () 0 0)
425 'output_halo' '' 'output_halo' 424 ((VARIABLE OUT UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0 DUMMY) (DERIVED 131 0 0 DERIVED ()) 0 0 () () 0 () ()
() 0 0)
73 'nsends' '' 'nsends' 71 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0
0 DIMENSION OPTIONAL DUMMY) (INTEGER 4 0 0 INTEGER ()) 0 0 () (1 0
EXPLICIT (CONSTANT (INTEGER 4 0 0 INTEGER ()) 0 '1') (FUNCTION (INTEGER
4 0 0 INTEGER ()) 0 384 (('' (VARIABLE (DERIVED 131 0 0 DERIVED ()) 0 72
()))) 'halo_proc_count' 1 384)) 0 () () () 0 0)
72 'halo' '' 'halo' 71 ((VARIABLE INOUT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DUMMY) (DERIVED 131 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
190 'halo' '' 'halo' 189 ((VARIABLE INOUT UNKNOWN-PROC UNKNOWN UNKNOWN 0
0 DUMMY) (DERIVED 131 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
74 'nreceives' '' 'nreceives' 71 ((VARIABLE IN UNKNOWN-PROC UNKNOWN
UNKNOWN 0 0 DIMENSION OPTIONAL DUMMY) (INTEGER 4 0 0 INTEGER ()) 0 0 ()
(1 0 EXPLICIT (CONSTANT (INTEGER 4 0 0 INTEGER ()) 0 '1') (FUNCTION (
INTEGER 4 0 0 INTEGER ()) 0 384 (('' (VARIABLE (DERIVED 131 0 0 DERIVED
()) 0 72 ()))) 'halo_proc_count' 1 384)) 0 () () () 0 0)
386 'halo' '' 'halo' 385 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DUMMY) (DERIVED 131 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
426 'base_halo' '' 'base_halo' 424 ((VARIABLE IN UNKNOWN-PROC UNKNOWN
UNKNOWN 0 0 DUMMY) (DERIVED 131 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
245 'nprocs' '' 'nprocs' 239 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN
0 0 OPTIONAL DUMMY) (INTEGER 4 0 0 INTEGER ()) 0 0 () () 0 () () () 0 0)
439 'halo' '' 'halo' 438 ((VARIABLE INOUT UNKNOWN-PROC UNKNOWN UNKNOWN 0
0 DUMMY) (DERIVED 131 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
772 'halo' '' 'halo' 771 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DUMMY) (DERIVED 131 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
761 'consistent' '' 'consistent' 759 ((VARIABLE UNKNOWN-INTENT
UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 RESULT) (LOGICAL 4 0 0 LOGICAL ()) 0 0
() () 0 () () () 0 0)
637 'valid' '' 'valid' 635 ((VARIABLE UNKNOWN-INTENT UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0 RESULT) (LOGICAL 4 0 0 LOGICAL ()) 0 0 () () 0 () ()
() 0 0)
780 'halo' '' 'halo' 779 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DUMMY) (DERIVED 131 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
94 'pending' '' 'pending' 92 ((VARIABLE UNKNOWN-INTENT UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0 RESULT) (LOGICAL 4 0 0 LOGICAL ()) 0 0 () () 0 () ()
() 0 0)
93 'halo' '' 'halo' 92 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DUMMY) (DERIVED 131 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
91 'pending' '' 'pending' 89 ((VARIABLE UNKNOWN-INTENT UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0 RESULT ALWAYS_EXPLICIT) (LOGICAL 4 0 0 LOGICAL ()) 0
0 () () 0 () () () 0 0)
90 'communicator' '' 'communicator' 89 ((VARIABLE IN UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0 OPTIONAL DUMMY) (INTEGER 4 0 0 INTEGER ()) 0 0 () ()
0 () () () 0 0)
781 'valid' '' 'valid' 779 ((VARIABLE UNKNOWN-INTENT UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0 RESULT) (LOGICAL 4 0 0 LOGICAL ()) 0 0 () () 0 () ()
() 0 0)
773 'valid' '' 'valid' 771 ((VARIABLE UNKNOWN-INTENT UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0 RESULT) (LOGICAL 4 0 0 LOGICAL ()) 0 0 () () 0 () ()
() 0 0)
636 'halo' '' 'halo' 635 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DUMMY) (DERIVED 131 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
681 'halo' '' 'halo' 680 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DUMMY) (DERIVED 131 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
760 'halo' '' 'halo' 759 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DUMMY) (DERIVED 131 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
776 'halo' '' 'halo' 775 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DUMMY) (DERIVED 131 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
777 'valid' '' 'valid' 775 ((VARIABLE UNKNOWN-INTENT UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0 RESULT) (LOGICAL 4 0 0 LOGICAL ()) 0 0 () () 0 () ()
() 0 0)
482 'start_indices' '' 'start_indices' 478 ((VARIABLE OUT UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0 DIMENSION OPTIONAL DUMMY) (INTEGER 4 0 0 INTEGER ())
0 0 () (1 0 EXPLICIT (CONSTANT (INTEGER 4 0 0 INTEGER ()) 0 '1') (
FUNCTION (INTEGER 4 0 0 INTEGER ()) 0 384 (('' (VARIABLE (DERIVED 131 0
0 DERIVED ()) 0 479 ()))) 'halo_proc_count' 1 384)) 0 () () () 0 0)
528 'halo' '' 'halo' 527 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DUMMY) (DERIVED 131 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
531 'halo' '' 'halo' 530 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DUMMY) (DERIVED 131 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
481 'nsends' '' 'nsends' 478 ((VARIABLE OUT UNKNOWN-PROC UNKNOWN UNKNOWN
0 0 DIMENSION OPTIONAL DUMMY) (INTEGER 4 0 0 INTEGER ()) 0 0 () (1 0
EXPLICIT (CONSTANT (INTEGER 4 0 0 INTEGER ()) 0 '1') (FUNCTION (INTEGER
4 0 0 INTEGER ()) 0 384 (('' (VARIABLE (DERIVED 131 0 0 DERIVED ()) 0
479 ()))) 'halo_proc_count' 1 384)) 0 () () () 0 0)
480 'sends' '' 'sends' 478 ((VARIABLE OUT UNKNOWN-PROC UNKNOWN UNKNOWN 0
0 DIMENSION DUMMY) (INTEGER 4 0 0 INTEGER ()) 0 0 () (1 0 EXPLICIT (
CONSTANT (INTEGER 4 0 0 INTEGER ()) 0 '1') (FUNCTION (INTEGER 4 0 0
INTEGER ()) 0 526 (('' (VARIABLE (DERIVED 131 0 0 DERIVED ()) 0 479 ())))
'halo_all_sends_count' 1 526)) 0 () () () 0 0)
537 'halo' '' 'halo' 539 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DUMMY) (DERIVED 131 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
536 'len_trim' '(intrinsic)' 'len_trim' 539 ((PROCEDURE UNKNOWN-INTENT
UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 FUNCTION) (INTEGER 4 0 0 INTEGER ()) 0
0 () () 536 () () () 0 0)
560 'halo' '' 'halo' 559 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DUMMY) (DERIVED 131 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
566 'process' '' 'process' 564 ((VARIABLE IN UNKNOWN-PROC UNKNOWN
UNKNOWN 0 0 DUMMY) (INTEGER 4 0 0 INTEGER ()) 0 0 () () 0 () () () 0 0)
567 'index' '' 'index' 564 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0
0 DUMMY) (INTEGER 4 0 0 INTEGER ()) 0 0 () () 0 () () () 0 0)
571 'process' '' 'process' 569 ((VARIABLE IN UNKNOWN-PROC UNKNOWN
UNKNOWN 0 0 DUMMY) (INTEGER 4 0 0 INTEGER ()) 0 0 () () 0 () () () 0 0)
565 'halo' '' 'halo' 564 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DUMMY) (DERIVED 131 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
575 'nreceives' '' 'nreceives' 573 ((VARIABLE OUT UNKNOWN-PROC UNKNOWN
UNKNOWN 0 0 DIMENSION DUMMY) (INTEGER 4 0 0 INTEGER ()) 0 0 () (1 0
EXPLICIT (CONSTANT (INTEGER 4 0 0 INTEGER ()) 0 '1') (FUNCTION (INTEGER
4 0 0 INTEGER ()) 0 384 (('' (VARIABLE (DERIVED 131 0 0 DERIVED ()) 0
574 ()))) 'halo_proc_count' 1 384)) 0 () () () 0 0)
579 'process' '' 'process' 577 ((VARIABLE IN UNKNOWN-PROC UNKNOWN
UNKNOWN 0 0 DUMMY) (INTEGER 4 0 0 INTEGER ()) 0 0 () () 0 () () () 0 0)
578 'halo' '' 'halo' 577 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DUMMY) (DERIVED 131 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
583 'process' '' 'process' 581 ((VARIABLE IN UNKNOWN-PROC UNKNOWN
UNKNOWN 0 0 DUMMY) (INTEGER 4 0 0 INTEGER ()) 0 0 () () 0 () () () 0 0)
591 'halo' '' 'halo' 590 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DUMMY) (DERIVED 131 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
595 'halo' '' 'halo' 594 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DUMMY) (DERIVED 131 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
592 'nsends' '' 'nsends' 590 ((VARIABLE OUT UNKNOWN-PROC UNKNOWN UNKNOWN
0 0 DIMENSION DUMMY) (INTEGER 4 0 0 INTEGER ()) 0 0 () (1 0 EXPLICIT (
CONSTANT (INTEGER 4 0 0 INTEGER ()) 0 '1') (FUNCTION (INTEGER 4 0 0
INTEGER ()) 0 384 (('' (VARIABLE (DERIVED 131 0 0 DERIVED ()) 0 591 ())))
'halo_proc_count' 1 384)) 0 () () () 0 0)
584 'index' '' 'index' 581 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0
0 DUMMY) (INTEGER 4 0 0 INTEGER ()) 0 0 () () 0 () () () 0 0)
645 'halo' '' 'halo' 644 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DUMMY) (DERIVED 131 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
617 'process' '' 'process' 615 ((VARIABLE IN UNKNOWN-PROC UNKNOWN
UNKNOWN 0 0 DUMMY) (INTEGER 4 0 0 INTEGER ()) 0 0 () () 0 () () () 0 0)
596 'process' '' 'process' 594 ((VARIABLE IN UNKNOWN-PROC UNKNOWN
UNKNOWN 0 0 DUMMY) (INTEGER 4 0 0 INTEGER ()) 0 0 () () 0 () () () 0 0)
658 'max_node' '' 'max_node' 656 ((VARIABLE UNKNOWN-INTENT UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0 RESULT) (INTEGER 4 0 0 INTEGER ()) 0 0 () () 0 () ()
() 0 0)
662 'max_node' '' 'max_node' 660 ((VARIABLE UNKNOWN-INTENT UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0 RESULT) (INTEGER 4 0 0 INTEGER ()) 0 0 () () 0 () ()
() 0 0)
661 'halo' '' 'halo' 660 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DUMMY) (DERIVED 131 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
657 'halo' '' 'halo' 656 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DUMMY) (DERIVED 131 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
582 'halo' '' 'halo' 581 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DUMMY) (DERIVED 131 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
574 'halo' '' 'halo' 573 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DUMMY) (DERIVED 131 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
666 'max_node' '' 'max_node' 664 ((VARIABLE UNKNOWN-INTENT UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0 RESULT) (INTEGER 4 0 0 INTEGER ()) 0 0 () () 0 () ()
() 0 0)
704 'halo' '' 'halo' 703 ((VARIABLE INOUT UNKNOWN-PROC UNKNOWN UNKNOWN 0
0 DUMMY) (DERIVED 131 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
708 'halo' '' 'halo' 707 ((VARIABLE INOUT UNKNOWN-PROC UNKNOWN UNKNOWN 0
0 DUMMY) (DERIVED 131 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
705 'receives' '' 'receives' 703 ((VARIABLE IN UNKNOWN-PROC UNKNOWN
UNKNOWN 0 0 DIMENSION DUMMY) (INTEGER 4 0 0 INTEGER ()) 0 0 () (1 0
ASSUMED_SHAPE (CONSTANT (INTEGER 4 0 0 INTEGER ()) 0 '1') ()) 0 () () ()
0 0)
713 'communicator' '' 'communicator' 711 ((VARIABLE IN UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0 DUMMY) (INTEGER 4 0 0 INTEGER ()) 0 0 () () 0 () ()
() 0 0)
717 'data_type' '' 'data_type' 715 ((VARIABLE IN UNKNOWN-PROC UNKNOWN
UNKNOWN 0 0 DUMMY) (INTEGER 4 0 0 INTEGER ()) 0 0 () () 0 () () () 0 0)
720 'halo' '' 'halo' 719 ((VARIABLE INOUT UNKNOWN-PROC UNKNOWN UNKNOWN 0
0 DUMMY) (DERIVED 131 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
721 'name' '' 'name' 719 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DUMMY) (CHARACTER 1 0 0 CHARACTER (())) 0 0 () () 0 () () () 0 0)
728 'halo' '' 'halo' 727 ((VARIABLE INOUT UNKNOWN-PROC UNKNOWN UNKNOWN 0
0 DUMMY) (DERIVED 131 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
716 'halo' '' 'halo' 715 ((VARIABLE INOUT UNKNOWN-PROC UNKNOWN UNKNOWN 0
0 DUMMY) (DERIVED 131 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
709 'sends' '' 'sends' 707 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0
0 DIMENSION DUMMY) (INTEGER 4 0 0 INTEGER ()) 0 0 () (1 0 ASSUMED_SHAPE
(CONSTANT (INTEGER 4 0 0 INTEGER ()) 0 '1') ()) 0 () () () 0 0)
733 'process' '' 'process' 731 ((VARIABLE IN UNKNOWN-PROC UNKNOWN
UNKNOWN 0 0 DUMMY) (INTEGER 4 0 0 INTEGER ()) 0 0 () () 0 () () () 0 0)
732 'halo' '' 'halo' 731 ((VARIABLE INOUT UNKNOWN-PROC UNKNOWN UNKNOWN 0
0 DUMMY) (DERIVED 131 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
735 'node' '' 'node' 731 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DUMMY) (INTEGER 4 0 0 INTEGER ()) 0 0 () () 0 () () () 0 0)
734 'index' '' 'index' 731 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0
0 DUMMY) (INTEGER 4 0 0 INTEGER ()) 0 0 () () 0 () () () 0 0)
739 'process' '' 'process' 737 ((VARIABLE IN UNKNOWN-PROC UNKNOWN
UNKNOWN 0 0 DUMMY) (INTEGER 4 0 0 INTEGER ()) 0 0 () () 0 () () () 0 0)
738 'halo' '' 'halo' 737 ((VARIABLE INOUT UNKNOWN-PROC UNKNOWN UNKNOWN 0
0 DUMMY) (DERIVED 131 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
740 'receives' '' 'receives' 737 ((VARIABLE IN UNKNOWN-PROC UNKNOWN
UNKNOWN 0 0 DIMENSION DUMMY) (INTEGER 4 0 0 INTEGER ()) 0 0 () (1 0
EXPLICIT (CONSTANT (INTEGER 4 0 0 INTEGER ()) 0 '1') (FUNCTION (INTEGER
4 0 0 INTEGER ()) 0 568 (('' (VARIABLE (DERIVED 131 0 0 DERIVED ()) 0
738 ())) ('' (VARIABLE (INTEGER 4 0 0 INTEGER ()) 0 739 ())))
'halo_receive_count' 1 568)) 0 () () () 0 0)
744 'process' '' 'process' 742 ((VARIABLE IN UNKNOWN-PROC UNKNOWN
UNKNOWN 0 0 DUMMY) (INTEGER 4 0 0 INTEGER ()) 0 0 () () 0 () () () 0 0)
746 'node' '' 'node' 742 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DUMMY) (INTEGER 4 0 0 INTEGER ()) 0 0 () () 0 () () () 0 0)
749 'halo' '' 'halo' 748 ((VARIABLE INOUT UNKNOWN-PROC UNKNOWN UNKNOWN 0
0 DUMMY) (DERIVED 131 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
751 'sends' '' 'sends' 748 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0
0 DIMENSION DUMMY) (INTEGER 4 0 0 INTEGER ()) 0 0 () (1 0 EXPLICIT (
CONSTANT (INTEGER 4 0 0 INTEGER ()) 0 '1') (FUNCTION (INTEGER 4 0 0
INTEGER ()) 0 585 (('' (VARIABLE (DERIVED 131 0 0 DERIVED ()) 0 749 ()))
('' (VARIABLE (INTEGER 4 0 0 INTEGER ()) 0 750 ()))) 'halo_send_count' 1
585)) 0 () () () 0 0)
750 'process' '' 'process' 748 ((VARIABLE IN UNKNOWN-PROC UNKNOWN
UNKNOWN 0 0 DUMMY) (INTEGER 4 0 0 INTEGER ()) 0 0 () () 0 () () () 0 0)
745 'index' '' 'index' 742 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0
0 DUMMY) (INTEGER 4 0 0 INTEGER ()) 0 0 () () 0 () () () 0 0)
743 'halo' '' 'halo' 742 ((VARIABLE INOUT UNKNOWN-PROC UNKNOWN UNKNOWN 0
0 DUMMY) (DERIVED 131 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
411 'halo' '' 'halo' 410 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DUMMY) (DERIVED 131 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
412 'communicator' '' 'communicator' 410 ((VARIABLE UNKNOWN-INTENT
UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 RESULT) (INTEGER 4 0 0 INTEGER ()) 0 0
() () 0 () () () 0 0)
389 'node_count' '' 'node_count' 387 ((VARIABLE UNKNOWN-INTENT
UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 RESULT) (INTEGER 4 0 0 INTEGER ()) 0 0
() () 0 () () () 0 0)
388 'halo' '' 'halo' 387 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DUMMY) (DERIVED 131 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
414 'halo' '' 'halo' 413 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DUMMY) (DERIVED 131 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
404 'size' '(intrinsic)' 'size' 402 ((PROCEDURE UNKNOWN-INTENT
UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 FUNCTION) (REAL 8 0 0 REAL ()) 0 0 () ()
404 () () () 0 0)
405 'serial' '' 'serial' 402 ((VARIABLE UNKNOWN-INTENT UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0 DIMENSION RESULT ALWAYS_EXPLICIT) (LOGICAL 4 0 0
LOGICAL ()) 0 0 () (1 0 EXPLICIT (CONSTANT (INTEGER 4 0 0 INTEGER ()) 0
'1') (FUNCTION (INTEGER 4 0 0 INTEGER ()) 0 404 (('' (VARIABLE (DERIVED
131 0 0 DERIVED ()) 1 403 ((ARRAY (FULL 0))))) ('' ()) ('' ())) '' 0
'size')) 0 () () () 0 0)
403 'halos' '' 'halos' 402 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0
0 DIMENSION DUMMY) (DERIVED 131 0 0 DERIVED ()) 0 0 () (1 0
ASSUMED_SHAPE (CONSTANT (INTEGER 4 0 0 INTEGER ()) 0 '1') ()) 0 () () ()
0 0)
415 'serial' '' 'serial' 413 ((VARIABLE UNKNOWN-INTENT UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0 RESULT) (LOGICAL 4 0 0 LOGICAL ()) 0 0 () () 0 () ()
() 0 0)
423 'halo' '' 'halo' 422 ((VARIABLE INOUT UNKNOWN-PROC UNKNOWN UNKNOWN 0
0 DUMMY) (DERIVED 131 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
729 'ordering_scheme' '' 'ordering_scheme' 727 ((VARIABLE IN
UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 DUMMY) (INTEGER 4 0 0 INTEGER ()) 0 0 ()
() 0 () () () 0 0)
665 'halo' '' 'halo' 664 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DUMMY) (DERIVED 131 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
474 'receives' '' 'receives' 472 ((VARIABLE OUT UNKNOWN-PROC UNKNOWN
UNKNOWN 0 0 DIMENSION DUMMY) (INTEGER 4 0 0 INTEGER ()) 0 0 () (1 0
EXPLICIT (CONSTANT (INTEGER 4 0 0 INTEGER ()) 0 '1') (FUNCTION (INTEGER
4 0 0 INTEGER ()) 0 523 (('' (VARIABLE (DERIVED 131 0 0 DERIVED ()) 0
473 ()))) 'halo_all_receives_count' 1 523)) 0 () () () 0 0)
475 'nreceives' '' 'nreceives' 472 ((VARIABLE OUT UNKNOWN-PROC UNKNOWN
UNKNOWN 0 0 DIMENSION OPTIONAL DUMMY) (INTEGER 4 0 0 INTEGER ()) 0 0 ()
(1 0 EXPLICIT (CONSTANT (INTEGER 4 0 0 INTEGER ()) 0 '1') (FUNCTION (
INTEGER 4 0 0 INTEGER ()) 0 384 (('' (VARIABLE (DERIVED 131 0 0 DERIVED
()) 0 473 ()))) 'halo_proc_count' 1 384)) 0 () () () 0 0)
476 'start_indices' '' 'start_indices' 472 ((VARIABLE OUT UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0 DIMENSION OPTIONAL DUMMY) (INTEGER 4 0 0 INTEGER ())
0 0 () (1 0 EXPLICIT (CONSTANT (INTEGER 4 0 0 INTEGER ()) 0 '1') (
FUNCTION (INTEGER 4 0 0 INTEGER ()) 0 384 (('' (VARIABLE (DERIVED 131 0
0 DERIVED ()) 0 473 ()))) 'halo_proc_count' 1 384)) 0 () () () 0 0)
473 'halo' '' 'halo' 472 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DUMMY) (DERIVED 131 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
479 'halo' '' 'halo' 478 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DUMMY) (DERIVED 131 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
525 'halo' '' 'halo' 524 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DUMMY) (DERIVED 131 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
534 'halo' '' 'halo' 533 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DUMMY) (DERIVED 131 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
555 'halo' '' 'halo' 554 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DUMMY) (DERIVED 131 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
570 'halo' '' 'halo' 569 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DUMMY) (DERIVED 131 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
587 'halo' '' 'halo' 586 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DUMMY) (DERIVED 131 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
588 'process' '' 'process' 586 ((VARIABLE IN UNKNOWN-PROC UNKNOWN
UNKNOWN 0 0 DUMMY) (INTEGER 4 0 0 INTEGER ()) 0 0 () () 0 () () () 0 0)
616 'halo' '' 'halo' 615 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DUMMY) (DERIVED 131 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
670 'min_node' '' 'min_node' 668 ((VARIABLE UNKNOWN-INTENT UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0 RESULT) (INTEGER 4 0 0 INTEGER ()) 0 0 () () 0 () ()
() 0 0)
669 'halo' '' 'halo' 668 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DUMMY) (DERIVED 131 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
674 'min_node' '' 'min_node' 672 ((VARIABLE UNKNOWN-INTENT UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0 RESULT) (INTEGER 4 0 0 INTEGER ()) 0 0 () () 0 () ()
() 0 0)
673 'halo' '' 'halo' 672 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DUMMY) (DERIVED 131 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
678 'min_node' '' 'min_node' 676 ((VARIABLE UNKNOWN-INTENT UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0 RESULT) (INTEGER 4 0 0 INTEGER ()) 0 0 () () 0 () ()
() 0 0)
677 'halo' '' 'halo' 676 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DUMMY) (DERIVED 131 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
712 'halo' '' 'halo' 711 ((VARIABLE INOUT UNKNOWN-PROC UNKNOWN UNKNOWN 0
0 DUMMY) (DERIVED 131 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
725 'nowned_nodes' '' 'nowned_nodes' 723 ((VARIABLE IN UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0 DUMMY) (INTEGER 4 0 0 INTEGER ()) 0 0 () () 0 () ()
() 0 0)
724 'halo' '' 'halo' 723 ((VARIABLE INOUT UNKNOWN-PROC UNKNOWN UNKNOWN 0
0 DUMMY) (DERIVED 131 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
794 'name' '' 'name' 790 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DUMMY) (CHARACTER 1 0 0 CHARACTER (())) 0 0 () () 0 () () () 0 0)
793 'position' '' 'position' 790 ((VARIABLE IN UNKNOWN-PROC UNKNOWN
UNKNOWN 0 0 DUMMY) (DERIVED 95 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
792 'mesh' '' 'mesh' 790 ((VARIABLE INOUT UNKNOWN-PROC UNKNOWN UNKNOWN 0
0 DUMMY) (DERIVED 105 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
791 'halo' '' 'halo' 790 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DUMMY) (DERIVED 131 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
443 'ordering_scheme' '' 'ordering_scheme' 441 ((VARIABLE IN
UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 OPTIONAL DUMMY) (INTEGER 4 0 0 INTEGER
()) 0 0 () () 0 () () () 0 0)
448 'element_halo' '' 'element_halo' 446 ((VARIABLE IN UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0 DUMMY) (DERIVED 131 0 0 DERIVED ()) 0 0 () () 0 () ()
() 0 0)
447 'mesh' '' 'mesh' 446 ((VARIABLE INOUT UNKNOWN-PROC UNKNOWN UNKNOWN 0
0 DUMMY) (DERIVED 105 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
444 'create_caches' '' 'create_caches' 441 ((VARIABLE IN UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0 OPTIONAL DUMMY) (LOGICAL 4 0 0 LOGICAL ()) 0 0 () ()
0 () () () 0 0)
450 'create_caches' '' 'create_caches' 446 ((VARIABLE IN UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0 OPTIONAL DUMMY) (LOGICAL 4 0 0 LOGICAL ()) 0 0 () ()
0 () () () 0 0)
449 'ordering_scheme' '' 'ordering_scheme' 446 ((VARIABLE IN
UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 OPTIONAL DUMMY) (INTEGER 4 0 0 INTEGER
()) 0 0 () () 0 () () () 0 0)
442 'mesh' '' 'mesh' 441 ((VARIABLE INOUT UNKNOWN-PROC UNKNOWN UNKNOWN 0
0 DUMMY) (DERIVED 105 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
451 'selement_halo' '' 'selement_halo' 446 ((VARIABLE UNKNOWN-INTENT
UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 RESULT ALWAYS_EXPLICIT) (DERIVED 131 0
0 DERIVED ()) 0 0 () () 0 () () () 0 0)
456 'aliased_to_new_node_number' '' 'aliased_to_new_node_number' 453 ((
VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 DUMMY) (DERIVED 294 0 0
DERIVED ()) 0 0 () () 0 () () () 0 0)
459 'halo' '' 'halo' 458 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DUMMY) (DERIVED 131 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
455 'model' '' 'model' 453 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0
0 DUMMY) (DERIVED 95 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
461 'sub_halo' '' 'sub_halo' 458 ((VARIABLE UNKNOWN-INTENT UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0 RESULT ALWAYS_EXPLICIT) (DERIVED 131 0 0 DERIVED ())
0 0 () () 0 () () () 0 0)
460 'sub_nodes' '' 'sub_nodes' 458 ((VARIABLE IN UNKNOWN-PROC UNKNOWN
UNKNOWN 0 0 DIMENSION DUMMY) (INTEGER 4 0 0 INTEGER ()) 0 0 () (1 0
ASSUMED_SHAPE (CONSTANT (INTEGER 4 0 0 INTEGER ()) 0 '1') ()) 0 () () ()
0 0)
291 'mesh' '' 'mesh' 290 ((VARIABLE INOUT UNKNOWN-PROC UNKNOWN UNKNOWN 0
0 DUMMY) (DERIVED 105 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
292 'ordering_scheme' '' 'ordering_scheme' 290 ((VARIABLE IN
UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 OPTIONAL DUMMY) (INTEGER 4 0 0 INTEGER
()) 0 0 () () 0 () () () 0 0)
293 'create_caches' '' 'create_caches' 290 ((VARIABLE IN UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0 OPTIONAL DUMMY) (LOGICAL 4 0 0 LOGICAL ()) 0 0 () ()
0 () () () 0 0)
654 'unknowns_sizes' '' 'unknowns_sizes' 650 ((VARIABLE UNKNOWN-INTENT
UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 DIMENSION RESULT ALWAYS_EXPLICIT) (
INTEGER 4 0 0 INTEGER ()) 0 0 () (1 0 EXPLICIT (CONSTANT (INTEGER 4 0 0
INTEGER ()) 0 '1') (FUNCTION (INTEGER 4 0 0 INTEGER ()) 0 653 (('' (
VARIABLE (INTEGER 4 0 0 INTEGER ()) 1 651 ((ARRAY (FULL 0))))) ('' ()) (
'' ())) '' 0 'size')) 0 () () () 0 0)
454 'new_positions' '' 'new_positions' 453 ((VARIABLE INOUT UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0 DUMMY) (DERIVED 95 0 0 DERIVED ()) 0 0 () () 0 () ()
() 0 0)
464 'ele' '' 'ele' 463 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DUMMY) (INTEGER 4 0 0 INTEGER ()) 0 0 () () 0 () () () 0 0)
466 'node_halo' '' 'node_halo' 463 ((VARIABLE IN UNKNOWN-PROC UNKNOWN
UNKNOWN 0 0 DUMMY) (DERIVED 131 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
651 'knowns_sizes' '' 'knowns_sizes' 650 ((VARIABLE IN UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0 DIMENSION DUMMY) (INTEGER 4 0 0 INTEGER ()) 0 0 () (
1 0 ASSUMED_SHAPE (CONSTANT (INTEGER 4 0 0 INTEGER ()) 0 '1') ()) 0 () ()
() 0 0)
653 'size' '(intrinsic)' 'size' 650 ((PROCEDURE UNKNOWN-INTENT
UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 FUNCTION) (REAL 8 0 0 REAL ()) 0 0 () ()
653 () () () 0 0)
652 'communicator' '' 'communicator' 650 ((VARIABLE IN UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0 OPTIONAL DUMMY) (INTEGER 4 0 0 INTEGER ()) 0 0 () ()
0 () () () 0 0)
465 'mesh' '' 'mesh' 463 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DUMMY) (DERIVED 105 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
328 'halo' '' 'halo' 327 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DUMMY) (DERIVED 131 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
329 'integer_data' '' 'integer_data' 327 ((VARIABLE INOUT UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0 DIMENSION DUMMY) (INTEGER 4 0 0 INTEGER ()) 0 0 () (
1 0 ASSUMED_SHAPE (CONSTANT (INTEGER 4 0 0 INTEGER ()) 0 '1') ()) 0 () ()
() 0 0)
332 'integer_data' '' 'integer_data' 330 ((VARIABLE INOUT UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0 DIMENSION DUMMY) (INTEGER 4 0 0 INTEGER ()) 0 0 () (
3 0 ASSUMED_SHAPE (CONSTANT (INTEGER 4 0 0 INTEGER ()) 0 '1') () (
CONSTANT (INTEGER 4 0 0 INTEGER ()) 0 '1') () (CONSTANT (INTEGER 4 0 0
INTEGER ()) 0 '1') ()) 0 () () () 0 0)
331 'halo' '' 'halo' 330 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DUMMY) (DERIVED 131 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
321 'halo' '' 'halo' 320 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DUMMY) (DERIVED 131 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
326 'integer_data' '' 'integer_data' 324 ((VARIABLE INOUT UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0 DIMENSION DUMMY) (INTEGER 4 0 0 INTEGER ()) 0 0 () (
2 0 ASSUMED_SHAPE (CONSTANT (INTEGER 4 0 0 INTEGER ()) 0 '1') () (
CONSTANT (INTEGER 4 0 0 INTEGER ()) 0 '1') ()) 0 () () () 0 0)
334 'halo' '' 'halo' 333 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DUMMY) (DERIVED 131 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
335 'real_data' '' 'real_data' 333 ((VARIABLE INOUT UNKNOWN-PROC UNKNOWN
UNKNOWN 0 0 DIMENSION DUMMY) (REAL 8 0 0 REAL ()) 0 0 () (1 0
ASSUMED_SHAPE (CONSTANT (INTEGER 4 0 0 INTEGER ()) 0 '1') ()) 0 () () ()
0 0)
323 'block_size' '' 'block_size' 320 ((VARIABLE IN UNKNOWN-PROC UNKNOWN
UNKNOWN 0 0 DUMMY) (INTEGER 4 0 0 INTEGER ()) 0 0 () () 0 () () () 0 0)
312 'halo' '' 'halo' 311 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DUMMY) (DERIVED 131 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
322 'integer_data' '' 'integer_data' 320 ((VARIABLE INOUT UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0 DIMENSION DUMMY) (INTEGER 4 0 0 INTEGER ()) 0 0 () (
1 0 ASSUMED_SIZE (CONSTANT (INTEGER 4 0 0 INTEGER ()) 0 '1') ()) 0 () ()
() 0 0)
325 'halo' '' 'halo' 324 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DUMMY) (DERIVED 131 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
316 'real_data' '' 'real_data' 314 ((VARIABLE INOUT UNKNOWN-PROC UNKNOWN
UNKNOWN 0 0 DIMENSION DUMMY) (REAL 8 0 0 REAL ()) 0 0 () (3 0
ASSUMED_SHAPE (CONSTANT (INTEGER 4 0 0 INTEGER ()) 0 '1') () (CONSTANT (
INTEGER 4 0 0 INTEGER ()) 0 '1') () (CONSTANT (INTEGER 4 0 0 INTEGER ())
0 '1') ()) 0 () () () 0 0)
310 's_field' '' 's_field' 308 ((VARIABLE INOUT UNKNOWN-PROC UNKNOWN
UNKNOWN 0 0 DUMMY) (DERIVED 227 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
309 'halo' '' 'halo' 308 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DUMMY) (DERIVED 131 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
315 'halo' '' 'halo' 314 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DUMMY) (DERIVED 131 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
319 'v_field' '' 'v_field' 317 ((VARIABLE INOUT UNKNOWN-PROC UNKNOWN
UNKNOWN 0 0 DUMMY) (DERIVED 95 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
318 'halo' '' 'halo' 317 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DUMMY) (DERIVED 131 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
304 't_field' '' 't_field' 302 ((VARIABLE INOUT UNKNOWN-PROC UNKNOWN
UNKNOWN 0 0 DUMMY) (DERIVED 366 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
297 's_field' '' 's_field' 296 ((VARIABLE INOUT UNKNOWN-PROC UNKNOWN
UNKNOWN 0 0 DUMMY) (DERIVED 227 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
303 'halo' '' 'halo' 302 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DUMMY) (DERIVED 131 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
313 'real_data' '' 'real_data' 311 ((VARIABLE INOUT UNKNOWN-PROC UNKNOWN
UNKNOWN 0 0 DIMENSION DUMMY) (REAL 8 0 0 REAL ()) 0 0 () (2 0
ASSUMED_SHAPE (CONSTANT (INTEGER 4 0 0 INTEGER ()) 0 '1') () (CONSTANT (
INTEGER 4 0 0 INTEGER ()) 0 '1') ()) 0 () () () 0 0)
301 'level' '' 'level' 299 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0
0 OPTIONAL DUMMY) (INTEGER 4 0 0 INTEGER ()) 0 0 () () 0 () () () 0 0)
300 'v_field' '' 'v_field' 299 ((VARIABLE INOUT UNKNOWN-PROC UNKNOWN
UNKNOWN 0 0 DUMMY) (DERIVED 95 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
307 'level' '' 'level' 305 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0
0 OPTIONAL DUMMY) (INTEGER 4 0 0 INTEGER ()) 0 0 () () 0 () () () 0 0)
355 'halo' '' 'halo' 354 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DUMMY) (DERIVED 131 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
356 's_field' '' 's_field' 354 ((VARIABLE INOUT UNKNOWN-PROC UNKNOWN
UNKNOWN 0 0 DUMMY) (DERIVED 227 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
345 's_field' '' 's_field' 344 ((VARIABLE INOUT UNKNOWN-PROC UNKNOWN
UNKNOWN 0 0 DUMMY) (DERIVED 227 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
346 'level' '' 'level' 344 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0
0 OPTIONAL DUMMY) (INTEGER 4 0 0 INTEGER ()) 0 0 () () 0 () () () 0 0)
353 'real_data' '' 'real_data' 351 ((VARIABLE INOUT UNKNOWN-PROC UNKNOWN
UNKNOWN 0 0 DIMENSION DUMMY) (REAL 8 0 0 REAL ()) 0 0 () (1 0
ASSUMED_SHAPE (CONSTANT (INTEGER 4 0 0 INTEGER ()) 0 '1') ()) 0 () () ()
0 0)
343 'verifies' '' 'verifies' 340 ((VARIABLE UNKNOWN-INTENT UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0 RESULT ALWAYS_EXPLICIT) (LOGICAL 4 0 0 LOGICAL ()) 0
0 () () 0 () () () 0 0)
342 'integer_array' '' 'integer_array' 340 ((VARIABLE IN UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0 DIMENSION DUMMY) (INTEGER 4 0 0 INTEGER ()) 0 0 () (
1 0 ASSUMED_SHAPE (CONSTANT (INTEGER 4 0 0 INTEGER ()) 0 '1') ()) 0 () ()
() 0 0)
350 'verifies' '' 'verifies' 347 ((VARIABLE UNKNOWN-INTENT UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0 RESULT ALWAYS_EXPLICIT) (LOGICAL 4 0 0 LOGICAL ()) 0
0 () () 0 () () () 0 0)
349 'real_array' '' 'real_array' 347 ((VARIABLE IN UNKNOWN-PROC UNKNOWN
UNKNOWN 0 0 DIMENSION DUMMY) (REAL 8 0 0 REAL ()) 0 0 () (1 0
ASSUMED_SHAPE (CONSTANT (INTEGER 4 0 0 INTEGER ()) 0 '1') ()) 0 () () ()
0 0)
348 'halo' '' 'halo' 347 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DUMMY) (DERIVED 131 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
341 'halo' '' 'halo' 340 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DUMMY) (DERIVED 131 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
352 'halo' '' 'halo' 351 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DUMMY) (DERIVED 131 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
339 'verifies' '' 'verifies' 336 ((VARIABLE UNKNOWN-INTENT UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0 RESULT) (LOGICAL 4 0 0 LOGICAL ()) 0 0 () () 0 () ()
() 0 0)
338 'sfield' '' 'sfield' 336 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN
0 0 DUMMY) (DERIVED 227 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
359 'vfield' '' 'vfield' 357 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN
0 0 DUMMY) (DERIVED 95 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
358 'halo' '' 'halo' 357 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DUMMY) (DERIVED 131 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
360 'dim' '' 'dim' 357 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DUMMY) (INTEGER 4 0 0 INTEGER ()) 0 0 () () 0 () () () 0 0)
337 'halo' '' 'halo' 336 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DUMMY) (DERIVED 131 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
306 't_field' '' 't_field' 305 ((VARIABLE INOUT UNKNOWN-PROC UNKNOWN
UNKNOWN 0 0 DUMMY) (DERIVED 366 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
363 'halo' '' 'halo' 362 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DUMMY) (DERIVED 131 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
361 'verifies' '' 'verifies' 357 ((VARIABLE UNKNOWN-INTENT UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0 RESULT) (LOGICAL 4 0 0 LOGICAL ()) 0 0 () () 0 () ()
() 0 0)
298 'level' '' 'level' 296 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0
0 OPTIONAL DUMMY) (INTEGER 4 0 0 INTEGER ()) 0 0 () () 0 () () () 0 0)
365 'verifies' '' 'verifies' 362 ((VARIABLE UNKNOWN-INTENT UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0 RESULT) (LOGICAL 4 0 0 LOGICAL ()) 0 0 () () 0 () ()
() 0 0)
364 'vfield' '' 'vfield' 362 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN
0 0 DUMMY) (DERIVED 95 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
469 'halo' '' 'halo' 468 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DUMMY) (DERIVED 131 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
518 'unns' '' 'unns' 516 ((VARIABLE OUT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DIMENSION DUMMY) (INTEGER 4 0 0 INTEGER ()) 0 0 () (1 0 EXPLICIT (
CONSTANT (INTEGER 4 0 0 INTEGER ()) 0 '1') (FUNCTION (INTEGER 4 0 0
INTEGER ()) 0 795 (('' (VARIABLE (DERIVED 131 0 0 DERIVED ()) 0 517 ())))
'node_count_halo' 1 32)) 0 () () () 0 0)
517 'halo' '' 'halo' 516 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DUMMY) (DERIVED 131 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
522 'gnns' '' 'gnns' 520 ((VARIABLE OUT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DUMMY) (DERIVED 294 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
626 'global_number' '' 'global_number' 624 ((VARIABLE IN UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0 DUMMY) (INTEGER 4 0 0 INTEGER ()) 0 0 () () 0 () ()
() 0 0)
625 'halo' '' 'halo' 624 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DUMMY) (DERIVED 131 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
633 'unns' '' 'unns' 629 ((VARIABLE UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN
UNKNOWN 0 0 DIMENSION RESULT ALWAYS_EXPLICIT) (INTEGER 4 0 0 INTEGER ())
0 0 () (1 0 EXPLICIT (CONSTANT (INTEGER 4 0 0 INTEGER ()) 0 '1') (
FUNCTION (INTEGER 4 0 0 INTEGER ()) 0 632 (('' (VARIABLE (INTEGER 4 0 0
INTEGER ()) 1 631 ((ARRAY (FULL 0))))) ('' ()) ('' ())) '' 0 'size')) 0
() () () 0 0)
627 'unn' '' 'unn' 624 ((VARIABLE UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN
UNKNOWN 0 0 RESULT) (INTEGER 4 0 0 INTEGER ()) 0 0 () () 0 () () () 0 0)
521 'halo' '' 'halo' 520 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DUMMY) (DERIVED 131 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
754 'halo' '' 'halo' 753 ((VARIABLE INOUT UNKNOWN-PROC UNKNOWN UNKNOWN 0
0 DUMMY) (DERIVED 131 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
755 'node' '' 'node' 753 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DUMMY) (INTEGER 4 0 0 INTEGER ()) 0 0 () () 0 () () () 0 0)
756 'universal_number' '' 'universal_number' 753 ((VARIABLE IN
UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 DUMMY) (INTEGER 4 0 0 INTEGER ()) 0 0 ()
() 0 () () () 0 0)
769 'valid' '' 'valid' 767 ((VARIABLE UNKNOWN-INTENT UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0 RESULT) (LOGICAL 4 0 0 LOGICAL ()) 0 0 () () 0 () ()
() 0 0)
757 'stat' '' 'stat' 753 ((VARIABLE OUT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
OPTIONAL DUMMY) (INTEGER 4 0 0 INTEGER ()) 0 0 () () 0 () () () 0 0)
382 'size' '(intrinsic)' 'size' 379 ((PROCEDURE UNKNOWN-INTENT
UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 FUNCTION) (REAL 8 0 0 REAL ()) 0 0 () ()
382 () () () 0 0)
381 'global_number' '' 'global_number' 379 ((VARIABLE IN UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0 DIMENSION DUMMY) (INTEGER 4 0 0 INTEGER ()) 0 0 () (
1 0 ASSUMED_SHAPE (CONSTANT (INTEGER 4 0 0 INTEGER ()) 0 '1') ()) 0 () ()
() 0 0)
380 'halo' '' 'halo' 379 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DUMMY) (DERIVED 131 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
377 'halo' '' 'halo' 376 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DUMMY) (DERIVED 131 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
378 'unns' '' 'unns' 376 ((VARIABLE OUT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DIMENSION DUMMY) (INTEGER 4 0 0 INTEGER ()) 0 0 () (2 0 ASSUMED_SHAPE (
CONSTANT (INTEGER 4 0 0 INTEGER ()) 0 '1') () (CONSTANT (INTEGER 4 0 0
INTEGER ()) 0 '1') ()) 0 () () () 0 0)
383 'unn' '' 'unn' 379 ((VARIABLE UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN
UNKNOWN 0 0 DIMENSION RESULT ALWAYS_EXPLICIT) (INTEGER 4 0 0 INTEGER ())
0 0 () (1 0 EXPLICIT (CONSTANT (INTEGER 4 0 0 INTEGER ()) 0 '1') (
FUNCTION (INTEGER 4 0 0 INTEGER ()) 0 382 (('' (VARIABLE (INTEGER 4 0 0
INTEGER ()) 1 381 ((ARRAY (FULL 0))))) ('' ()) ('' ())) '' 0 'size')) 0
() () () 0 0)
470 'debug_level' '' 'debug_level' 468 ((VARIABLE IN UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0 DUMMY) (INTEGER 4 0 0 INTEGER ()) 0 0 () () 0 () ()
() 0 0)
430 'local_only' '' 'local_only' 428 ((VARIABLE IN UNKNOWN-PROC UNKNOWN
UNKNOWN 0 0 OPTIONAL DUMMY) (LOGICAL 4 0 0 LOGICAL ()) 0 0 () () 0 () ()
() 0 0)
429 'halo' '' 'halo' 428 ((VARIABLE INOUT UNKNOWN-PROC UNKNOWN UNKNOWN 0
0 DUMMY) (DERIVED 131 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
631 'global_numbers' '' 'global_numbers' 629 ((VARIABLE IN UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0 DIMENSION DUMMY) (INTEGER 4 0 0 INTEGER ()) 0 0 () (
1 0 ASSUMED_SHAPE (CONSTANT (INTEGER 4 0 0 INTEGER ()) 0 '1') ()) 0 () ()
() 0 0)
630 'halo' '' 'halo' 629 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DUMMY) (DERIVED 131 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
632 'size' '(intrinsic)' 'size' 629 ((PROCEDURE UNKNOWN-INTENT
UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 FUNCTION) (REAL 8 0 0 REAL ()) 0 0 () ()
632 () () () 0 0)
642 'has_gnn_to_unn' '' 'has_gnn_to_unn' 640 ((VARIABLE UNKNOWN-INTENT
UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 RESULT) (LOGICAL 4 0 0 LOGICAL ()) 0 0
() () 0 () () () 0 0)
641 'halo' '' 'halo' 640 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DUMMY) (DERIVED 131 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
764 'halo' '' 'halo' 763 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DUMMY) (DERIVED 131 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
765 'unn_count' '' 'unn_count' 763 ((VARIABLE UNKNOWN-INTENT
UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 RESULT) (INTEGER 4 0 0 INTEGER ()) 0 0
() () 0 () () () 0 0)
768 'halo' '' 'halo' 767 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DUMMY) (DERIVED 131 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
542 'halo' '' 'halo' 541 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DUMMY) (DERIVED 131 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
543 'node' '' 'node' 541 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DUMMY) (INTEGER 4 0 0 INTEGER ()) 0 0 () () 0 () () () 0 0)
545 'node_owner' '' 'node_owner' 541 ((VARIABLE UNKNOWN-INTENT
UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 RESULT ALWAYS_EXPLICIT) (INTEGER 4 0 0
INTEGER ()) 0 0 () () 0 () () () 0 0)
549 'nodes' '' 'nodes' 547 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0
0 DIMENSION DUMMY) (INTEGER 4 0 0 INTEGER ()) 0 0 () (1 0 ASSUMED_SHAPE
(CONSTANT (INTEGER 4 0 0 INTEGER ()) 0 '1') ()) 0 () () () 0 0)
544 'permit_extended' '' 'permit_extended' 541 ((VARIABLE IN
UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 OPTIONAL DUMMY) (LOGICAL 4 0 0 LOGICAL
()) 0 0 () () 0 () () () 0 0)
433 'halo' '' 'halo' 432 ((VARIABLE INOUT UNKNOWN-PROC UNKNOWN UNKNOWN 0
0 DUMMY) (DERIVED 131 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
551 'size' '(intrinsic)' 'size' 547 ((PROCEDURE UNKNOWN-INTENT
UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 FUNCTION) (REAL 8 0 0 REAL ()) 0 0 () ()
551 () () () 0 0)
552 'node_owners' '' 'node_owners' 547 ((VARIABLE UNKNOWN-INTENT
UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 DIMENSION RESULT ALWAYS_EXPLICIT) (
INTEGER 4 0 0 INTEGER ()) 0 0 () (1 0 EXPLICIT (CONSTANT (INTEGER 4 0 0
INTEGER ()) 0 '1') (FUNCTION (INTEGER 4 0 0 INTEGER ()) 0 551 (('' (
VARIABLE (INTEGER 4 0 0 INTEGER ()) 1 549 ((ARRAY (FULL 0))))) ('' ()) (
'' ())) '' 0 'size')) 0 () () () 0 0)
622 'size' '(intrinsic)' 'size' 619 ((PROCEDURE UNKNOWN-INTENT
UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 FUNCTION) (REAL 8 0 0 REAL ()) 0 0 () ()
622 () () () 0 0)
623 'node_owners' '' 'node_owners' 619 ((VARIABLE UNKNOWN-INTENT
UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 DIMENSION RESULT ALWAYS_EXPLICIT) (
INTEGER 4 0 0 INTEGER ()) 0 0 () (1 0 EXPLICIT (CONSTANT (INTEGER 4 0 0
INTEGER ()) 0 '1') (FUNCTION (INTEGER 4 0 0 INTEGER ()) 0 622 (('' (
VARIABLE (INTEGER 4 0 0 INTEGER ()) 1 621 ((ARRAY (FULL 0))))) ('' ()) (
'' ())) '' 0 'size')) 0 () () () 0 0)
648 'halo' '' 'halo' 647 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DUMMY) (DERIVED 131 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
396 'halo' '' 'halo' 395 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DUMMY) (DERIVED 131 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
621 'unns' '' 'unns' 619 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DIMENSION DUMMY) (INTEGER 4 0 0 INTEGER ()) 0 0 () (1 0 ASSUMED_SHAPE (
CONSTANT (INTEGER 4 0 0 INTEGER ()) 0 '1') ()) 0 () () () 0 0)
394 'owned' '' 'owned' 390 ((VARIABLE UNKNOWN-INTENT UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0 DIMENSION RESULT ALWAYS_EXPLICIT) (LOGICAL 4 0 0
LOGICAL ()) 0 0 () (1 0 EXPLICIT (CONSTANT (INTEGER 4 0 0 INTEGER ()) 0
'1') (FUNCTION (INTEGER 4 0 0 INTEGER ()) 0 393 (('' (VARIABLE (INTEGER
4 0 0 INTEGER ()) 1 392 ((ARRAY (FULL 0))))) ('' ()) ('' ())) '' 0 'size'))
0 () () () 0 0)
393 'size' '(intrinsic)' 'size' 390 ((PROCEDURE UNKNOWN-INTENT
UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 FUNCTION) (REAL 8 0 0 REAL ()) 0 0 () ()
393 () () () 0 0)
392 'nodes' '' 'nodes' 390 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0
0 DIMENSION DUMMY) (INTEGER 4 0 0 INTEGER ()) 0 0 () (1 0 ASSUMED_SHAPE
(CONSTANT (INTEGER 4 0 0 INTEGER ()) 0 '1') ()) 0 () () () 0 0)
391 'halo' '' 'halo' 390 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DUMMY) (DERIVED 131 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
397 'node' '' 'node' 395 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DUMMY) (INTEGER 4 0 0 INTEGER ()) 0 0 () () 0 () () () 0 0)
550 'permit_extended' '' 'permit_extended' 547 ((VARIABLE IN
UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 OPTIONAL DUMMY) (LOGICAL 4 0 0 LOGICAL
()) 0 0 () () 0 () () () 0 0)
510 'halo' '' 'halo' 509 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DUMMY) (DERIVED 131 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
511 'owners' '' 'owners' 509 ((VARIABLE OUT UNKNOWN-PROC UNKNOWN UNKNOWN
0 0 DIMENSION DUMMY) (INTEGER 4 0 0 INTEGER ()) 0 0 () (1 0
ASSUMED_SHAPE (CONSTANT (INTEGER 4 0 0 INTEGER ()) 0 '1') ()) 0 () () ()
0 0)
515 'owned_nodes' '' 'owned_nodes' 513 ((VARIABLE OUT UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0 DIMENSION DUMMY) (INTEGER 4 0 0 INTEGER ()) 0 0 () (
1 0 ASSUMED_SHAPE (CONSTANT (INTEGER 4 0 0 INTEGER ()) 0 '1') ()) 0 () ()
() 0 0)
514 'halo' '' 'halo' 513 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DUMMY) (DERIVED 131 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
548 'halo' '' 'halo' 547 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DUMMY) (DERIVED 131 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
620 'halo' '' 'halo' 619 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DUMMY) (DERIVED 131 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
489 'receive_starts' '' 'receive_starts' 484 ((VARIABLE OUT UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0 DIMENSION DUMMY) (INTEGER 4 0 0 INTEGER ()) 0 0 () (
1 0 ASSUMED_SHAPE (CONSTANT (INTEGER 4 0 0 INTEGER ()) 0 '1') ()) 0 () ()
() 0 0)
490 'nowned_nodes' '' 'nowned_nodes' 484 ((VARIABLE OUT UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0 OPTIONAL DUMMY) (INTEGER 4 0 0 INTEGER ()) 0 0 () ()
0 () () () 0 0)
495 'sends' '' 'sends' 492 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0
0 DIMENSION DUMMY) (INTEGER 4 0 0 INTEGER ()) 0 0 () (1 0 ASSUMED_SHAPE
(CONSTANT (INTEGER 4 0 0 INTEGER ()) 0 '1') ()) 0 () () () 0 0)
497 'receives' '' 'receives' 492 ((VARIABLE IN UNKNOWN-PROC UNKNOWN
UNKNOWN 0 0 DIMENSION DUMMY) (INTEGER 4 0 0 INTEGER ()) 0 0 () (1 0
ASSUMED_SHAPE (CONSTANT (INTEGER 4 0 0 INTEGER ()) 0 '1') ()) 0 () () ()
0 0)
496 'send_starts' '' 'send_starts' 492 ((VARIABLE IN UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0 DIMENSION DUMMY) (INTEGER 4 0 0 INTEGER ()) 0 0 () (
1 0 ASSUMED_SHAPE (CONSTANT (INTEGER 4 0 0 INTEGER ()) 0 '1') ()) 0 () ()
() 0 0)
494 'nprocs' '' 'nprocs' 492 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN
0 0 DUMMY) (INTEGER 4 0 0 INTEGER ()) 0 0 () () 0 () () () 0 0)
493 'halo' '' 'halo' 492 ((VARIABLE INOUT UNKNOWN-PROC UNKNOWN UNKNOWN 0
0 DUMMY) (DERIVED 131 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
501 'create_caches' '' 'create_caches' 492 ((VARIABLE IN UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0 OPTIONAL DUMMY) (LOGICAL 4 0 0 LOGICAL ()) 0 0 () ()
0 () () () 0 0)
500 'ordering_scheme' '' 'ordering_scheme' 492 ((VARIABLE IN
UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 OPTIONAL DUMMY) (INTEGER 4 0 0 INTEGER
()) 0 0 () () 0 () () () 0 0)
506 'node_list' '' 'node_list' 503 ((VARIABLE UNKNOWN-INTENT
UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 DIMENSION DUMMY) (INTEGER 4 0 0 INTEGER
()) 0 0 () (1 0 ASSUMED_SHAPE (CONSTANT (INTEGER 4 0 0 INTEGER ()) 0 '1')
()) 0 () () () 0 0)
505 'subdomain_mesh' '' 'subdomain_mesh' 503 ((VARIABLE INOUT
UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 DUMMY) (DERIVED 105 0 0 DERIVED ()) 0 0
() () 0 () () () 0 0)
788 'mesh' '' 'mesh' 786 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DUMMY) (DERIVED 105 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
784 'positions' '' 'positions' 783 ((VARIABLE INOUT UNKNOWN-PROC UNKNOWN
UNKNOWN 0 0 DUMMY) (DERIVED 95 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
507 'inverse_node_list' '' 'inverse_node_list' 503 ((VARIABLE
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 DIMENSION DUMMY) (
INTEGER 4 0 0 INTEGER ()) 0 0 () (1 0 ASSUMED_SHAPE (CONSTANT (INTEGER 4
0 0 INTEGER ()) 0 '1') ()) 0 () () () 0 0)
499 'nowned_nodes' '' 'nowned_nodes' 492 ((VARIABLE IN UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0 OPTIONAL DUMMY) (INTEGER 4 0 0 INTEGER ()) 0 0 () ()
0 () () () 0 0)
498 'receive_starts' '' 'receive_starts' 492 ((VARIABLE IN UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0 DIMENSION DUMMY) (INTEGER 4 0 0 INTEGER ()) 0 0 () (
1 0 ASSUMED_SHAPE (CONSTANT (INTEGER 4 0 0 INTEGER ()) 0 '1') ()) 0 () ()
() 0 0)
488 'receives' '' 'receives' 484 ((VARIABLE OUT UNKNOWN-PROC UNKNOWN
UNKNOWN 0 0 DIMENSION DUMMY) (INTEGER 4 0 0 INTEGER ()) 0 0 () (1 0
ASSUMED_SHAPE (CONSTANT (INTEGER 4 0 0 INTEGER ()) 0 '1') ()) 0 () () ()
0 0)
400 'mesh' '' 'mesh' 398 ((VARIABLE INOUT UNKNOWN-PROC UNKNOWN UNKNOWN 0
0 DUMMY) (DERIVED 105 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
401 'communicator' '' 'communicator' 398 ((VARIABLE IN UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0 OPTIONAL DUMMY) (INTEGER 4 0 0 INTEGER ()) 0 0 () ()
0 () () () 0 0)
408 'positions' '' 'positions' 406 ((VARIABLE INOUT UNKNOWN-PROC UNKNOWN
UNKNOWN 0 0 DUMMY) (DERIVED 95 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
407 'filename' '' 'filename' 406 ((VARIABLE IN UNKNOWN-PROC UNKNOWN
UNKNOWN 0 0 DUMMY) (CHARACTER 1 0 0 CHARACTER (())) 0 0 () () 0 () () ()
0 0)
409 'communicator' '' 'communicator' 406 ((VARIABLE IN UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0 OPTIONAL DUMMY) (INTEGER 4 0 0 INTEGER ()) 0 0 () ()
0 () () () 0 0)
399 'filename' '' 'filename' 398 ((VARIABLE IN UNKNOWN-PROC UNKNOWN
UNKNOWN 0 0 DUMMY) (CHARACTER 1 0 0 CHARACTER (())) 0 0 () () 0 () () ()
0 0)
485 'halo' '' 'halo' 484 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DUMMY) (DERIVED 131 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
486 'sends' '' 'sends' 484 ((VARIABLE OUT UNKNOWN-PROC UNKNOWN UNKNOWN 0
0 DIMENSION DUMMY) (INTEGER 4 0 0 INTEGER ()) 0 0 () (1 0 ASSUMED_SHAPE
(CONSTANT (INTEGER 4 0 0 INTEGER ()) 0 '1') ()) 0 () () () 0 0)
487 'send_starts' '' 'send_starts' 484 ((VARIABLE OUT UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0 DIMENSION DUMMY) (INTEGER 4 0 0 INTEGER ()) 0 0 () (
1 0 ASSUMED_SHAPE (CONSTANT (INTEGER 4 0 0 INTEGER ()) 0 '1') ()) 0 () ()
() 0 0)
504 'external_mesh' '' 'external_mesh' 503 ((VARIABLE IN UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0 DUMMY) (DERIVED 105 0 0 DERIVED ()) 0 0 () () 0 () ()
() 0 0)
787 'filename' '' 'filename' 786 ((VARIABLE IN UNKNOWN-PROC UNKNOWN
UNKNOWN 0 0 DUMMY) (CHARACTER 1 0 0 CHARACTER (())) 0 0 () () 0 () () ()
0 0)
687 'mesh' '' 'mesh' 684 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DUMMY) (DERIVED 105 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
686 'node_halo' '' 'node_halo' 684 ((VARIABLE IN UNKNOWN-PROC UNKNOWN
UNKNOWN 0 0 DUMMY) (DERIVED 131 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
692 'mesh' '' 'mesh' 689 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DUMMY) (DERIVED 105 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
420 'halo' '' 'halo' 419 ((VARIABLE INOUT UNKNOWN-PROC UNKNOWN UNKNOWN 0
0 DUMMY) (DERIVED 131 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
701 'sorted_l1_halo' '' 'sorted_l1_halo' 698 ((VARIABLE IN UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0 OPTIONAL DUMMY) (LOGICAL 4 0 0 LOGICAL ()) 0 0 () ()
0 () () () 0 0)
418 'repair_halo' '' 'repair_halo' 416 ((VARIABLE IN UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0 DUMMY) (DERIVED 131 0 0 DERIVED ()) 0 0 () () 0 () ()
() 0 0)
417 'halo' '' 'halo' 416 ((VARIABLE INOUT UNKNOWN-PROC UNKNOWN UNKNOWN 0
0 DUMMY) (DERIVED 131 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
421 'repair_field' '' 'repair_field' 419 ((VARIABLE IN UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0 DUMMY) (DERIVED 95 0 0 DERIVED ()) 0 0 () () 0 () ()
() 0 0)
691 'element_halo' '' 'element_halo' 689 ((VARIABLE IN UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0 DUMMY) (DERIVED 131 0 0 DERIVED ()) 0 0 () () 0 () ()
() 0 0)
690 'node_halo' '' 'node_halo' 689 ((VARIABLE INOUT UNKNOWN-PROC UNKNOWN
UNKNOWN 0 0 DUMMY) (DERIVED 131 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
685 'element_halo' '' 'element_halo' 684 ((VARIABLE INOUT UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0 DUMMY) (DERIVED 131 0 0 DERIVED ()) 0 0 () () 0 () ()
() 0 0)
695 'halo' '' 'halo' 694 ((VARIABLE INOUT UNKNOWN-PROC UNKNOWN UNKNOWN 0
0 DUMMY) (DERIVED 131 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
700 'l2_halo' '' 'l2_halo' 698 ((VARIABLE IN UNKNOWN-PROC UNKNOWN
UNKNOWN 0 0 DUMMY) (DERIVED 131 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
699 'l1_halo' '' 'l1_halo' 698 ((VARIABLE INOUT UNKNOWN-PROC UNKNOWN
UNKNOWN 0 0 DUMMY) (DERIVED 131 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
696 'repair_field' '' 'repair_field' 694 ((VARIABLE IN UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0 DUMMY) (DERIVED 95 0 0 DERIVED ()) 0 0 () () 0 () ()
() 0 0)
795 'node_count' 'halos_base' 'node_count' 1 ((PROCEDURE UNKNOWN-INTENT
UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 GENERIC) (UNKNOWN 0 0 0 UNKNOWN ()) 0 0
() () 0 () () () 0 0)
)

('create_global_to_universal_numbering' 0 427 'create_ownership' 0 431
'deallocate_ownership_cache' 0 434 'deallocate_universal_numbering_cache'
0 437 'derive_element_halo_from_node_halo' 0 440
'derive_maximal_surface_element_halo' 0 445
'derive_nonperiodic_halos_from_periodic_halos' 0 452 'derive_sub_halo' 0
457 'ele_owner' 0 462 'ewrite_universal_numbers' 0 467
'extract_all_halo_receives' 0 471 'extract_all_halo_sends' 0 477
'extract_raw_halo_data' 0 483 'form_halo_from_raw_data' 0 491
'generate_substate_halos' 0 502 'get_node_owners' 0 508 'get_owned_nodes'
0 512 'get_universal_numbering' 0 6 'get_universal_numbering_inverse' 0
519 'halo_all_receives_count' 0 523 'halo_all_sends_count' 0 526
'halo_all_unique_receives_count' 0 529 'halo_data_type' 0 532 'halo_name'
0 535 'halo_node_owner' 0 540 'halo_node_owners' 0 546 'halo_nowned_nodes'
0 553 'halo_order_general' 0 556 'halo_order_trailing_receives' 0 557
'halo_ordering_scheme' 0 558 'halo_pointer' 0 561 'halo_proc_count' 0
384 'halo_receive' 0 563 'halo_receive_count' 0 568 'halo_receive_counts'
0 572 'halo_receives' 0 576 'halo_send' 0 580 'halo_send_count' 0 585
'halo_send_counts' 0 589 'halo_sends' 0 593 'halo_type' 0 131
'halo_type_cg_node' 0 611 'halo_type_dg_node' 0 612 'halo_type_element'
0 613 'halo_unique_receive_count' 0 614 'halo_universal_node_owners' 0
618 'halo_universal_number' 0 45 'halo_universal_numbers' 0 628
'halo_valid_for_communication' 0 634 'halos' 0 638
'has_global_to_universal_numbering' 0 639 'has_nowned_nodes' 0 643
'has_ownership' 0 646 'invert_comms_sizes' 0 649 'max_halo_node' 0 655
'max_halo_receive_node' 0 659 'max_halo_send_node' 0 663 'min_halo_node'
0 667 'min_halo_receive_node' 0 671 'min_halo_send_node' 0 675
'print_halo' 0 679 'reorder_element_halo' 0 683
'reorder_halo_from_element_halo' 0 688 'reorder_halo_receives' 0 693
'reorder_l1_from_l2_halo' 0 697 'set_all_halo_receives' 0 702
'set_all_halo_sends' 0 706 'set_halo_communicator' 0 710
'set_halo_data_type' 0 714 'set_halo_name' 0 718 'set_halo_nowned_nodes'
0 722 'set_halo_ordering_scheme' 0 726 'set_halo_receive' 0 730
'set_halo_receives' 0 736 'set_halo_send' 0 741 'set_halo_sends' 0 747
'set_halo_universal_number' 0 752 'trailing_receives_consistent' 0 758
'universal_numbering_count' 0 762 'valid_global_to_universal_numbering'
0 766 'valid_halo_communicator' 0 770 'valid_halo_node_counts' 0 774
'valid_serial_halo' 0 778 'verify_halos' 0 782 'write_halos' 0 785
'write_universal_numbering' 0 789)
