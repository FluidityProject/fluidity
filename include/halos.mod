GFORTRAN module version '6' created from Halos.F90 on Fri Nov  2 16:08:29 2012
MD5:5feef56533d908e4903b236df9f5255f -- If you edit this, you'll get what you deserve.

(() () () () () () () () () () () () () () () () () () () () () () ()
() () () ())

()

(('allocate' 'halos_allocates' 2 3) ('deallocate' 'halos_allocates' 4 5)
('get_universal_numbering' 'halos_numbering' 6 7) ('halo_communicator'
'halos_base' 8) ('halo_max' 'halos_communications' 9 10 11) (
'halo_universal_number' 'halos_numbering' 12 13) ('halo_verifies'
'halos_communications' 14 15 16 17 18) ('halo_update'
'halos_communications' 19 20 21 22 23 24 25 26 27 28 29 30 31) ('incref'
'halos_allocates' 32) ('node_owned' 'halos_ownership' 33) ('node_count'
'halos_base' 34) ('pending_communication' 'halos_debug' 35 36) ('nullify'
'halos_allocates' 37) ('read_halos' 'halos_registration' 38 39) (
'reallocate' 'halos_allocates' 40) ('nodes_owned' 'halos_ownership' 41)
('reorder_halo' 'halos_repair' 42 43) ('has_references' 'halos_allocates'
44) ('serial_storage_halo' 'halos_base' 45 46) ('derive_l1_from_l2_halo'
'halos_derivation' 47) ('zero' 'halos_base' 48))

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
2 'allocate_halo_halo' 'halos_allocates' 'allocate_halo_halo' 1 ((
PROCEDURE UNKNOWN-INTENT MODULE-PROC DECL UNKNOWN 0 0 SUBROUTINE) (
UNKNOWN 0 0 0 UNKNOWN ()) 66 0 (67 68) () 0 () () () 0 0)
5 'deallocate_halo' 'halos_allocates' 'deallocate_halo' 1 ((PROCEDURE
UNKNOWN-INTENT MODULE-PROC DECL UNKNOWN 0 0 SUBROUTINE) (UNKNOWN 0 0 0
UNKNOWN ()) 69 0 (70) () 0 () () () 0 0)
32 'incref_halo_type' 'halos_allocates' 'incref_halo_type' 1 ((
PROCEDURE UNKNOWN-INTENT MODULE-PROC DECL UNKNOWN 0 0 SUBROUTINE
ALWAYS_EXPLICIT) (UNKNOWN 0 0 0 UNKNOWN ()) 71 0 (72) () 0 () () () 0 0)
44 'has_references_halo_type' 'halos_allocates' 'has_references_halo_type'
1 ((PROCEDURE UNKNOWN-INTENT MODULE-PROC DECL UNKNOWN 0 0 FUNCTION PURE)
(LOGICAL 4 0 0 LOGICAL ()) 73 0 (74) () 75 () () () 0 0)
76 'refcount_type' 'reference_counting' 'refcount_type' 1 ((DERIVED
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 POINTER_COMP) (UNKNOWN 0
0 0 UNKNOWN ()) 0 0 () () 0 ((77 'prev' (DERIVED 76 0 0 DERIVED ()) () (
UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 POINTER)
UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0)) (78 'next' (DERIVED
76 0 0 DERIVED ()) () (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN
UNKNOWN 0 0 POINTER) UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0))
(79 'count' (INTEGER 4 0 0 INTEGER ()) () (UNKNOWN-FL UNKNOWN-INTENT
UNKNOWN-PROC UNKNOWN UNKNOWN 0 0) UNKNOWN-ACCESS (CONSTANT (INTEGER 4 0
0 INTEGER ()) 0 '0')) (80 'id' (INTEGER 4 0 0 INTEGER ()) () (
UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0)
UNKNOWN-ACCESS ()) (81 'name' (CHARACTER 1 0 0 CHARACTER ((CONSTANT (
INTEGER 4 0 0 INTEGER ()) 0 '101'))) () (UNKNOWN-FL UNKNOWN-INTENT
UNKNOWN-PROC UNKNOWN UNKNOWN 0 0) UNKNOWN-ACCESS ()) (82 'type' (
CHARACTER 1 0 0 CHARACTER ((CONSTANT (INTEGER 4 0 0 INTEGER ()) 0 '101')))
() (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0)
UNKNOWN-ACCESS ()) (83 'tagged' (LOGICAL 4 0 0 LOGICAL ()) () (
UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0)
UNKNOWN-ACCESS (CONSTANT (LOGICAL 4 0 0 LOGICAL ()) 0 0))) PUBLIC (() ()
() ()) () 0 0 25948645)
84 'integer_vector' 'futils' 'integer_vector' 1 ((DERIVED UNKNOWN-INTENT
UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 POINTER_COMP) (UNKNOWN 0 0 0 UNKNOWN ())
0 0 () () 0 ((85 'ptr' (INTEGER 4 0 0 INTEGER ()) (1 0 DEFERRED () ()) (
UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 DIMENSION
POINTER) UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0))) PUBLIC (()
() () ()) () 0 0 9661976)
40 'reallocate_halo' 'halos_allocates' 'reallocate_halo' 1 ((PROCEDURE
UNKNOWN-INTENT MODULE-PROC DECL UNKNOWN 0 0 SUBROUTINE ALWAYS_EXPLICIT)
(UNKNOWN 0 0 0 UNKNOWN ()) 86 0 (87 88 89) () 0 () () () 0 0)
51 'mpi_status_ignore' 'mpi_interfaces' 'mpi_status_ignore' 1 ((
VARIABLE UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 DIMENSION
IN_COMMON) (INTEGER 4 0 0 INTEGER ()) 0 0 () (1 0 EXPLICIT (CONSTANT (
INTEGER 4 0 0 INTEGER ()) 0 '1') (CONSTANT (INTEGER 4 0 0 INTEGER ()) 0
'5')) 0 () () () 0 0)
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
37 'nullify_halo' 'halos_allocates' 'nullify_halo' 1 ((PROCEDURE
UNKNOWN-INTENT MODULE-PROC DECL UNKNOWN 0 0 SUBROUTINE) (UNKNOWN 0 0 0
UNKNOWN ()) 90 0 (91) () 0 () () () 0 0)
52 'mpi_statuses_ignore' 'mpi_interfaces' 'mpi_statuses_ignore' 1 ((
VARIABLE UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 IN_COMMON) (
REAL 8 0 0 REAL ()) 0 0 () () 0 () () () 0 0)
36 'pending_communication_communicator' 'parallel_tools'
'pending_communication_communicator' 1 ((PROCEDURE UNKNOWN-INTENT
MODULE-PROC DECL UNKNOWN 0 0 FUNCTION ALWAYS_EXPLICIT) (LOGICAL 4 0 0
LOGICAL ()) 92 0 (93) () 94 () () () 0 0)
35 'pending_communication_halo' 'halos_debug' 'pending_communication_halo'
1 ((PROCEDURE UNKNOWN-INTENT MODULE-PROC DECL UNKNOWN 0 0 FUNCTION) (
LOGICAL 4 0 0 LOGICAL ()) 95 0 (96) () 97 () () () 0 0)
3 'allocate_halo' 'halos_allocates' 'allocate_halo' 1 ((PROCEDURE
UNKNOWN-INTENT MODULE-PROC DECL UNKNOWN 0 0 SUBROUTINE ALWAYS_EXPLICIT)
(UNKNOWN 0 0 0 UNKNOWN ()) 98 0 (99 100 101 102 103 104 105 106 107) ()
0 () () () 0 0)
59 'petsc_null_real' 'petscsysdef' 'petsc_null_real' 1 ((VARIABLE
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 IN_COMMON) (REAL 8 0 0
REAL ()) 0 0 () () 0 () () () 0 0)
58 'petsc_null_double' 'petscsysdef' 'petsc_null_double' 1 ((VARIABLE
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 IN_COMMON) (REAL 8 0 0
REAL ()) 0 0 () () 0 () () () 0 0)
61 'petsc_null_truth' 'petscsysdef' 'petsc_null_truth' 1 ((VARIABLE
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 IN_COMMON) (LOGICAL 4 0
0 LOGICAL ()) 0 0 () () 0 () () () 0 0)
57 'petsc_null_scalar' 'petscsysdef' 'petsc_null_scalar' 1 ((VARIABLE
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 IN_COMMON) (REAL 8 0 0
REAL ()) 0 0 () () 0 () () () 0 0)
63 'petsc_comm_world' 'petscsysdef' 'petsc_comm_world' 1 ((VARIABLE
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 IN_COMMON) (INTEGER 4 0
0 INTEGER ()) 0 0 () () 0 () () () 0 0)
62 'petsc_null_object' 'petscsysdef' 'petsc_null_object' 1 ((VARIABLE
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 IN_COMMON) (INTEGER 8 0
0 INTEGER ()) 0 0 () () 0 () () () 0 0)
53 'petsc_null_character' 'petscsysdef' 'petsc_null_character' 1 ((
VARIABLE UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN IMPLICIT-SAVE 0 0 IN_COMMON)
(CHARACTER 1 0 0 CHARACTER ((CONSTANT (INTEGER 4 0 0 INTEGER ()) 0 '80')))
0 0 () () 0 () () () 0 0)
54 'petsc_comm_self' 'petscsysdef' 'petsc_comm_self' 1 ((VARIABLE
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 IN_COMMON) (INTEGER 4 0
0 INTEGER ()) 0 0 () () 0 () () () 0 0)
55 'petsc_null_integer' 'petscsysdef' 'petsc_null_integer' 1 ((VARIABLE
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 IN_COMMON) (INTEGER 4 0
0 INTEGER ()) 0 0 () () 0 () () () 0 0)
56 'petsc_null' 'petscsysdef' 'petsc_null' 1 ((VARIABLE UNKNOWN-INTENT
UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 IN_COMMON) (INTEGER 4 0 0 INTEGER ()) 0
0 () () 0 () () () 0 0)
108 'vector_field' 'fields_data_types' 'vector_field' 1 ((DERIVED
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 POINTER_COMP) (UNKNOWN 0
0 0 UNKNOWN ()) 0 0 () () 0 ((109 'val' (REAL 8 0 0 REAL ()) (2 0
DEFERRED () () () ()) (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN
UNKNOWN 0 0 DIMENSION POINTER) UNKNOWN-ACCESS ()) (110 'wrapped' (
LOGICAL 4 0 0 LOGICAL ()) () (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0) UNKNOWN-ACCESS (CONSTANT (LOGICAL 4 0 0 LOGICAL ())
0 1)) (111 'field_type' (INTEGER 4 0 0 INTEGER ()) () (UNKNOWN-FL
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0) UNKNOWN-ACCESS (
CONSTANT (INTEGER 4 0 0 INTEGER ()) 0 '0')) (112 'bc' (DERIVED 113 0 0
DERIVED ()) () (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0
0 POINTER) UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0)) (114 'name'
(CHARACTER 1 0 0 CHARACTER ((CONSTANT (INTEGER 4 0 0 INTEGER ()) 0 '101')))
() (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0)
UNKNOWN-ACCESS ()) (115 'dim' (INTEGER 4 0 0 INTEGER ()) () (UNKNOWN-FL
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0) UNKNOWN-ACCESS ()) (
116 'option_path' (CHARACTER 1 0 0 CHARACTER ((CONSTANT (INTEGER 4 0 0
INTEGER ()) 0 '8192'))) () (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0) UNKNOWN-ACCESS (CONSTANT (CHARACTER 1 0 0 CHARACTER
(())) 0 8192
'/uninitialised_path/                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            '))
(117 'mesh' (DERIVED 118 0 0 DERIVED ()) () (UNKNOWN-FL UNKNOWN-INTENT
UNKNOWN-PROC UNKNOWN UNKNOWN 0 0) UNKNOWN-ACCESS (STRUCTURE (DERIVED 118
0 0 DERIVED ()) 0 ((() ()) ((CONSTANT (LOGICAL 4 0 0 LOGICAL ()) 0 1) ())
((STRUCTURE (DERIVED 119 0 0 DERIVED ()) 0 ((() ()) (() ()) (() ()) (()
()) ((NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0) ()) ((NULL (UNKNOWN 0 0 0
UNKNOWN ()) 0) ()) ((NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0) ()) ((NULL (
UNKNOWN 0 0 0 UNKNOWN ()) 0) ()) ((NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0) ())
((NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0) ()) ((NULL (UNKNOWN 0 0 0 UNKNOWN ())
0) ()) ((STRUCTURE (DERIVED 120 0 0 DERIVED ()) 0 ((() ()) (() ()) (() ())
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
UNKNOWN ()) 0) ()) ((NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0) ()) ((CONSTANT (
LOGICAL 4 0 0 LOGICAL ()) 0 0) ())) ())) (121 'refcount' (DERIVED 76 0 0
DERIVED ()) () (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0
0 POINTER) UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0)) (122
'aliased' (LOGICAL 4 0 0 LOGICAL ()) () (UNKNOWN-FL UNKNOWN-INTENT
UNKNOWN-PROC UNKNOWN UNKNOWN 0 0) UNKNOWN-ACCESS (CONSTANT (LOGICAL 4 0
0 LOGICAL ()) 0 0)) (123 'picker' (DERIVED 124 0 0 DERIVED ()) () (
UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 POINTER)
UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0))) PUBLIC (() () () ())
() 0 0 43598963)
118 'mesh_type' 'fields_data_types' 'mesh_type' 1 ((DERIVED
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 POINTER_COMP) (UNKNOWN 0
0 0 UNKNOWN ()) 0 0 () () 0 ((125 'ndglno' (INTEGER 4 0 0 INTEGER ()) (
1 0 DEFERRED () ()) (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN
UNKNOWN 0 0 DIMENSION POINTER) UNKNOWN-ACCESS ()) (126 'wrapped' (
LOGICAL 4 0 0 LOGICAL ()) () (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0) UNKNOWN-ACCESS (CONSTANT (LOGICAL 4 0 0 LOGICAL ())
0 1)) (127 'shape' (DERIVED 119 0 0 DERIVED ()) () (UNKNOWN-FL
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0) UNKNOWN-ACCESS (
STRUCTURE (DERIVED 119 0 0 DERIVED ()) 0 ((() ()) (() ()) (() ()) (() ())
((NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0) ()) ((NULL (UNKNOWN 0 0 0 UNKNOWN ())
0) ()) ((NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0) ()) ((NULL (UNKNOWN 0 0 0
UNKNOWN ()) 0) ()) ((NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0) ()) ((NULL (
UNKNOWN 0 0 0 UNKNOWN ()) 0) ()) ((NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0) ())
((STRUCTURE (DERIVED 120 0 0 DERIVED ()) 0 ((() ()) (() ()) (() ()) (()
()) ((NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0) ()) ((NULL (UNKNOWN 0 0 0
UNKNOWN ()) 0) ()) (() ()) ((NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0) ()) (()
())) ()) ()) ((NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0) ()) ((NULL (UNKNOWN 0
0 0 UNKNOWN ()) 0) ()) ((NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0) ()) ((NULL (
UNKNOWN 0 0 0 UNKNOWN ()) 0) ()) (() ())) ())) (128 'elements' (INTEGER
4 0 0 INTEGER ()) () (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN
UNKNOWN 0 0) UNKNOWN-ACCESS ()) (129 'nodes' (INTEGER 4 0 0 INTEGER ())
() (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0)
UNKNOWN-ACCESS ()) (130 'name' (CHARACTER 1 0 0 CHARACTER ((CONSTANT (
INTEGER 4 0 0 INTEGER ()) 0 '101'))) () (UNKNOWN-FL UNKNOWN-INTENT
UNKNOWN-PROC UNKNOWN UNKNOWN 0 0) UNKNOWN-ACCESS ()) (131 'option_path'
(CHARACTER 1 0 0 CHARACTER ((CONSTANT (INTEGER 4 0 0 INTEGER ()) 0 '8192')))
() (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0)
UNKNOWN-ACCESS (CONSTANT (CHARACTER 1 0 0 CHARACTER (())) 0 8192
'/uninitialised_path/                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            '))
(132 'continuity' (INTEGER 4 0 0 INTEGER ()) () (UNKNOWN-FL
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0) UNKNOWN-ACCESS (
CONSTANT (INTEGER 4 0 0 INTEGER ()) 0 '0')) (133 'refcount' (DERIVED 76
0 0 DERIVED ()) () (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN
UNKNOWN 0 0 POINTER) UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0))
(134 'faces' (DERIVED 135 0 0 DERIVED ()) () (UNKNOWN-FL UNKNOWN-INTENT
UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 POINTER) UNKNOWN-ACCESS (NULL (UNKNOWN
0 0 0 UNKNOWN ()) 0)) (136 'subdomain_mesh' (DERIVED 137 0 0 DERIVED ())
() (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 POINTER)
UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0)) (138 'adj_lists' (
DERIVED 139 0 0 DERIVED ()) () (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0 POINTER) UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN
()) 0)) (140 'columns' (INTEGER 4 0 0 INTEGER ()) (1 0 DEFERRED () ()) (
UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 DIMENSION
POINTER) UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0)) (141
'element_columns' (INTEGER 4 0 0 INTEGER ()) (1 0 DEFERRED () ()) (
UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 DIMENSION
POINTER) UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0)) (142
'region_ids' (INTEGER 4 0 0 INTEGER ()) (1 0 DEFERRED () ()) (
UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 DIMENSION
POINTER) UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0)) (143 'halos'
(DERIVED 144 0 0 DERIVED ()) (1 0 DEFERRED () ()) (UNKNOWN-FL
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 DIMENSION POINTER)
UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0)) (145 'element_halos'
(DERIVED 144 0 0 DERIVED ()) (1 0 DEFERRED () ()) (UNKNOWN-FL
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 DIMENSION POINTER)
UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0)) (146 'periodic' (
LOGICAL 4 0 0 LOGICAL ()) () (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0) UNKNOWN-ACCESS (CONSTANT (LOGICAL 4 0 0 LOGICAL ())
0 0))) PUBLIC (() () () ()) () 0 0 2572791)
113 'vector_boundary_conditions_ptr' 'fields_data_types'
'vector_boundary_conditions_ptr' 1 ((DERIVED UNKNOWN-INTENT UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0 POINTER_COMP) (UNKNOWN 0 0 0 UNKNOWN ()) 0 0 () () 0
((147 'boundary_condition' (DERIVED 148 0 0 DERIVED ()) (1 0 DEFERRED ()
()) (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DIMENSION POINTER) UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0)))
PUBLIC (() () () ()) () 0 0 84476181)
120 'quadrature_type' 'quadrature' 'quadrature_type' 1 ((DERIVED
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 POINTER_COMP) (UNKNOWN 0
0 0 UNKNOWN ()) 0 0 () () 0 ((149 'dim' (INTEGER 4 0 0 INTEGER ()) () (
UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0)
UNKNOWN-ACCESS ()) (150 'degree' (INTEGER 4 0 0 INTEGER ()) () (
UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0)
UNKNOWN-ACCESS ()) (151 'vertices' (INTEGER 4 0 0 INTEGER ()) () (
UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0)
UNKNOWN-ACCESS ()) (152 'ngi' (INTEGER 4 0 0 INTEGER ()) () (UNKNOWN-FL
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0) UNKNOWN-ACCESS ()) (
153 'weight' (REAL 8 0 0 REAL ()) (1 0 DEFERRED () ()) (UNKNOWN-FL
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 DIMENSION POINTER)
UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0)) (154 'l' (REAL 8 0 0
REAL ()) (2 0 DEFERRED () () () ()) (UNKNOWN-FL UNKNOWN-INTENT
UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 DIMENSION POINTER) UNKNOWN-ACCESS (
NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0)) (155 'name' (CHARACTER 1 0 0
CHARACTER ((CONSTANT (INTEGER 4 0 0 INTEGER ()) 0 '0'))) () (UNKNOWN-FL
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0) UNKNOWN-ACCESS ()) (
156 'refcount' (DERIVED 76 0 0 DERIVED ()) () (UNKNOWN-FL UNKNOWN-INTENT
UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 POINTER) UNKNOWN-ACCESS (NULL (UNKNOWN
0 0 0 UNKNOWN ()) 0)) (157 'family' (INTEGER 4 0 0 INTEGER ()) () (
UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0)
UNKNOWN-ACCESS ())) PUBLIC (() () () ()) () 0 0 59837722)
119 'element_type' 'elements' 'element_type' 1 ((DERIVED UNKNOWN-INTENT
UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 POINTER_COMP) (UNKNOWN 0 0 0 UNKNOWN ())
0 0 () () 0 ((158 'dim' (INTEGER 4 0 0 INTEGER ()) () (UNKNOWN-FL
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0) UNKNOWN-ACCESS ()) (
159 'loc' (INTEGER 4 0 0 INTEGER ()) () (UNKNOWN-FL UNKNOWN-INTENT
UNKNOWN-PROC UNKNOWN UNKNOWN 0 0) UNKNOWN-ACCESS ()) (160 'ngi' (
INTEGER 4 0 0 INTEGER ()) () (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0) UNKNOWN-ACCESS ()) (161 'degree' (INTEGER 4 0 0
INTEGER ()) () (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0
0) UNKNOWN-ACCESS ()) (162 'n' (REAL 8 0 0 REAL ()) (2 0 DEFERRED () ()
() ()) (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DIMENSION POINTER) UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0)) (
163 'dn' (REAL 8 0 0 REAL ()) (3 0 DEFERRED () () () () () ()) (
UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 DIMENSION
POINTER) UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0)) (164 'n_s'
(REAL 8 0 0 REAL ()) (3 0 DEFERRED () () () () () ()) (UNKNOWN-FL
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 DIMENSION POINTER)
UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0)) (165 'dn_s' (REAL 8
0 0 REAL ()) (4 0 DEFERRED () () () () () () () ()) (UNKNOWN-FL
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 DIMENSION POINTER)
UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0)) (166 'spoly' (
DERIVED 167 0 0 DERIVED ()) (2 0 DEFERRED () () () ()) (UNKNOWN-FL
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 DIMENSION POINTER)
UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0)) (168 'dspoly' (
DERIVED 167 0 0 DERIVED ()) (2 0 DEFERRED () () () ()) (UNKNOWN-FL
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 DIMENSION POINTER)
UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0)) (169 'numbering' (
DERIVED 170 0 0 DERIVED ()) () (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0 POINTER) UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN
()) 0)) (171 'quadrature' (DERIVED 120 0 0 DERIVED ()) () (UNKNOWN-FL
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0) UNKNOWN-ACCESS (
STRUCTURE (DERIVED 120 0 0 DERIVED ()) 0 ((() ()) (() ()) (() ()) (() ())
((NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0) ()) ((NULL (UNKNOWN 0 0 0 UNKNOWN ())
0) ()) (() ()) ((NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0) ()) (() ())) ())) (
172 'surface_quadrature' (DERIVED 120 0 0 DERIVED ()) () (UNKNOWN-FL
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 POINTER) UNKNOWN-ACCESS
(NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0)) (173 'superconvergence' (DERIVED
174 0 0 DERIVED ()) () (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN
UNKNOWN 0 0 POINTER) UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0))
(175 'constraints' (DERIVED 176 0 0 DERIVED ()) () (UNKNOWN-FL
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 POINTER) UNKNOWN-ACCESS
(NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0)) (177 'refcount' (DERIVED 76 0 0
DERIVED ()) () (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0
0 POINTER) UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0)) (178 'name'
(CHARACTER 1 0 0 CHARACTER ((CONSTANT (INTEGER 4 0 0 INTEGER ()) 0 '0')))
() (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0)
UNKNOWN-ACCESS ())) PUBLIC (() () () ()) () 0 0 79461029)
135 'mesh_faces' 'fields_data_types' 'mesh_faces' 1 ((DERIVED
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 POINTER_COMP) (UNKNOWN 0
0 0 UNKNOWN ()) 0 0 () () 0 ((179 'shape' (DERIVED 119 0 0 DERIVED ()) ()
(UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 POINTER)
UNKNOWN-ACCESS ()) (180 'face_list' (DERIVED 181 0 0 DERIVED ()) () (
UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0)
UNKNOWN-ACCESS (STRUCTURE (DERIVED 181 0 0 DERIVED ()) 0 (((STRUCTURE (
DERIVED 182 0 0 DERIVED ()) 0 (((NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0) ())
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
())) ())) (183 'face_lno' (INTEGER 4 0 0 INTEGER ()) (1 0 DEFERRED () ())
(UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 DIMENSION
POINTER) UNKNOWN-ACCESS ()) (184 'surface_mesh' (DERIVED 118 0 0 DERIVED
()) () (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0)
UNKNOWN-ACCESS (STRUCTURE (DERIVED 118 0 0 DERIVED ()) 0 ((() ()) ((
CONSTANT (LOGICAL 4 0 0 LOGICAL ()) 0 1) ()) ((STRUCTURE (DERIVED 119 0
0 DERIVED ()) 0 ((() ()) (() ()) (() ()) (() ()) ((NULL (UNKNOWN 0 0 0
UNKNOWN ()) 0) ()) ((NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0) ()) ((NULL (
UNKNOWN 0 0 0 UNKNOWN ()) 0) ()) ((NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0) ())
((NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0) ()) ((NULL (UNKNOWN 0 0 0 UNKNOWN ())
0) ()) ((NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0) ()) ((STRUCTURE (DERIVED 120
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
LOGICAL 4 0 0 LOGICAL ()) 0 0) ())) ())) (185 'surface_node_list' (
INTEGER 4 0 0 INTEGER ()) (1 0 DEFERRED () ()) (UNKNOWN-FL
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 DIMENSION POINTER)
UNKNOWN-ACCESS ()) (186 'face_element_list' (INTEGER 4 0 0 INTEGER ()) (
1 0 DEFERRED () ()) (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN
UNKNOWN 0 0 DIMENSION POINTER) UNKNOWN-ACCESS ()) (187 'boundary_ids' (
INTEGER 4 0 0 INTEGER ()) (1 0 DEFERRED () ()) (UNKNOWN-FL
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 DIMENSION POINTER)
UNKNOWN-ACCESS ()) (188 'coplanar_ids' (INTEGER 4 0 0 INTEGER ()) (1 0
DEFERRED () ()) (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN
0 0 DIMENSION POINTER) UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0))
(189 'dg_surface_mesh' (DERIVED 118 0 0 DERIVED ()) () (UNKNOWN-FL
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 POINTER) UNKNOWN-ACCESS
(NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0)) (190 'has_internal_boundaries' (
LOGICAL 4 0 0 LOGICAL ()) () (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0) UNKNOWN-ACCESS (CONSTANT (LOGICAL 4 0 0 LOGICAL ())
0 0))) PUBLIC (() () () ()) () 0 0 11936185)
124 'picker_ptr' 'picker_data_types' 'picker_ptr' 1 ((DERIVED
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 POINTER_COMP) (UNKNOWN 0
0 0 UNKNOWN ()) 0 0 () () 0 ((191 'ptr' (DERIVED 192 0 0 DERIVED ()) ()
(UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 POINTER)
UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0))) PUBLIC (() () () ())
() 0 0 71120007)
139 'adjacency_cache' 'fields_data_types' 'adjacency_cache' 1 ((DERIVED
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 POINTER_COMP) (UNKNOWN 0
0 0 UNKNOWN ()) 0 0 () () 0 ((193 'nnlist' (DERIVED 182 0 0 DERIVED ())
() (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 POINTER)
UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0)) (194 'nelist' (
DERIVED 182 0 0 DERIVED ()) () (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0 POINTER) UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN
()) 0)) (195 'eelist' (DERIVED 182 0 0 DERIVED ()) () (UNKNOWN-FL
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 POINTER) UNKNOWN-ACCESS
(NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0))) PUBLIC (() () () ()) () 0 0
37419158)
137 'mesh_subdomain_mesh' 'fields_data_types' 'mesh_subdomain_mesh' 1 (
(DERIVED UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 POINTER_COMP) (
UNKNOWN 0 0 0 UNKNOWN ()) 0 0 () () 0 ((196 'element_list' (INTEGER 4 0
0 INTEGER ()) (1 0 DEFERRED () ()) (UNKNOWN-FL UNKNOWN-INTENT
UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 DIMENSION POINTER) UNKNOWN-ACCESS ()) (
197 'node_list' (INTEGER 4 0 0 INTEGER ()) (1 0 DEFERRED () ()) (
UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 DIMENSION
POINTER) UNKNOWN-ACCESS ())) PUBLIC (() () () ()) () 0 0 85302181)
167 'polynomial' 'polynomials' 'polynomial' 1 ((DERIVED UNKNOWN-INTENT
UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 POINTER_COMP) (UNKNOWN 0 0 0 UNKNOWN ())
0 0 () () 0 ((198 'coefs' (REAL 8 0 0 REAL ()) (1 0 DEFERRED () ()) (
UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 DIMENSION
POINTER) UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0)) (199 'degree'
(INTEGER 4 0 0 INTEGER ()) () (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0) UNKNOWN-ACCESS (CONSTANT (INTEGER 4 0 0 INTEGER ())
0 '-1'))) PUBLIC (() () () ()) () 0 0 87989236)
170 'ele_numbering_type' 'element_numbering' 'ele_numbering_type' 1 ((
DERIVED UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 POINTER_COMP) (
UNKNOWN 0 0 0 UNKNOWN ()) 0 0 () () 0 ((200 'faces' (INTEGER 4 0 0
INTEGER ()) () (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0
0) UNKNOWN-ACCESS ()) (201 'vertices' (INTEGER 4 0 0 INTEGER ()) () (
UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0)
UNKNOWN-ACCESS ()) (202 'edges' (INTEGER 4 0 0 INTEGER ()) () (
UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0)
UNKNOWN-ACCESS ()) (203 'boundaries' (INTEGER 4 0 0 INTEGER ()) () (
UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0)
UNKNOWN-ACCESS ()) (204 'degree' (INTEGER 4 0 0 INTEGER ()) () (
UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0)
UNKNOWN-ACCESS ()) (205 'dimension' (INTEGER 4 0 0 INTEGER ()) () (
UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0)
UNKNOWN-ACCESS ()) (206 'nodes' (INTEGER 4 0 0 INTEGER ()) () (
UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0)
UNKNOWN-ACCESS ()) (207 'type' (INTEGER 4 0 0 INTEGER ()) () (
UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0)
UNKNOWN-ACCESS (CONSTANT (INTEGER 4 0 0 INTEGER ()) 0 '1')) (208 'family'
(INTEGER 4 0 0 INTEGER ()) () (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0) UNKNOWN-ACCESS ()) (209 'count2number' (INTEGER 4 0
0 INTEGER ()) (3 0 DEFERRED () () () () () ()) (UNKNOWN-FL
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 DIMENSION POINTER)
UNKNOWN-ACCESS ()) (210 'number2count' (INTEGER 4 0 0 INTEGER ()) (2 0
DEFERRED () () () ()) (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN
UNKNOWN 0 0 DIMENSION POINTER) UNKNOWN-ACCESS ()) (211 'boundary_coord'
(INTEGER 4 0 0 INTEGER ()) (1 0 DEFERRED () ()) (UNKNOWN-FL
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 DIMENSION POINTER)
UNKNOWN-ACCESS ()) (212 'boundary_val' (INTEGER 4 0 0 INTEGER ()) (1 0
DEFERRED () ()) (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN
0 0 DIMENSION POINTER) UNKNOWN-ACCESS ())) PUBLIC (() () () ()) () 0 0
96431082)
174 'superconvergence_type' 'elements' 'superconvergence_type' 1 ((
DERIVED UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 POINTER_COMP) (
UNKNOWN 0 0 0 UNKNOWN ()) 0 0 () () 0 ((213 'nsp' (INTEGER 4 0 0 INTEGER
()) () (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0)
UNKNOWN-ACCESS ()) (214 'l' (REAL 8 0 0 REAL ()) (2 0 DEFERRED () () ()
()) (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DIMENSION POINTER) UNKNOWN-ACCESS ()) (215 'n' (REAL 8 0 0 REAL ()) (2 0
DEFERRED () () () ()) (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN
UNKNOWN 0 0 DIMENSION POINTER) UNKNOWN-ACCESS ()) (216 'dn' (REAL 8 0 0
REAL ()) (3 0 DEFERRED () () () () () ()) (UNKNOWN-FL UNKNOWN-INTENT
UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 DIMENSION POINTER) UNKNOWN-ACCESS ()))
PUBLIC (() () () ()) () 0 0 18282395)
176 'constraints_type' 'elements' 'constraints_type' 1 ((DERIVED
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 POINTER_COMP) (UNKNOWN 0
0 0 UNKNOWN ()) 0 0 () () 0 ((217 'type' (INTEGER 4 0 0 INTEGER ()) () (
UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0)
UNKNOWN-ACCESS ()) (218 'dim' (INTEGER 4 0 0 INTEGER ()) () (UNKNOWN-FL
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0) UNKNOWN-ACCESS ()) (
219 'degree' (INTEGER 4 0 0 INTEGER ()) () (UNKNOWN-FL UNKNOWN-INTENT
UNKNOWN-PROC UNKNOWN UNKNOWN 0 0) UNKNOWN-ACCESS ()) (220 'loc' (
INTEGER 4 0 0 INTEGER ()) () (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0) UNKNOWN-ACCESS ()) (221 'n_constraints' (INTEGER 4
0 0 INTEGER ()) () (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN
UNKNOWN 0 0) UNKNOWN-ACCESS ()) (222 'orthogonal' (REAL 8 0 0 REAL ()) (
3 0 DEFERRED () () () () () ()) (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0 DIMENSION POINTER) UNKNOWN-ACCESS (NULL (UNKNOWN 0 0
0 UNKNOWN ()) 0))) PUBLIC (() () () ()) () 0 0 57548971)
181 'csr_matrix' 'sparse_tools' 'csr_matrix' 1 ((DERIVED UNKNOWN-INTENT
UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 POINTER_COMP) (UNKNOWN 0 0 0 UNKNOWN ())
0 0 () () 0 ((223 'sparsity' (DERIVED 182 0 0 DERIVED ()) () (
UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0)
UNKNOWN-ACCESS (STRUCTURE (DERIVED 182 0 0 DERIVED ()) 0 (((NULL (
UNKNOWN 0 0 0 UNKNOWN ()) 0) ()) ((NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0) ())
((NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0) ()) (() ()) ((NULL (UNKNOWN 0 0 0
UNKNOWN ()) 0) ()) ((NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0) ()) ((NULL (
UNKNOWN 0 0 0 UNKNOWN ()) 0) ()) ((CONSTANT (CHARACTER 1 0 0 CHARACTER (
())) 0 101
'                                                                                                     ')
()) ((CONSTANT (LOGICAL 4 0 0 LOGICAL ()) 0 0) ()) ((CONSTANT (LOGICAL 4
0 0 LOGICAL ()) 0 0) ())) ())) (224 'val' (REAL 8 0 0 REAL ()) (1 0
DEFERRED () ()) (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN
0 0 DIMENSION POINTER) UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0))
(225 'ival' (INTEGER 4 0 0 INTEGER ()) (1 0 DEFERRED () ()) (UNKNOWN-FL
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 DIMENSION POINTER)
UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0)) (226 'clone' (
LOGICAL 4 0 0 LOGICAL ()) () (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0) UNKNOWN-ACCESS (CONSTANT (LOGICAL 4 0 0 LOGICAL ())
0 0)) (227 'external_val' (LOGICAL 4 0 0 LOGICAL ()) () (UNKNOWN-FL
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0) UNKNOWN-ACCESS (
CONSTANT (LOGICAL 4 0 0 LOGICAL ()) 0 0)) (228 'inactive' (DERIVED 229 0
0 DERIVED ()) () (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN
0 0 POINTER) UNKNOWN-ACCESS ()) (230 'ksp' (INTEGER 8 0 0 INTEGER ()) ()
(UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 POINTER)
UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0)) (231 'refcount' (
DERIVED 76 0 0 DERIVED ()) () (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0 POINTER) UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN
()) 0)) (232 'name' (CHARACTER 1 0 0 CHARACTER ((CONSTANT (INTEGER 4 0 0
INTEGER ()) 0 '101'))) () (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0) UNKNOWN-ACCESS (CONSTANT (CHARACTER 1 0 0 CHARACTER
(())) 0 101
'                                                                                                     ')))
PUBLIC (() () () ()) () 0 0 66371681)
182 'csr_sparsity' 'sparse_tools' 'csr_sparsity' 1 ((DERIVED
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 POINTER_COMP) (UNKNOWN 0
0 0 UNKNOWN ()) 0 0 () () 0 ((233 'findrm' (INTEGER 4 0 0 INTEGER ()) (
1 0 DEFERRED () ()) (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN
UNKNOWN 0 0 DIMENSION POINTER) UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0
UNKNOWN ()) 0)) (234 'centrm' (INTEGER 4 0 0 INTEGER ()) (1 0 DEFERRED ()
()) (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DIMENSION POINTER) UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0)) (
235 'colm' (INTEGER 4 0 0 INTEGER ()) (1 0 DEFERRED () ()) (UNKNOWN-FL
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 DIMENSION POINTER)
UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0)) (236 'columns' (
INTEGER 4 0 0 INTEGER ()) () (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0) UNKNOWN-ACCESS ()) (237 'row_halo' (DERIVED 144 0 0
DERIVED ()) () (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0
0 POINTER) UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0)) (238
'column_halo' (DERIVED 144 0 0 DERIVED ()) () (UNKNOWN-FL UNKNOWN-INTENT
UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 POINTER) UNKNOWN-ACCESS (NULL (UNKNOWN
0 0 0 UNKNOWN ()) 0)) (239 'refcount' (DERIVED 76 0 0 DERIVED ()) () (
UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 POINTER)
UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0)) (240 'name' (
CHARACTER 1 0 0 CHARACTER ((CONSTANT (INTEGER 4 0 0 INTEGER ()) 0 '101')))
() (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0)
UNKNOWN-ACCESS (CONSTANT (CHARACTER 1 0 0 CHARACTER (())) 0 101
'                                                                                                     '))
(241 'wrapped' (LOGICAL 4 0 0 LOGICAL ()) () (UNKNOWN-FL UNKNOWN-INTENT
UNKNOWN-PROC UNKNOWN UNKNOWN 0 0) UNKNOWN-ACCESS (CONSTANT (LOGICAL 4 0
0 LOGICAL ()) 0 0)) (242 'sorted_rows' (LOGICAL 4 0 0 LOGICAL ()) () (
UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0)
UNKNOWN-ACCESS (CONSTANT (LOGICAL 4 0 0 LOGICAL ()) 0 0))) PUBLIC (() ()
() ()) () 0 0 78882345)
192 'picker_type' 'picker_data_types' 'picker_type' 1 ((DERIVED
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 POINTER_COMP) (UNKNOWN 0
0 0 UNKNOWN ()) 0 0 () () 0 ((243 'name' (CHARACTER 1 0 0 CHARACTER ((
CONSTANT (INTEGER 4 0 0 INTEGER ()) 0 '101'))) () (UNKNOWN-FL
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0) UNKNOWN-ACCESS ()) (
244 'refcount' (DERIVED 76 0 0 DERIVED ()) () (UNKNOWN-FL UNKNOWN-INTENT
UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 POINTER) UNKNOWN-ACCESS (NULL (UNKNOWN
0 0 0 UNKNOWN ()) 0)) (245 'picker_id' (INTEGER 4 0 0 INTEGER ()) () (
UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0)
UNKNOWN-ACCESS (CONSTANT (INTEGER 4 0 0 INTEGER ()) 0 '0')) (246
'last_mesh_movement' (INTEGER 4 0 0 INTEGER ()) () (UNKNOWN-FL
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0) UNKNOWN-ACCESS (
CONSTANT (INTEGER 4 0 0 INTEGER ()) 0 '0'))) PUBLIC (() () () ()) () 0 0
8821665)
229 'logical_array_ptr' 'sparse_tools' 'logical_array_ptr' 1 ((DERIVED
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 POINTER_COMP) (UNKNOWN 0
0 0 UNKNOWN ()) 0 0 () () 0 ((247 'ptr' (LOGICAL 4 0 0 LOGICAL ()) (1 0
DEFERRED () ()) (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN
0 0 DIMENSION POINTER) UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0)))
PUBLIC (() () () ()) () 0 0 30974511)
148 'vector_boundary_condition' 'fields_data_types'
'vector_boundary_condition' 1 ((DERIVED UNKNOWN-INTENT UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0 POINTER_COMP) (UNKNOWN 0 0 0 UNKNOWN ()) 0 0 () () 0
((248 'name' (CHARACTER 1 0 0 CHARACTER ((CONSTANT (INTEGER 4 0 0
INTEGER ()) 0 '101'))) () (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0) UNKNOWN-ACCESS ()) (249 'type' (CHARACTER 1 0 0
CHARACTER ((CONSTANT (INTEGER 4 0 0 INTEGER ()) 0 '101'))) () (
UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0)
UNKNOWN-ACCESS (CONSTANT (CHARACTER 1 0 0 CHARACTER (())) 0 101
'                                                                                                     '))
(250 'applies' (LOGICAL 4 0 0 LOGICAL ()) (1 0 EXPLICIT (CONSTANT (
INTEGER 4 0 0 INTEGER ()) 0 '1') (CONSTANT (INTEGER 4 0 0 INTEGER ()) 0
'3')) (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DIMENSION) UNKNOWN-ACCESS ()) (251 'surface_element_list' (INTEGER 4 0 0
INTEGER ()) (1 0 DEFERRED () ()) (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0 DIMENSION POINTER) UNKNOWN-ACCESS (NULL (UNKNOWN 0 0
0 UNKNOWN ()) 0)) (252 'surface_node_list' (INTEGER 4 0 0 INTEGER ()) (
1 0 DEFERRED () ()) (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN
UNKNOWN 0 0 DIMENSION POINTER) UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0
UNKNOWN ()) 0)) (253 'surface_mesh' (DERIVED 118 0 0 DERIVED ()) () (
UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 POINTER)
UNKNOWN-ACCESS ()) (254 'surface_fields' (DERIVED 108 0 0 DERIVED ()) (
1 0 DEFERRED () ()) (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN
UNKNOWN 0 0 DIMENSION POINTER) UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0
UNKNOWN ()) 0)) (255 'scalar_surface_fields' (DERIVED 256 0 0 DERIVED ())
(1 0 DEFERRED () ()) (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN
UNKNOWN 0 0 DIMENSION POINTER) UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0
UNKNOWN ()) 0)) (257 'option_path' (CHARACTER 1 0 0 CHARACTER ((
CONSTANT (INTEGER 4 0 0 INTEGER ()) 0 '8192'))) () (UNKNOWN-FL
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0) UNKNOWN-ACCESS (
CONSTANT (CHARACTER 1 0 0 CHARACTER (())) 0 8192
'/uninitialised_path/                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            ')))
PUBLIC (() () () ()) () 0 0 42205165)
256 'scalar_field' 'fields_data_types' 'scalar_field' 1 ((DERIVED
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 POINTER_COMP) (UNKNOWN 0
0 0 UNKNOWN ()) 0 0 () () 0 ((258 'val' (REAL 8 0 0 REAL ()) (1 0
DEFERRED () ()) (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN
0 0 DIMENSION POINTER) UNKNOWN-ACCESS ()) (259 'val_stride' (INTEGER 4 0
0 INTEGER ()) () (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN
0 0) UNKNOWN-ACCESS (CONSTANT (INTEGER 4 0 0 INTEGER ()) 0 '1')) (260
'wrapped' (LOGICAL 4 0 0 LOGICAL ()) () (UNKNOWN-FL UNKNOWN-INTENT
UNKNOWN-PROC UNKNOWN UNKNOWN 0 0) UNKNOWN-ACCESS (CONSTANT (LOGICAL 4 0
0 LOGICAL ()) 0 1)) (261 'field_type' (INTEGER 4 0 0 INTEGER ()) () (
UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0)
UNKNOWN-ACCESS (CONSTANT (INTEGER 4 0 0 INTEGER ()) 0 '0')) (262 'bc' (
DERIVED 263 0 0 DERIVED ()) () (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0 POINTER) UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN
()) 0)) (264 'name' (CHARACTER 1 0 0 CHARACTER ((CONSTANT (INTEGER 4 0 0
INTEGER ()) 0 '101'))) () (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0) UNKNOWN-ACCESS ()) (265 'option_path' (CHARACTER 1
0 0 CHARACTER ((CONSTANT (INTEGER 4 0 0 INTEGER ()) 0 '8192'))) () (
UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0)
UNKNOWN-ACCESS (CONSTANT (CHARACTER 1 0 0 CHARACTER (())) 0 8192
'/uninitialised_path/                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            '))
(266 'mesh' (DERIVED 118 0 0 DERIVED ()) () (UNKNOWN-FL UNKNOWN-INTENT
UNKNOWN-PROC UNKNOWN UNKNOWN 0 0) UNKNOWN-ACCESS (STRUCTURE (DERIVED 118
0 0 DERIVED ()) 0 ((() ()) ((CONSTANT (LOGICAL 4 0 0 LOGICAL ()) 0 1) ())
((STRUCTURE (DERIVED 119 0 0 DERIVED ()) 0 ((() ()) (() ()) (() ()) (()
()) ((NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0) ()) ((NULL (UNKNOWN 0 0 0
UNKNOWN ()) 0) ()) ((NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0) ()) ((NULL (
UNKNOWN 0 0 0 UNKNOWN ()) 0) ()) ((NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0) ())
((NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0) ()) ((NULL (UNKNOWN 0 0 0 UNKNOWN ())
0) ()) ((STRUCTURE (DERIVED 120 0 0 DERIVED ()) 0 ((() ()) (() ()) (() ())
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
UNKNOWN ()) 0) ()) ((NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0) ()) ((CONSTANT (
LOGICAL 4 0 0 LOGICAL ()) 0 0) ())) ())) (267 'refcount' (DERIVED 76 0 0
DERIVED ()) () (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0
0 POINTER) UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0)) (268
'aliased' (LOGICAL 4 0 0 LOGICAL ()) () (UNKNOWN-FL UNKNOWN-INTENT
UNKNOWN-PROC UNKNOWN UNKNOWN 0 0) UNKNOWN-ACCESS (CONSTANT (LOGICAL 4 0
0 LOGICAL ()) 0 0)) (269 'py_locweight' (REAL 8 0 0 REAL ()) (2 0
DEFERRED () () () ()) (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN
UNKNOWN 0 0 DIMENSION POINTER) UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0
UNKNOWN ()) 0)) (270 'py_func' (CHARACTER 1 0 0 CHARACTER ((CONSTANT (
INTEGER 4 0 0 INTEGER ()) 0 '8192'))) () (UNKNOWN-FL UNKNOWN-INTENT
UNKNOWN-PROC UNKNOWN UNKNOWN 0 0) UNKNOWN-ACCESS ()) (271 'py_positions'
(DERIVED 108 0 0 DERIVED ()) () (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0 POINTER) UNKNOWN-ACCESS ()) (272
'py_positions_same_mesh' (LOGICAL 4 0 0 LOGICAL ()) () (UNKNOWN-FL
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0) UNKNOWN-ACCESS ()) (
273 'py_dim' (INTEGER 4 0 0 INTEGER ()) () (UNKNOWN-FL UNKNOWN-INTENT
UNKNOWN-PROC UNKNOWN UNKNOWN 0 0) UNKNOWN-ACCESS ()) (274
'py_positions_shape' (DERIVED 119 0 0 DERIVED ()) () (UNKNOWN-FL
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 POINTER) UNKNOWN-ACCESS
(NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0))) PUBLIC (() () () ()) () 0 0
64912956)
263 'scalar_boundary_conditions_ptr' 'fields_data_types'
'scalar_boundary_conditions_ptr' 1 ((DERIVED UNKNOWN-INTENT UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0 POINTER_COMP) (UNKNOWN 0 0 0 UNKNOWN ()) 0 0 () () 0
((275 'boundary_condition' (DERIVED 276 0 0 DERIVED ()) (1 0 DEFERRED ()
()) (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DIMENSION POINTER) UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0)))
PUBLIC (() () () ()) () 0 0 47502942)
276 'scalar_boundary_condition' 'fields_data_types'
'scalar_boundary_condition' 1 ((DERIVED UNKNOWN-INTENT UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0 POINTER_COMP) (UNKNOWN 0 0 0 UNKNOWN ()) 0 0 () () 0
((277 'name' (CHARACTER 1 0 0 CHARACTER ((CONSTANT (INTEGER 4 0 0
INTEGER ()) 0 '101'))) () (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0) UNKNOWN-ACCESS ()) (278 'type' (CHARACTER 1 0 0
CHARACTER ((CONSTANT (INTEGER 4 0 0 INTEGER ()) 0 '101'))) () (
UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0)
UNKNOWN-ACCESS (CONSTANT (CHARACTER 1 0 0 CHARACTER (())) 0 101
'                                                                                                     '))
(279 'surface_element_list' (INTEGER 4 0 0 INTEGER ()) (1 0 DEFERRED ()
()) (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DIMENSION POINTER) UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0)) (
280 'surface_node_list' (INTEGER 4 0 0 INTEGER ()) (1 0 DEFERRED () ())
(UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 DIMENSION
POINTER) UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0)) (281
'surface_mesh' (DERIVED 118 0 0 DERIVED ()) () (UNKNOWN-FL
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 POINTER) UNKNOWN-ACCESS
()) (282 'surface_fields' (DERIVED 256 0 0 DERIVED ()) (1 0 DEFERRED ()
()) (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DIMENSION POINTER) UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0)) (
283 'option_path' (CHARACTER 1 0 0 CHARACTER ((CONSTANT (INTEGER 4 0 0
INTEGER ()) 0 '8192'))) () (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0) UNKNOWN-ACCESS (CONSTANT (CHARACTER 1 0 0 CHARACTER
(())) 0 8192
'/uninitialised_path/                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            ')))
PUBLIC (() () () ()) () 0 0 50740068)
47 'derive_l1_from_l2_halo_mesh' 'halos_derivation'
'derive_l1_from_l2_halo_mesh' 1 ((PROCEDURE UNKNOWN-INTENT MODULE-PROC
DECL UNKNOWN 0 0 SUBROUTINE ALWAYS_EXPLICIT) (UNKNOWN 0 0 0 UNKNOWN ())
284 0 (285 286 287) () 0 () () () 0 0)
288 'integer_hash_table' 'integer_hash_table_module' 'integer_hash_table'
1 ((DERIVED UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0) (UNKNOWN 0
0 0 UNKNOWN ()) 0 0 () () 0 ((289 'address' (DERIVED 290 0 0 DERIVED ())
() (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0)
UNKNOWN-ACCESS ())) PUBLIC (() () () ()) () 0 0 67448272)
290 'c_ptr' '__iso_c_binding' '' 1 ((DERIVED UNKNOWN-INTENT UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0 IS_BIND_C IS_C_INTEROP IS_ISO_C) (DERIVED 290 1 1
UNKNOWN ()) 0 0 () () 0 ((291 '__c_ptr_c_address' (INTEGER 8 1 0 INTEGER
()) () (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0)
UNKNOWN-ACCESS ())) UNKNOWN-ACCESS () () 2 39 0)
9 'halo_max_scalar' 'halos_communications' 'halo_max_scalar' 1 ((
PROCEDURE UNKNOWN-INTENT MODULE-PROC DECL UNKNOWN 0 0 SUBROUTINE
ALWAYS_EXPLICIT) (UNKNOWN 0 0 0 UNKNOWN ()) 292 0 (293 294) () 0 () () ()
0 0)
20 'halo_update_vector' 'halos_communications' 'halo_update_vector' 1 (
(PROCEDURE UNKNOWN-INTENT MODULE-PROC DECL UNKNOWN 0 0 SUBROUTINE
ALWAYS_EXPLICIT) (UNKNOWN 0 0 0 UNKNOWN ()) 295 0 (296 297) () 0 () () ()
0 0)
19 'halo_update_tensor' 'halos_communications' 'halo_update_tensor' 1 (
(PROCEDURE UNKNOWN-INTENT MODULE-PROC DECL UNKNOWN 0 0 SUBROUTINE
ALWAYS_EXPLICIT) (UNKNOWN 0 0 0 UNKNOWN ()) 298 0 (299 300) () 0 () () ()
0 0)
11 'halo_max_array_real' 'halos_communications' 'halo_max_array_real' 1
((PROCEDURE UNKNOWN-INTENT MODULE-PROC DECL UNKNOWN 0 0 SUBROUTINE
ALWAYS_EXPLICIT) (UNKNOWN 0 0 0 UNKNOWN ()) 301 0 (302 303) () 0 () () ()
0 0)
10 'halo_max_scalar_on_halo' 'halos_communications'
'halo_max_scalar_on_halo' 1 ((PROCEDURE UNKNOWN-INTENT MODULE-PROC DECL
UNKNOWN 0 0 SUBROUTINE) (UNKNOWN 0 0 0 UNKNOWN ()) 304 0 (305 306) () 0
() () () 0 0)
22 'halo_update_tensor_on_halo' 'halos_communications'
'halo_update_tensor_on_halo' 1 ((PROCEDURE UNKNOWN-INTENT MODULE-PROC
DECL UNKNOWN 0 0 SUBROUTINE) (UNKNOWN 0 0 0 UNKNOWN ()) 307 0 (308 309)
() 0 () () () 0 0)
23 'halo_update_vector_on_halo' 'halos_communications'
'halo_update_vector_on_halo' 1 ((PROCEDURE UNKNOWN-INTENT MODULE-PROC
DECL UNKNOWN 0 0 SUBROUTINE) (UNKNOWN 0 0 0 UNKNOWN ()) 310 0 (311 312)
() 0 () () () 0 0)
26 'halo_update_array_real_block' 'halos_communications'
'halo_update_array_real_block' 1 ((PROCEDURE UNKNOWN-INTENT MODULE-PROC
DECL UNKNOWN 0 0 SUBROUTINE ALWAYS_EXPLICIT) (UNKNOWN 0 0 0 UNKNOWN ())
313 0 (314 315) () 0 () () () 0 0)
27 'halo_update_array_real' 'halos_communications'
'halo_update_array_real' 1 ((PROCEDURE UNKNOWN-INTENT MODULE-PROC DECL
UNKNOWN 0 0 SUBROUTINE ALWAYS_EXPLICIT) (UNKNOWN 0 0 0 UNKNOWN ()) 316 0
(317 318) () 0 () () () 0 0)
25 'halo_update_array_real_block2' 'halos_communications'
'halo_update_array_real_block2' 1 ((PROCEDURE UNKNOWN-INTENT MODULE-PROC
DECL UNKNOWN 0 0 SUBROUTINE ALWAYS_EXPLICIT) (UNKNOWN 0 0 0 UNKNOWN ())
319 0 (320 321) () 0 () () () 0 0)
24 'halo_update_scalar_on_halo' 'halos_communications'
'halo_update_scalar_on_halo' 1 ((PROCEDURE UNKNOWN-INTENT MODULE-PROC
DECL UNKNOWN 0 0 SUBROUTINE) (UNKNOWN 0 0 0 UNKNOWN ()) 322 0 (323 324)
() 0 () () () 0 0)
28 'halo_update_array_integer_star' 'halos_communications'
'halo_update_array_integer_star' 1 ((PROCEDURE UNKNOWN-INTENT
MODULE-PROC DECL UNKNOWN 0 0 SUBROUTINE) (UNKNOWN 0 0 0 UNKNOWN ()) 325
0 (326 327 328) () 0 () () () 0 0)
30 'halo_update_array_integer_block' 'halos_communications'
'halo_update_array_integer_block' 1 ((PROCEDURE UNKNOWN-INTENT
MODULE-PROC DECL UNKNOWN 0 0 SUBROUTINE ALWAYS_EXPLICIT) (UNKNOWN 0 0 0
UNKNOWN ()) 329 0 (330 331) () 0 () () () 0 0)
15 'halo_verifies_vector_dim' 'halos_communications'
'halo_verifies_vector_dim' 1 ((PROCEDURE UNKNOWN-INTENT MODULE-PROC DECL
UNKNOWN 0 0 FUNCTION) (LOGICAL 4 0 0 LOGICAL ()) 332 0 (333 334 335) ()
336 () () () 0 0)
17 'halo_verifies_array_real' 'halos_communications'
'halo_verifies_array_real' 1 ((PROCEDURE UNKNOWN-INTENT MODULE-PROC DECL
UNKNOWN 0 0 FUNCTION ALWAYS_EXPLICIT) (LOGICAL 4 0 0 LOGICAL ()) 337 0 (
338 339) () 340 () () () 0 0)
16 'halo_verifies_scalar' 'halos_communications' 'halo_verifies_scalar'
1 ((PROCEDURE UNKNOWN-INTENT MODULE-PROC DECL UNKNOWN 0 0 FUNCTION) (
LOGICAL 4 0 0 LOGICAL ()) 341 0 (342 343) () 344 () () () 0 0)
14 'halo_verifies_vector' 'halos_communications' 'halo_verifies_vector'
1 ((PROCEDURE UNKNOWN-INTENT MODULE-PROC DECL UNKNOWN 0 0 FUNCTION) (
LOGICAL 4 0 0 LOGICAL ()) 345 0 (346 347) () 348 () () () 0 0)
18 'halo_verifies_array_integer' 'halos_communications'
'halo_verifies_array_integer' 1 ((PROCEDURE UNKNOWN-INTENT MODULE-PROC
DECL UNKNOWN 0 0 FUNCTION ALWAYS_EXPLICIT) (LOGICAL 4 0 0 LOGICAL ())
349 0 (350 351) () 352 () () () 0 0)
31 'halo_update_array_integer' 'halos_communications'
'halo_update_array_integer' 1 ((PROCEDURE UNKNOWN-INTENT MODULE-PROC
DECL UNKNOWN 0 0 SUBROUTINE ALWAYS_EXPLICIT) (UNKNOWN 0 0 0 UNKNOWN ())
353 0 (354 355) () 0 () () () 0 0)
29 'halo_update_array_integer_block2' 'halos_communications'
'halo_update_array_integer_block2' 1 ((PROCEDURE UNKNOWN-INTENT
MODULE-PROC DECL UNKNOWN 0 0 SUBROUTINE ALWAYS_EXPLICIT) (UNKNOWN 0 0 0
UNKNOWN ()) 356 0 (357 358) () 0 () () () 0 0)
359 'tensor_field' 'fields_data_types' 'tensor_field' 1 ((DERIVED
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 POINTER_COMP) (UNKNOWN 0
0 0 UNKNOWN ()) 0 0 () () 0 ((360 'val' (REAL 8 0 0 REAL ()) (3 0
DEFERRED () () () () () ()) (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0 DIMENSION POINTER) UNKNOWN-ACCESS ()) (361 'wrapped'
(LOGICAL 4 0 0 LOGICAL ()) () (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0) UNKNOWN-ACCESS (CONSTANT (LOGICAL 4 0 0 LOGICAL ())
0 1)) (362 'field_type' (INTEGER 4 0 0 INTEGER ()) () (UNKNOWN-FL
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0) UNKNOWN-ACCESS (
CONSTANT (INTEGER 4 0 0 INTEGER ()) 0 '0')) (363 'name' (CHARACTER 1 0 0
CHARACTER ((CONSTANT (INTEGER 4 0 0 INTEGER ()) 0 '101'))) () (
UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0)
UNKNOWN-ACCESS ()) (364 'dim' (INTEGER 4 0 0 INTEGER ()) (1 0 EXPLICIT (
CONSTANT (INTEGER 4 0 0 INTEGER ()) 0 '1') (CONSTANT (INTEGER 4 0 0
INTEGER ()) 0 '2')) (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN
UNKNOWN 0 0 DIMENSION) UNKNOWN-ACCESS ()) (365 'option_path' (CHARACTER
1 0 0 CHARACTER ((CONSTANT (INTEGER 4 0 0 INTEGER ()) 0 '8192'))) () (
UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0)
UNKNOWN-ACCESS (CONSTANT (CHARACTER 1 0 0 CHARACTER (())) 0 8192
'/uninitialised_path/                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            '))
(366 'mesh' (DERIVED 118 0 0 DERIVED ()) () (UNKNOWN-FL UNKNOWN-INTENT
UNKNOWN-PROC UNKNOWN UNKNOWN 0 0) UNKNOWN-ACCESS (STRUCTURE (DERIVED 118
0 0 DERIVED ()) 0 ((() ()) ((CONSTANT (LOGICAL 4 0 0 LOGICAL ()) 0 1) ())
((STRUCTURE (DERIVED 119 0 0 DERIVED ()) 0 ((() ()) (() ()) (() ()) (()
()) ((NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0) ()) ((NULL (UNKNOWN 0 0 0
UNKNOWN ()) 0) ()) ((NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0) ()) ((NULL (
UNKNOWN 0 0 0 UNKNOWN ()) 0) ()) ((NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0) ())
((NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0) ()) ((NULL (UNKNOWN 0 0 0 UNKNOWN ())
0) ()) ((STRUCTURE (DERIVED 120 0 0 DERIVED ()) 0 ((() ()) (() ()) (() ())
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
UNKNOWN ()) 0) ()) ((NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0) ()) ((CONSTANT (
LOGICAL 4 0 0 LOGICAL ()) 0 0) ())) ())) (367 'refcount' (DERIVED 76 0 0
DERIVED ()) () (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0
0 POINTER) UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0)) (368
'aliased' (LOGICAL 4 0 0 LOGICAL ()) () (UNKNOWN-FL UNKNOWN-INTENT
UNKNOWN-PROC UNKNOWN UNKNOWN 0 0) UNKNOWN-ACCESS (CONSTANT (LOGICAL 4 0
0 LOGICAL ()) 0 0))) PUBLIC (() () () ()) () 0 0 25110185)
21 'halo_update_scalar' 'halos_communications' 'halo_update_scalar' 1 (
(PROCEDURE UNKNOWN-INTENT MODULE-PROC DECL UNKNOWN 0 0 SUBROUTINE
ALWAYS_EXPLICIT) (UNKNOWN 0 0 0 UNKNOWN ()) 369 0 (370 371) () 0 () () ()
0 0)
12 'halo_universal_number_vector' 'halos_numbering'
'halo_universal_number_vector' 1 ((PROCEDURE UNKNOWN-INTENT MODULE-PROC
DECL UNKNOWN 0 0 DIMENSION FUNCTION ALWAYS_EXPLICIT) (INTEGER 4 0 0
INTEGER ()) 372 0 (373 374) (1 0 EXPLICIT (CONSTANT (INTEGER 4 0 0
INTEGER ()) 0 '1') (FUNCTION (INTEGER 4 0 0 INTEGER ()) 0 375 (('' (
VARIABLE (INTEGER 4 0 0 INTEGER ()) 1 374 ((ARRAY (FULL 0))))) ('' ()) (
'' ())) '' 0 'size')) 376 () () () 0 0)
6 'get_universal_numbering_multiple_components' 'halos_numbering'
'get_universal_numbering_multiple_components' 1 ((PROCEDURE
UNKNOWN-INTENT MODULE-PROC DECL UNKNOWN 0 0 SUBROUTINE ALWAYS_EXPLICIT)
(UNKNOWN 0 0 0 UNKNOWN ()) 377 0 (378 379) () 0 () () () 0 0)
34 'node_count_halo' 'halos_base' 'node_count_halo' 1 ((PROCEDURE
UNKNOWN-INTENT MODULE-PROC DECL UNKNOWN 0 0 FUNCTION PURE) (INTEGER 4 0
0 INTEGER ()) 380 0 (381) () 382 () () () 0 0)
33 'node_owned_halo' 'halos_ownership' 'node_owned_halo' 1 ((PROCEDURE
UNKNOWN-INTENT MODULE-PROC DECL UNKNOWN 0 0 FUNCTION) (LOGICAL 4 0 0
LOGICAL ()) 383 0 (384 385) () 33 () () () 0 0)
41 'nodes_owned_halo' 'halos_ownership' 'nodes_owned_halo' 1 ((
PROCEDURE UNKNOWN-INTENT MODULE-PROC DECL UNKNOWN 0 0 DIMENSION FUNCTION
ALWAYS_EXPLICIT) (LOGICAL 4 0 0 LOGICAL ()) 386 0 (387 388) (1 0
EXPLICIT (CONSTANT (INTEGER 4 0 0 INTEGER ()) 0 '1') (FUNCTION (INTEGER
4 0 0 INTEGER ()) 0 389 (('' (VARIABLE (INTEGER 4 0 0 INTEGER ()) 1 388
((ARRAY (FULL 0))))) ('' ()) ('' ())) '' 0 'size')) 390 () () () 0 0)
391 'halo_proc_count' 'halos_base' 'halo_proc_count' 1 ((PROCEDURE
UNKNOWN-INTENT MODULE-PROC DECL UNKNOWN 0 0 FUNCTION PURE) (INTEGER 4 0
0 INTEGER ()) 392 0 (393) () 391 () () () 0 0)
46 'serial_storage_halo_single' 'halos_base' 'serial_storage_halo_single'
1 ((PROCEDURE UNKNOWN-INTENT MODULE-PROC DECL UNKNOWN 0 0 FUNCTION) (
LOGICAL 4 0 0 LOGICAL ()) 394 0 (395) () 396 () () () 0 0)
39 'read_halos_mesh' 'halos_registration' 'read_halos_mesh' 1 ((
PROCEDURE UNKNOWN-INTENT MODULE-PROC DECL UNKNOWN 0 0 SUBROUTINE
ALWAYS_EXPLICIT) (UNKNOWN 0 0 0 UNKNOWN ()) 397 0 (398 399 400) () 0 ()
() () 0 0)
38 'read_halos_positions' 'halos_registration' 'read_halos_positions' 1
((PROCEDURE UNKNOWN-INTENT MODULE-PROC DECL UNKNOWN 0 0 SUBROUTINE
ALWAYS_EXPLICIT) (UNKNOWN 0 0 0 UNKNOWN ()) 401 0 (402 403 404) () 0 ()
() () 0 0)
45 'serial_storage_halo_multiple' 'halos_base'
'serial_storage_halo_multiple' 1 ((PROCEDURE UNKNOWN-INTENT MODULE-PROC
DECL UNKNOWN 0 0 DIMENSION FUNCTION ALWAYS_EXPLICIT) (LOGICAL 4 0 0
LOGICAL ()) 405 0 (406) (1 0 EXPLICIT (CONSTANT (INTEGER 4 0 0 INTEGER ())
0 '1') (FUNCTION (INTEGER 4 0 0 INTEGER ()) 0 407 (('' (VARIABLE (
DERIVED 144 0 0 DERIVED ()) 1 406 ((ARRAY (FULL 0))))) ('' ()) ('' ())) ''
0 'size')) 408 () () () 0 0)
8 'halo_communicator_halo' 'halos_base' 'halo_communicator_halo' 1 ((
PROCEDURE UNKNOWN-INTENT MODULE-PROC DECL UNKNOWN 0 0 FUNCTION PURE) (
INTEGER 4 0 0 INTEGER ()) 409 0 (410) () 411 () () () 0 0)
43 'reorder_halo_vector' 'halos_repair' 'reorder_halo_vector' 1 ((
PROCEDURE UNKNOWN-INTENT MODULE-PROC DECL UNKNOWN 0 0 SUBROUTINE) (
UNKNOWN 0 0 0 UNKNOWN ()) 412 0 (413 414) () 0 () () () 0 0)
42 'reorder_halo_halo' 'halos_repair' 'reorder_halo_halo' 1 ((PROCEDURE
UNKNOWN-INTENT MODULE-PROC DECL UNKNOWN 0 0 SUBROUTINE) (UNKNOWN 0 0 0
UNKNOWN ()) 415 0 (416 417) () 0 () () () 0 0)
48 'zero_halo' 'halos_base' 'zero_halo' 1 ((PROCEDURE UNKNOWN-INTENT
MODULE-PROC DECL UNKNOWN 0 0 SUBROUTINE) (UNKNOWN 0 0 0 UNKNOWN ()) 418
0 (419) () 0 () () () 0 0)
4 'deallocate_halo_vector' 'halos_allocates' 'deallocate_halo_vector' 1
((PROCEDURE UNKNOWN-INTENT MODULE-PROC DECL UNKNOWN 0 0 SUBROUTINE
ALWAYS_EXPLICIT) (UNKNOWN 0 0 0 UNKNOWN ()) 420 0 (421) () 0 () () () 0
0)
422 'create_global_to_universal_numbering' 'halos_numbering'
'create_global_to_universal_numbering' 1 ((PROCEDURE UNKNOWN-INTENT
MODULE-PROC DECL UNKNOWN 0 0 SUBROUTINE ALWAYS_EXPLICIT) (UNKNOWN 0 0 0
UNKNOWN ()) 423 0 (424 425) () 0 () () () 0 0)
426 'create_ownership' 'halos_ownership' 'create_ownership' 1 ((
PROCEDURE UNKNOWN-INTENT MODULE-PROC DECL UNKNOWN 0 0 SUBROUTINE) (
UNKNOWN 0 0 0 UNKNOWN ()) 427 0 (428) () 0 () () () 0 0)
429 'deallocate_ownership_cache' 'halos_allocates'
'deallocate_ownership_cache' 1 ((PROCEDURE UNKNOWN-INTENT MODULE-PROC
DECL UNKNOWN 0 0 SUBROUTINE IMPLICIT_PURE) (UNKNOWN 0 0 0 UNKNOWN ())
430 0 (431) () 0 () () () 0 0)
432 'deallocate_universal_numbering_cache' 'halos_allocates'
'deallocate_universal_numbering_cache' 1 ((PROCEDURE UNKNOWN-INTENT
MODULE-PROC DECL UNKNOWN 0 0 SUBROUTINE IMPLICIT_PURE) (UNKNOWN 0 0 0
UNKNOWN ()) 433 0 (434) () 0 () () () 0 0)
435 'derive_element_halo_from_node_halo' 'halos_derivation'
'derive_element_halo_from_node_halo' 1 ((PROCEDURE UNKNOWN-INTENT
MODULE-PROC DECL UNKNOWN 0 0 SUBROUTINE ALWAYS_EXPLICIT) (UNKNOWN 0 0 0
UNKNOWN ()) 436 0 (437 438 439) () 0 () () () 0 0)
440 'derive_maximal_surface_element_halo' 'halos_derivation'
'derive_maximal_surface_element_halo' 1 ((PROCEDURE UNKNOWN-INTENT
MODULE-PROC DECL UNKNOWN 0 0 FUNCTION ALWAYS_EXPLICIT) (DERIVED 144 0 0
DERIVED ()) 441 0 (442 443 444 445) () 446 () () () 0 0)
447 'derive_nonperiodic_halos_from_periodic_halos' 'halos_derivation'
'derive_nonperiodic_halos_from_periodic_halos' 1 ((PROCEDURE
UNKNOWN-INTENT MODULE-PROC DECL UNKNOWN 0 0 SUBROUTINE) (UNKNOWN 0 0 0
UNKNOWN ()) 448 0 (449 450 451) () 0 () () () 0 0)
452 'derive_sub_halo' 'halos_derivation' 'derive_sub_halo' 1 ((
PROCEDURE UNKNOWN-INTENT MODULE-PROC DECL UNKNOWN 0 0 FUNCTION
ALWAYS_EXPLICIT) (DERIVED 144 0 0 DERIVED ()) 453 0 (454 455) () 456 ()
() () 0 0)
457 'ele_owner' 'halos_derivation' 'ele_owner' 1 ((PROCEDURE
UNKNOWN-INTENT MODULE-PROC DECL UNKNOWN 0 0 FUNCTION) (INTEGER 4 0 0
INTEGER ()) 458 0 (459 460 461) () 457 () () () 0 0)
462 'ewrite_universal_numbers' 'halos_numbering'
'ewrite_universal_numbers' 1 ((PROCEDURE UNKNOWN-INTENT MODULE-PROC DECL
UNKNOWN 0 0 SUBROUTINE) (UNKNOWN 0 0 0 UNKNOWN ()) 463 0 (464 465) () 0
() () () 0 0)
466 'extract_all_halo_receives' 'halos_base' 'extract_all_halo_receives'
1 ((PROCEDURE UNKNOWN-INTENT MODULE-PROC DECL UNKNOWN 0 0 SUBROUTINE
ALWAYS_EXPLICIT) (UNKNOWN 0 0 0 UNKNOWN ()) 467 0 (468 469 470 471) () 0
() () () 0 0)
472 'extract_all_halo_sends' 'halos_base' 'extract_all_halo_sends' 1 ((
PROCEDURE UNKNOWN-INTENT MODULE-PROC DECL UNKNOWN 0 0 SUBROUTINE
ALWAYS_EXPLICIT) (UNKNOWN 0 0 0 UNKNOWN ()) 473 0 (474 475 476 477) () 0
() () () 0 0)
478 'extract_raw_halo_data' 'halos_registration' 'extract_raw_halo_data'
1 ((PROCEDURE UNKNOWN-INTENT MODULE-PROC DECL UNKNOWN 0 0 SUBROUTINE
ALWAYS_EXPLICIT) (UNKNOWN 0 0 0 UNKNOWN ()) 479 0 (480 481 482 483 484
485) () 0 () () () 0 0)
486 'form_halo_from_raw_data' 'halos_registration'
'form_halo_from_raw_data' 1 ((PROCEDURE UNKNOWN-INTENT MODULE-PROC DECL
UNKNOWN 0 0 SUBROUTINE ALWAYS_EXPLICIT) (UNKNOWN 0 0 0 UNKNOWN ()) 487 0
(488 489 490 491 492 493 494 495 496) () 0 () () () 0 0)
497 'generate_substate_halos' 'halos_registration'
'generate_substate_halos' 1 ((PROCEDURE UNKNOWN-INTENT MODULE-PROC DECL
UNKNOWN 0 0 SUBROUTINE ALWAYS_EXPLICIT) (UNKNOWN 0 0 0 UNKNOWN ()) 498 0
(499 500 501 502) () 0 () () () 0 0)
503 'get_node_owners' 'halos_ownership' 'get_node_owners' 1 ((PROCEDURE
UNKNOWN-INTENT MODULE-PROC DECL UNKNOWN 0 0 SUBROUTINE ALWAYS_EXPLICIT)
(UNKNOWN 0 0 0 UNKNOWN ()) 504 0 (505 506) () 0 () () () 0 0)
507 'get_owned_nodes' 'halos_ownership' 'get_owned_nodes' 1 ((PROCEDURE
UNKNOWN-INTENT MODULE-PROC DECL UNKNOWN 0 0 SUBROUTINE ALWAYS_EXPLICIT)
(UNKNOWN 0 0 0 UNKNOWN ()) 508 0 (509 510) () 0 () () () 0 0)
7 'get_universal_numbering' 'halos_numbering' 'get_universal_numbering'
1 ((PROCEDURE UNKNOWN-INTENT MODULE-PROC DECL UNKNOWN 0 0 SUBROUTINE
GENERIC) (UNKNOWN 0 0 0 UNKNOWN ()) 511 0 (512 513) () 0 () () () 0 0)
514 'get_universal_numbering_inverse' 'halos_numbering'
'get_universal_numbering_inverse' 1 ((PROCEDURE UNKNOWN-INTENT
MODULE-PROC DECL UNKNOWN 0 0 SUBROUTINE) (UNKNOWN 0 0 0 UNKNOWN ()) 515
0 (516 517) () 0 () () () 0 0)
518 'halo_all_receives_count' 'halos_base' 'halo_all_receives_count' 1 (
(PROCEDURE UNKNOWN-INTENT MODULE-PROC DECL UNKNOWN 0 0 FUNCTION PURE) (
INTEGER 4 0 0 INTEGER ()) 519 0 (520) () 518 () () () 0 0)
521 'halo_all_sends_count' 'halos_base' 'halo_all_sends_count' 1 ((
PROCEDURE UNKNOWN-INTENT MODULE-PROC DECL UNKNOWN 0 0 FUNCTION PURE) (
INTEGER 4 0 0 INTEGER ()) 522 0 (523) () 521 () () () 0 0)
524 'halo_all_unique_receives_count' 'halos_base'
'halo_all_unique_receives_count' 1 ((PROCEDURE UNKNOWN-INTENT
MODULE-PROC DECL UNKNOWN 0 0 FUNCTION) (INTEGER 4 0 0 INTEGER ()) 525 0
(526) () 524 () () () 0 0)
527 'halo_data_type' 'halos_base' 'halo_data_type' 1 ((PROCEDURE
UNKNOWN-INTENT MODULE-PROC DECL UNKNOWN 0 0 FUNCTION PURE) (INTEGER 4 0
0 INTEGER ()) 528 0 (529) () 527 () () () 0 0)
530 'halo_name' 'halos_base' 'halo_name' 1 ((PROCEDURE UNKNOWN-INTENT
MODULE-PROC DECL UNKNOWN 0 0 FUNCTION PURE) (CHARACTER 1 0 0 CHARACTER (
(FUNCTION (INTEGER 4 0 0 INTEGER ()) 0 531 (('' (VARIABLE (CHARACTER 1 0
0 CHARACTER ((CONSTANT (INTEGER 4 0 0 INTEGER ()) 0 '101'))) 0 532 ((
COMPONENT 144 533 'name')))) ('' ())) '__len_trim1' 0 'lnblnk'))) 534 0
(532) () 530 () () () 0 0)
535 'halo_node_owner' 'halos_ownership' 'halo_node_owner' 1 ((PROCEDURE
UNKNOWN-INTENT MODULE-PROC DECL UNKNOWN 0 0 FUNCTION ALWAYS_EXPLICIT) (
INTEGER 4 0 0 INTEGER ()) 536 0 (537 538 539) () 540 () () () 0 0)
541 'halo_node_owners' 'halos_ownership' 'halo_node_owners' 1 ((
PROCEDURE UNKNOWN-INTENT MODULE-PROC DECL UNKNOWN 0 0 DIMENSION FUNCTION
ALWAYS_EXPLICIT) (INTEGER 4 0 0 INTEGER ()) 542 0 (543 544 545) (1 0
EXPLICIT (CONSTANT (INTEGER 4 0 0 INTEGER ()) 0 '1') (FUNCTION (INTEGER
4 0 0 INTEGER ()) 0 546 (('' (VARIABLE (INTEGER 4 0 0 INTEGER ()) 1 544
((ARRAY (FULL 0))))) ('' ()) ('' ())) '' 0 'size')) 547 () () () 0 0)
548 'halo_nowned_nodes' 'halos_base' 'halo_nowned_nodes' 1 ((PROCEDURE
UNKNOWN-INTENT MODULE-PROC DECL UNKNOWN 0 0 FUNCTION PURE) (INTEGER 4 0
0 INTEGER ()) 549 0 (550) () 548 () () () 0 0)
551 'halo_order_general' 'halo_data_types' 'halo_order_general' 1 ((
PARAMETER UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN IMPLICIT-SAVE 0 0) (
INTEGER 4 0 0 INTEGER ()) 0 0 () (CONSTANT (INTEGER 4 0 0 INTEGER ()) 0
'1') () 0 () () () 0 0)
552 'halo_order_trailing_receives' 'halo_data_types'
'halo_order_trailing_receives' 1 ((PARAMETER UNKNOWN-INTENT UNKNOWN-PROC
UNKNOWN IMPLICIT-SAVE 0 0) (INTEGER 4 0 0 INTEGER ()) 0 0 () (CONSTANT (
INTEGER 4 0 0 INTEGER ()) 0 '2') () 0 () () () 0 0)
553 'halo_ordering_scheme' 'halos_base' 'halo_ordering_scheme' 1 ((
PROCEDURE UNKNOWN-INTENT MODULE-PROC DECL UNKNOWN 0 0 FUNCTION PURE) (
INTEGER 4 0 0 INTEGER ()) 554 0 (555) () 553 () () () 0 0)
556 'halo_pointer' 'halo_data_types' 'halo_pointer' 1 ((DERIVED
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 POINTER_COMP) (UNKNOWN 0
0 0 UNKNOWN ()) 0 0 () () 0 ((557 'ptr' (DERIVED 144 0 0 DERIVED ()) ()
(UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 POINTER)
UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0))) PUBLIC (() () () ())
() 0 0 72023282)
558 'halo_receive' 'halos_base' 'halo_receive' 1 ((PROCEDURE
UNKNOWN-INTENT MODULE-PROC DECL UNKNOWN 0 0 FUNCTION) (INTEGER 4 0 0
INTEGER ()) 559 0 (560 561 562) () 558 () () () 0 0)
563 'halo_receive_count' 'halos_base' 'halo_receive_count' 1 ((
PROCEDURE UNKNOWN-INTENT MODULE-PROC DECL UNKNOWN 0 0 FUNCTION PURE) (
INTEGER 4 0 0 INTEGER ()) 564 0 (565 566) () 563 () () () 0 0)
567 'halo_receive_counts' 'halos_base' 'halo_receive_counts' 1 ((
PROCEDURE UNKNOWN-INTENT MODULE-PROC DECL UNKNOWN 0 0 SUBROUTINE
IMPLICIT_PURE) (UNKNOWN 0 0 0 UNKNOWN ()) 568 0 (569 570) () 0 () () ()
0 0)
571 'halo_receives' 'halos_base' 'halo_receives' 1 ((PROCEDURE
UNKNOWN-INTENT MODULE-PROC DECL UNKNOWN 0 0 DIMENSION POINTER FUNCTION
ALWAYS_EXPLICIT) (INTEGER 4 0 0 INTEGER ()) 572 0 (573 574) (1 0
DEFERRED () ()) 571 () () () 0 0)
575 'halo_send' 'halos_base' 'halo_send' 1 ((PROCEDURE UNKNOWN-INTENT
MODULE-PROC DECL UNKNOWN 0 0 FUNCTION) (INTEGER 4 0 0 INTEGER ()) 576 0
(577 578 579) () 575 () () () 0 0)
580 'halo_send_count' 'halos_base' 'halo_send_count' 1 ((PROCEDURE
UNKNOWN-INTENT MODULE-PROC DECL UNKNOWN 0 0 FUNCTION PURE) (INTEGER 4 0
0 INTEGER ()) 581 0 (582 583) () 580 () () () 0 0)
584 'halo_send_counts' 'halos_base' 'halo_send_counts' 1 ((PROCEDURE
UNKNOWN-INTENT MODULE-PROC DECL UNKNOWN 0 0 SUBROUTINE IMPLICIT_PURE) (
UNKNOWN 0 0 0 UNKNOWN ()) 585 0 (586 587) () 0 () () () 0 0)
588 'halo_sends' 'halos_base' 'halo_sends' 1 ((PROCEDURE UNKNOWN-INTENT
MODULE-PROC DECL UNKNOWN 0 0 DIMENSION POINTER FUNCTION ALWAYS_EXPLICIT)
(INTEGER 4 0 0 INTEGER ()) 589 0 (590 591) (1 0 DEFERRED () ()) 588 () ()
() 0 0)
144 'halo_type' 'halo_data_types' 'halo_type' 1 ((DERIVED UNKNOWN-INTENT
UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 POINTER_COMP) (UNKNOWN 0 0 0 UNKNOWN ())
0 0 () () 0 ((533 'name' (CHARACTER 1 0 0 CHARACTER ((CONSTANT (INTEGER
4 0 0 INTEGER ()) 0 '101'))) () (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0) UNKNOWN-ACCESS ()) (592 'refcount' (DERIVED 76 0 0
DERIVED ()) () (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0
0 POINTER) UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0)) (593
'data_type' (INTEGER 4 0 0 INTEGER ()) () (UNKNOWN-FL UNKNOWN-INTENT
UNKNOWN-PROC UNKNOWN UNKNOWN 0 0) UNKNOWN-ACCESS (CONSTANT (INTEGER 4 0
0 INTEGER ()) 0 '0')) (594 'ordering_scheme' (INTEGER 4 0 0 INTEGER ())
() (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0)
UNKNOWN-ACCESS (CONSTANT (INTEGER 4 0 0 INTEGER ()) 0 '0')) (595
'communicator' (INTEGER 4 0 0 INTEGER ()) () (UNKNOWN-FL UNKNOWN-INTENT
UNKNOWN-PROC UNKNOWN UNKNOWN 0 0) UNKNOWN-ACCESS ()) (596 'nprocs' (
INTEGER 4 0 0 INTEGER ()) () (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0) UNKNOWN-ACCESS (CONSTANT (INTEGER 4 0 0 INTEGER ())
0 '0')) (597 'sends' (DERIVED 84 0 0 DERIVED ()) (1 0 DEFERRED () ()) (
UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 DIMENSION
POINTER) UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0)) (598
'receives' (DERIVED 84 0 0 DERIVED ()) (1 0 DEFERRED () ()) (UNKNOWN-FL
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 DIMENSION POINTER)
UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0)) (599 'nowned_nodes'
(INTEGER 4 0 0 INTEGER ()) () (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0) UNKNOWN-ACCESS (CONSTANT (INTEGER 4 0 0 INTEGER ())
0 '-1')) (600 'owners' (INTEGER 4 0 0 INTEGER ()) (1 0 DEFERRED () ()) (
UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 DIMENSION
POINTER) UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0)) (601
'unn_count' (INTEGER 4 0 0 INTEGER ()) () (UNKNOWN-FL UNKNOWN-INTENT
UNKNOWN-PROC UNKNOWN UNKNOWN 0 0) UNKNOWN-ACCESS (CONSTANT (INTEGER 4 0
0 INTEGER ()) 0 '-1')) (602 'owned_nodes_unn_base' (INTEGER 4 0 0
INTEGER ()) (1 0 DEFERRED () ()) (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0 DIMENSION POINTER) UNKNOWN-ACCESS (NULL (UNKNOWN 0 0
0 UNKNOWN ()) 0)) (603 'my_owned_nodes_unn_base' (INTEGER 4 0 0 INTEGER
()) () (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0)
UNKNOWN-ACCESS (CONSTANT (INTEGER 4 0 0 INTEGER ()) 0 '-1')) (604
'receives_gnn_to_unn' (INTEGER 4 0 0 INTEGER ()) (1 0 DEFERRED () ()) (
UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 DIMENSION
POINTER) UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0)) (605
'gnn_to_unn' (INTEGER 4 0 0 INTEGER ()) (1 0 DEFERRED () ()) (
UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 DIMENSION
POINTER) UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0))) PUBLIC (()
() () ()) () 0 0 63994053)
606 'halo_type_cg_node' 'halo_data_types' 'halo_type_cg_node' 1 ((
PARAMETER UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN IMPLICIT-SAVE 0 0) (
INTEGER 4 0 0 INTEGER ()) 0 0 () (CONSTANT (INTEGER 4 0 0 INTEGER ()) 0
'1') () 0 () () () 0 0)
607 'halo_type_dg_node' 'halo_data_types' 'halo_type_dg_node' 1 ((
PARAMETER UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN IMPLICIT-SAVE 0 0) (
INTEGER 4 0 0 INTEGER ()) 0 0 () (CONSTANT (INTEGER 4 0 0 INTEGER ()) 0
'2') () 0 () () () 0 0)
608 'halo_type_element' 'halo_data_types' 'halo_type_element' 1 ((
PARAMETER UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN IMPLICIT-SAVE 0 0) (
INTEGER 4 0 0 INTEGER ()) 0 0 () (CONSTANT (INTEGER 4 0 0 INTEGER ()) 0
'3') () 0 () () () 0 0)
609 'halo_unique_receive_count' 'halos_base' 'halo_unique_receive_count'
1 ((PROCEDURE UNKNOWN-INTENT MODULE-PROC DECL UNKNOWN 0 0 FUNCTION) (
INTEGER 4 0 0 INTEGER ()) 610 0 (611 612) () 609 () () () 0 0)
613 'halo_universal_node_owners' 'halos_ownership'
'halo_universal_node_owners' 1 ((PROCEDURE UNKNOWN-INTENT MODULE-PROC
DECL UNKNOWN 0 0 DIMENSION FUNCTION ALWAYS_EXPLICIT) (INTEGER 4 0 0
INTEGER ()) 614 0 (615 616) (1 0 EXPLICIT (CONSTANT (INTEGER 4 0 0
INTEGER ()) 0 '1') (FUNCTION (INTEGER 4 0 0 INTEGER ()) 0 617 (('' (
VARIABLE (INTEGER 4 0 0 INTEGER ()) 1 616 ((ARRAY (FULL 0))))) ('' ()) (
'' ())) '' 0 'size')) 618 () () () 0 0)
13 'halo_universal_number' 'halos_numbering' 'halo_universal_number' 1 (
(PROCEDURE UNKNOWN-INTENT MODULE-PROC DECL UNKNOWN 0 0 FUNCTION GENERIC)
(INTEGER 4 0 0 INTEGER ()) 619 0 (620 621) () 622 () () () 0 0)
623 'halo_universal_numbers' 'halos_numbering' 'halo_universal_numbers'
1 ((PROCEDURE UNKNOWN-INTENT MODULE-PROC DECL UNKNOWN 0 0 DIMENSION
FUNCTION ALWAYS_EXPLICIT) (INTEGER 4 0 0 INTEGER ()) 624 0 (625 626) (1
0 EXPLICIT (CONSTANT (INTEGER 4 0 0 INTEGER ()) 0 '1') (FUNCTION (
INTEGER 4 0 0 INTEGER ()) 0 627 (('' (VARIABLE (INTEGER 4 0 0 INTEGER ())
1 626 ((ARRAY (FULL 0))))) ('' ()) ('' ())) '' 0 'size')) 628 () () () 0
0)
629 'halo_valid_for_communication' 'halos_debug'
'halo_valid_for_communication' 1 ((PROCEDURE UNKNOWN-INTENT MODULE-PROC
DECL UNKNOWN 0 0 FUNCTION) (LOGICAL 4 0 0 LOGICAL ()) 630 0 (631) () 632
() () () 0 0)
633 'halos' 'halos' 'halos' 1 ((MODULE UNKNOWN-INTENT UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0) (UNKNOWN 0 0 0 UNKNOWN ()) 0 0 () () 0 () () () 0 0)
634 'has_global_to_universal_numbering' 'halos_numbering'
'has_global_to_universal_numbering' 1 ((PROCEDURE UNKNOWN-INTENT
MODULE-PROC DECL UNKNOWN 0 0 FUNCTION) (LOGICAL 4 0 0 LOGICAL ()) 635 0
(636) () 637 () () () 0 0)
638 'has_nowned_nodes' 'halos_base' 'has_nowned_nodes' 1 ((PROCEDURE
UNKNOWN-INTENT MODULE-PROC DECL UNKNOWN 0 0 FUNCTION PURE) (LOGICAL 4 0
0 LOGICAL ()) 639 0 (640) () 638 () () () 0 0)
641 'has_ownership' 'halos_ownership' 'has_ownership' 1 ((PROCEDURE
UNKNOWN-INTENT MODULE-PROC DECL UNKNOWN 0 0 FUNCTION PURE) (LOGICAL 4 0
0 LOGICAL ()) 642 0 (643) () 641 () () () 0 0)
644 'invert_comms_sizes' 'halos_derivation' 'invert_comms_sizes' 1 ((
PROCEDURE UNKNOWN-INTENT MODULE-PROC DECL UNKNOWN 0 0 DIMENSION FUNCTION
ALWAYS_EXPLICIT) (INTEGER 4 0 0 INTEGER ()) 645 0 (646 647) (1 0
EXPLICIT (CONSTANT (INTEGER 4 0 0 INTEGER ()) 0 '1') (FUNCTION (INTEGER
4 0 0 INTEGER ()) 0 648 (('' (VARIABLE (INTEGER 4 0 0 INTEGER ()) 1 646
((ARRAY (FULL 0))))) ('' ()) ('' ())) '' 0 'size')) 649 () () () 0 0)
650 'max_halo_node' 'halos_base' 'max_halo_node' 1 ((PROCEDURE
UNKNOWN-INTENT MODULE-PROC DECL UNKNOWN 0 0 FUNCTION) (INTEGER 4 0 0
INTEGER ()) 651 0 (652) () 653 () () () 0 0)
654 'max_halo_receive_node' 'halos_base' 'max_halo_receive_node' 1 ((
PROCEDURE UNKNOWN-INTENT MODULE-PROC DECL UNKNOWN 0 0 FUNCTION) (
INTEGER 4 0 0 INTEGER ()) 655 0 (656) () 657 () () () 0 0)
658 'max_halo_send_node' 'halos_base' 'max_halo_send_node' 1 ((
PROCEDURE UNKNOWN-INTENT MODULE-PROC DECL UNKNOWN 0 0 FUNCTION) (
INTEGER 4 0 0 INTEGER ()) 659 0 (660) () 661 () () () 0 0)
662 'min_halo_node' 'halos_base' 'min_halo_node' 1 ((PROCEDURE
UNKNOWN-INTENT MODULE-PROC DECL UNKNOWN 0 0 FUNCTION) (INTEGER 4 0 0
INTEGER ()) 663 0 (664) () 665 () () () 0 0)
666 'min_halo_receive_node' 'halos_base' 'min_halo_receive_node' 1 ((
PROCEDURE UNKNOWN-INTENT MODULE-PROC DECL UNKNOWN 0 0 FUNCTION) (
INTEGER 4 0 0 INTEGER ()) 667 0 (668) () 669 () () () 0 0)
670 'min_halo_send_node' 'halos_base' 'min_halo_send_node' 1 ((
PROCEDURE UNKNOWN-INTENT MODULE-PROC DECL UNKNOWN 0 0 FUNCTION) (
INTEGER 4 0 0 INTEGER ()) 671 0 (672) () 673 () () () 0 0)
674 'print_halo' 'halos_debug' 'print_halo' 1 ((PROCEDURE UNKNOWN-INTENT
MODULE-PROC DECL UNKNOWN 0 0 SUBROUTINE) (UNKNOWN 0 0 0 UNKNOWN ()) 675
0 (676 677) () 0 () () () 0 0)
678 'reorder_element_halo' 'halos_repair' 'reorder_element_halo' 1 ((
PROCEDURE UNKNOWN-INTENT MODULE-PROC DECL UNKNOWN 0 0 SUBROUTINE) (
UNKNOWN 0 0 0 UNKNOWN ()) 679 0 (680 681 682) () 0 () () () 0 0)
683 'reorder_halo_from_element_halo' 'halos_repair'
'reorder_halo_from_element_halo' 1 ((PROCEDURE UNKNOWN-INTENT
MODULE-PROC DECL UNKNOWN 0 0 SUBROUTINE) (UNKNOWN 0 0 0 UNKNOWN ()) 684
0 (685 686 687) () 0 () () () 0 0)
688 'reorder_halo_receives' 'halos_repair' 'reorder_halo_receives' 1 ((
PROCEDURE UNKNOWN-INTENT MODULE-PROC DECL UNKNOWN 0 0 SUBROUTINE) (
UNKNOWN 0 0 0 UNKNOWN ()) 689 0 (690 691) () 0 () () () 0 0)
692 'reorder_l1_from_l2_halo' 'halos_repair' 'reorder_l1_from_l2_halo' 1
((PROCEDURE UNKNOWN-INTENT MODULE-PROC DECL UNKNOWN 0 0 SUBROUTINE
ALWAYS_EXPLICIT) (UNKNOWN 0 0 0 UNKNOWN ()) 693 0 (694 695 696) () 0 ()
() () 0 0)
697 'set_all_halo_receives' 'halos_base' 'set_all_halo_receives' 1 ((
PROCEDURE UNKNOWN-INTENT MODULE-PROC DECL UNKNOWN 0 0 SUBROUTINE
ALWAYS_EXPLICIT) (UNKNOWN 0 0 0 UNKNOWN ()) 698 0 (699 700) () 0 () () ()
0 0)
701 'set_all_halo_sends' 'halos_base' 'set_all_halo_sends' 1 ((
PROCEDURE UNKNOWN-INTENT MODULE-PROC DECL UNKNOWN 0 0 SUBROUTINE
ALWAYS_EXPLICIT) (UNKNOWN 0 0 0 UNKNOWN ()) 702 0 (703 704) () 0 () () ()
0 0)
705 'set_halo_communicator' 'halos_base' 'set_halo_communicator' 1 ((
PROCEDURE UNKNOWN-INTENT MODULE-PROC DECL UNKNOWN 0 0 SUBROUTINE) (
UNKNOWN 0 0 0 UNKNOWN ()) 706 0 (707 708) () 0 () () () 0 0)
709 'set_halo_data_type' 'halos_base' 'set_halo_data_type' 1 ((
PROCEDURE UNKNOWN-INTENT MODULE-PROC DECL UNKNOWN 0 0 SUBROUTINE) (
UNKNOWN 0 0 0 UNKNOWN ()) 710 0 (711 712) () 0 () () () 0 0)
713 'set_halo_name' 'halos_base' 'set_halo_name' 1 ((PROCEDURE
UNKNOWN-INTENT MODULE-PROC DECL UNKNOWN 0 0 SUBROUTINE IMPLICIT_PURE) (
UNKNOWN 0 0 0 UNKNOWN ()) 714 0 (715 716) () 0 () () () 0 0)
717 'set_halo_nowned_nodes' 'halos_base' 'set_halo_nowned_nodes' 1 ((
PROCEDURE UNKNOWN-INTENT MODULE-PROC DECL UNKNOWN 0 0 SUBROUTINE) (
UNKNOWN 0 0 0 UNKNOWN ()) 718 0 (719 720) () 0 () () () 0 0)
721 'set_halo_ordering_scheme' 'halos_base' 'set_halo_ordering_scheme' 1
((PROCEDURE UNKNOWN-INTENT MODULE-PROC DECL UNKNOWN 0 0 SUBROUTINE) (
UNKNOWN 0 0 0 UNKNOWN ()) 722 0 (723 724) () 0 () () () 0 0)
725 'set_halo_receive' 'halos_base' 'set_halo_receive' 1 ((PROCEDURE
UNKNOWN-INTENT MODULE-PROC DECL UNKNOWN 0 0 SUBROUTINE) (UNKNOWN 0 0 0
UNKNOWN ()) 726 0 (727 728 729 730) () 0 () () () 0 0)
731 'set_halo_receives' 'halos_base' 'set_halo_receives' 1 ((PROCEDURE
UNKNOWN-INTENT MODULE-PROC DECL UNKNOWN 0 0 SUBROUTINE) (UNKNOWN 0 0 0
UNKNOWN ()) 732 0 (733 734 735) () 0 () () () 0 0)
736 'set_halo_send' 'halos_base' 'set_halo_send' 1 ((PROCEDURE
UNKNOWN-INTENT MODULE-PROC DECL UNKNOWN 0 0 SUBROUTINE) (UNKNOWN 0 0 0
UNKNOWN ()) 737 0 (738 739 740 741) () 0 () () () 0 0)
742 'set_halo_sends' 'halos_base' 'set_halo_sends' 1 ((PROCEDURE
UNKNOWN-INTENT MODULE-PROC DECL UNKNOWN 0 0 SUBROUTINE) (UNKNOWN 0 0 0
UNKNOWN ()) 743 0 (744 745 746) () 0 () () () 0 0)
747 'set_halo_universal_number' 'halos_numbering'
'set_halo_universal_number' 1 ((PROCEDURE UNKNOWN-INTENT MODULE-PROC
DECL UNKNOWN 0 0 SUBROUTINE ALWAYS_EXPLICIT) (UNKNOWN 0 0 0 UNKNOWN ())
748 0 (749 750 751 752) () 0 () () () 0 0)
753 'trailing_receives_consistent' 'halos_debug'
'trailing_receives_consistent' 1 ((PROCEDURE UNKNOWN-INTENT MODULE-PROC
DECL UNKNOWN 0 0 FUNCTION) (LOGICAL 4 0 0 LOGICAL ()) 754 0 (755) () 756
() () () 0 0)
757 'universal_numbering_count' 'halos_numbering'
'universal_numbering_count' 1 ((PROCEDURE UNKNOWN-INTENT MODULE-PROC
DECL UNKNOWN 0 0 FUNCTION PURE) (INTEGER 4 0 0 INTEGER ()) 758 0 (759) ()
760 () () () 0 0)
761 'valid_global_to_universal_numbering' 'halos_numbering'
'valid_global_to_universal_numbering' 1 ((PROCEDURE UNKNOWN-INTENT
MODULE-PROC DECL UNKNOWN 0 0 FUNCTION) (LOGICAL 4 0 0 LOGICAL ()) 762 0
(763) () 764 () () () 0 0)
765 'valid_halo_communicator' 'halos_debug' 'valid_halo_communicator' 1
((PROCEDURE UNKNOWN-INTENT MODULE-PROC DECL UNKNOWN 0 0 FUNCTION) (
LOGICAL 4 0 0 LOGICAL ()) 766 0 (767) () 768 () () () 0 0)
769 'valid_halo_node_counts' 'halos_debug' 'valid_halo_node_counts' 1 (
(PROCEDURE UNKNOWN-INTENT MODULE-PROC DECL UNKNOWN 0 0 FUNCTION) (
LOGICAL 4 0 0 LOGICAL ()) 770 0 (771) () 772 () () () 0 0)
773 'valid_serial_halo' 'halos_debug' 'valid_serial_halo' 1 ((PROCEDURE
UNKNOWN-INTENT MODULE-PROC DECL UNKNOWN 0 0 FUNCTION IMPLICIT_PURE) (
LOGICAL 4 0 0 LOGICAL ()) 774 0 (775) () 776 () () () 0 0)
777 'verify_halos' 'halos_registration' 'verify_halos' 1 ((PROCEDURE
UNKNOWN-INTENT MODULE-PROC DECL UNKNOWN 0 0 SUBROUTINE) (UNKNOWN 0 0 0
UNKNOWN ()) 778 0 (779) () 0 () () () 0 0)
780 'write_halos' 'halos_registration' 'write_halos' 1 ((PROCEDURE
UNKNOWN-INTENT MODULE-PROC DECL UNKNOWN 0 0 SUBROUTINE) (UNKNOWN 0 0 0
UNKNOWN ()) 781 0 (782 783) () 0 () () () 0 0)
784 'write_universal_numbering' 'halos_diagnostics'
'write_universal_numbering' 1 ((PROCEDURE UNKNOWN-INTENT MODULE-PROC
DECL UNKNOWN 0 0 SUBROUTINE) (UNKNOWN 0 0 0 UNKNOWN ()) 785 0 (786 787
788 789) () 0 () () () 0 0)
632 'valid' '' 'valid' 630 ((VARIABLE UNKNOWN-INTENT UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0 RESULT) (LOGICAL 4 0 0 LOGICAL ()) 0 0 () () 0 () ()
() 0 0)
677 'priority' '' 'priority' 675 ((VARIABLE IN UNKNOWN-PROC UNKNOWN
UNKNOWN 0 0 DUMMY) (INTEGER 4 0 0 INTEGER ()) 0 0 () () 0 () () () 0 0)
755 'halo' '' 'halo' 754 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DUMMY) (DERIVED 144 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
431 'halo' '' 'halo' 430 ((VARIABLE INOUT UNKNOWN-PROC UNKNOWN UNKNOWN 0
0 DUMMY) (DERIVED 144 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
99 'halo' '' 'halo' 98 ((VARIABLE OUT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DUMMY) (DERIVED 144 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
100 'nsends' '' 'nsends' 98 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0
0 DIMENSION DUMMY) (INTEGER 4 0 0 INTEGER ()) 0 0 () (1 0 ASSUMED_SHAPE
(CONSTANT (INTEGER 4 0 0 INTEGER ()) 0 '1') ()) 0 () () () 0 0)
103 'communicator' '' 'communicator' 98 ((VARIABLE IN UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0 OPTIONAL DUMMY) (INTEGER 4 0 0 INTEGER ()) 0 0 () ()
0 () () () 0 0)
106 'data_type' '' 'data_type' 98 ((VARIABLE IN UNKNOWN-PROC UNKNOWN
UNKNOWN 0 0 OPTIONAL DUMMY) (INTEGER 4 0 0 INTEGER ()) 0 0 () () 0 () ()
() 0 0)
105 'nowned_nodes' '' 'nowned_nodes' 98 ((VARIABLE IN UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0 OPTIONAL DUMMY) (INTEGER 4 0 0 INTEGER ()) 0 0 () ()
0 () () () 0 0)
104 'nprocs' '' 'nprocs' 98 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0
0 OPTIONAL DUMMY) (INTEGER 4 0 0 INTEGER ()) 0 0 () () 0 () () () 0 0)
102 'name' '' 'name' 98 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
OPTIONAL DUMMY) (CHARACTER 1 0 0 CHARACTER (())) 0 0 () () 0 () () () 0
0)
101 'nreceives' '' 'nreceives' 98 ((VARIABLE IN UNKNOWN-PROC UNKNOWN
UNKNOWN 0 0 DIMENSION DUMMY) (INTEGER 4 0 0 INTEGER ()) 0 0 () (1 0
ASSUMED_SHAPE (CONSTANT (INTEGER 4 0 0 INTEGER ()) 0 '1') ()) 0 () () ()
0 0)
88 'nsends' '' 'nsends' 86 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0
0 DIMENSION OPTIONAL DUMMY) (INTEGER 4 0 0 INTEGER ()) 0 0 () (1 0
EXPLICIT (CONSTANT (INTEGER 4 0 0 INTEGER ()) 0 '1') (FUNCTION (INTEGER
4 0 0 INTEGER ()) 0 391 (('' (VARIABLE (DERIVED 144 0 0 DERIVED ()) 0 87
()))) 'halo_proc_count' 1 391)) 0 () () () 0 0)
87 'halo' '' 'halo' 86 ((VARIABLE INOUT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DUMMY) (DERIVED 144 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
68 'base_halo' '' 'base_halo' 66 ((VARIABLE IN UNKNOWN-PROC UNKNOWN
UNKNOWN 0 0 DUMMY) (DERIVED 144 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
67 'output_halo' '' 'output_halo' 66 ((VARIABLE OUT UNKNOWN-PROC UNKNOWN
UNKNOWN 0 0 DUMMY) (DERIVED 144 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
70 'halo' '' 'halo' 69 ((VARIABLE INOUT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DUMMY) (DERIVED 144 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
89 'nreceives' '' 'nreceives' 86 ((VARIABLE IN UNKNOWN-PROC UNKNOWN
UNKNOWN 0 0 DIMENSION OPTIONAL DUMMY) (INTEGER 4 0 0 INTEGER ()) 0 0 ()
(1 0 EXPLICIT (CONSTANT (INTEGER 4 0 0 INTEGER ()) 0 '1') (FUNCTION (
INTEGER 4 0 0 INTEGER ()) 0 391 (('' (VARIABLE (DERIVED 144 0 0 DERIVED
()) 0 87 ()))) 'halo_proc_count' 1 391)) 0 () () () 0 0)
421 'halos' '' 'halos' 420 ((VARIABLE INOUT UNKNOWN-PROC UNKNOWN UNKNOWN
0 0 DIMENSION DUMMY) (DERIVED 144 0 0 DERIVED ()) 0 0 () (1 0
ASSUMED_SHAPE (CONSTANT (INTEGER 4 0 0 INTEGER ()) 0 '1') ()) 0 () () ()
0 0)
91 'halo' '' 'halo' 90 ((VARIABLE INOUT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DUMMY) (DERIVED 144 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
72 'object' '' 'object' 71 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0
0 TARGET DUMMY) (DERIVED 144 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
74 'object' '' 'object' 73 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0
0 DUMMY) (DERIVED 144 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
75 'has_references' '' 'has_references' 73 ((VARIABLE UNKNOWN-INTENT
UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 RESULT) (LOGICAL 4 0 0 LOGICAL ()) 0 0
() () 0 () () () 0 0)
107 'ordering_scheme' '' 'ordering_scheme' 98 ((VARIABLE IN UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0 OPTIONAL DUMMY) (INTEGER 4 0 0 INTEGER ()) 0 0 () ()
0 () () () 0 0)
434 'halo' '' 'halo' 433 ((VARIABLE INOUT UNKNOWN-PROC UNKNOWN UNKNOWN 0
0 DUMMY) (DERIVED 144 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
393 'halo' '' 'halo' 392 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DUMMY) (DERIVED 144 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
756 'consistent' '' 'consistent' 754 ((VARIABLE UNKNOWN-INTENT
UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 RESULT) (LOGICAL 4 0 0 LOGICAL ()) 0 0
() () 0 () () () 0 0)
772 'valid' '' 'valid' 770 ((VARIABLE UNKNOWN-INTENT UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0 RESULT) (LOGICAL 4 0 0 LOGICAL ()) 0 0 () () 0 () ()
() 0 0)
771 'halo' '' 'halo' 770 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DUMMY) (DERIVED 144 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
776 'valid' '' 'valid' 774 ((VARIABLE UNKNOWN-INTENT UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0 RESULT) (LOGICAL 4 0 0 LOGICAL ()) 0 0 () () 0 () ()
() 0 0)
775 'halo' '' 'halo' 774 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DUMMY) (DERIVED 144 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
97 'pending' '' 'pending' 95 ((VARIABLE UNKNOWN-INTENT UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0 RESULT) (LOGICAL 4 0 0 LOGICAL ()) 0 0 () () 0 () ()
() 0 0)
96 'halo' '' 'halo' 95 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DUMMY) (DERIVED 144 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
768 'valid' '' 'valid' 766 ((VARIABLE UNKNOWN-INTENT UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0 RESULT) (LOGICAL 4 0 0 LOGICAL ()) 0 0 () () 0 () ()
() 0 0)
94 'pending' '' 'pending' 92 ((VARIABLE UNKNOWN-INTENT UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0 RESULT ALWAYS_EXPLICIT) (LOGICAL 4 0 0 LOGICAL ()) 0
0 () () 0 () () () 0 0)
93 'communicator' '' 'communicator' 92 ((VARIABLE IN UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0 OPTIONAL DUMMY) (INTEGER 4 0 0 INTEGER ()) 0 0 () ()
0 () () () 0 0)
767 'halo' '' 'halo' 766 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DUMMY) (DERIVED 144 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
631 'halo' '' 'halo' 630 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DUMMY) (DERIVED 144 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
676 'halo' '' 'halo' 675 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DUMMY) (DERIVED 144 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
471 'start_indices' '' 'start_indices' 467 ((VARIABLE OUT UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0 DIMENSION OPTIONAL DUMMY) (INTEGER 4 0 0 INTEGER ())
0 0 () (1 0 EXPLICIT (CONSTANT (INTEGER 4 0 0 INTEGER ()) 0 '1') (
FUNCTION (INTEGER 4 0 0 INTEGER ()) 0 391 (('' (VARIABLE (DERIVED 144 0
0 DERIVED ()) 0 468 ()))) 'halo_proc_count' 1 391)) 0 () () () 0 0)
476 'nsends' '' 'nsends' 473 ((VARIABLE OUT UNKNOWN-PROC UNKNOWN UNKNOWN
0 0 DIMENSION OPTIONAL DUMMY) (INTEGER 4 0 0 INTEGER ()) 0 0 () (1 0
EXPLICIT (CONSTANT (INTEGER 4 0 0 INTEGER ()) 0 '1') (FUNCTION (INTEGER
4 0 0 INTEGER ()) 0 391 (('' (VARIABLE (DERIVED 144 0 0 DERIVED ()) 0
474 ()))) 'halo_proc_count' 1 391)) 0 () () () 0 0)
475 'sends' '' 'sends' 473 ((VARIABLE OUT UNKNOWN-PROC UNKNOWN UNKNOWN 0
0 DIMENSION DUMMY) (INTEGER 4 0 0 INTEGER ()) 0 0 () (1 0 EXPLICIT (
CONSTANT (INTEGER 4 0 0 INTEGER ()) 0 '1') (FUNCTION (INTEGER 4 0 0
INTEGER ()) 0 521 (('' (VARIABLE (DERIVED 144 0 0 DERIVED ()) 0 474 ())))
'halo_all_sends_count' 1 521)) 0 () () () 0 0)
474 'halo' '' 'halo' 473 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DUMMY) (DERIVED 144 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
470 'nreceives' '' 'nreceives' 467 ((VARIABLE OUT UNKNOWN-PROC UNKNOWN
UNKNOWN 0 0 DIMENSION OPTIONAL DUMMY) (INTEGER 4 0 0 INTEGER ()) 0 0 ()
(1 0 EXPLICIT (CONSTANT (INTEGER 4 0 0 INTEGER ()) 0 '1') (FUNCTION (
INTEGER 4 0 0 INTEGER ()) 0 391 (('' (VARIABLE (DERIVED 144 0 0 DERIVED
()) 0 468 ()))) 'halo_proc_count' 1 391)) 0 () () () 0 0)
523 'halo' '' 'halo' 522 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DUMMY) (DERIVED 144 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
477 'start_indices' '' 'start_indices' 473 ((VARIABLE OUT UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0 DIMENSION OPTIONAL DUMMY) (INTEGER 4 0 0 INTEGER ())
0 0 () (1 0 EXPLICIT (CONSTANT (INTEGER 4 0 0 INTEGER ()) 0 '1') (
FUNCTION (INTEGER 4 0 0 INTEGER ()) 0 391 (('' (VARIABLE (DERIVED 144 0
0 DERIVED ()) 0 474 ()))) 'halo_proc_count' 1 391)) 0 () () () 0 0)
526 'halo' '' 'halo' 525 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DUMMY) (DERIVED 144 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
531 'len_trim' '(intrinsic)' 'len_trim' 534 ((PROCEDURE UNKNOWN-INTENT
UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 FUNCTION) (INTEGER 4 0 0 INTEGER ()) 0
0 () () 531 () () () 0 0)
561 'process' '' 'process' 559 ((VARIABLE IN UNKNOWN-PROC UNKNOWN
UNKNOWN 0 0 DUMMY) (INTEGER 4 0 0 INTEGER ()) 0 0 () () 0 () () () 0 0)
560 'halo' '' 'halo' 559 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DUMMY) (DERIVED 144 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
562 'index' '' 'index' 559 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0
0 DUMMY) (INTEGER 4 0 0 INTEGER ()) 0 0 () () 0 () () () 0 0)
555 'halo' '' 'halo' 554 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DUMMY) (DERIVED 144 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
565 'halo' '' 'halo' 564 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DUMMY) (DERIVED 144 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
566 'process' '' 'process' 564 ((VARIABLE IN UNKNOWN-PROC UNKNOWN
UNKNOWN 0 0 DUMMY) (INTEGER 4 0 0 INTEGER ()) 0 0 () () 0 () () () 0 0)
550 'halo' '' 'halo' 549 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DUMMY) (DERIVED 144 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
574 'process' '' 'process' 572 ((VARIABLE IN UNKNOWN-PROC UNKNOWN
UNKNOWN 0 0 DUMMY) (INTEGER 4 0 0 INTEGER ()) 0 0 () () 0 () () () 0 0)
579 'index' '' 'index' 576 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0
0 DUMMY) (INTEGER 4 0 0 INTEGER ()) 0 0 () () 0 () () () 0 0)
578 'process' '' 'process' 576 ((VARIABLE IN UNKNOWN-PROC UNKNOWN
UNKNOWN 0 0 DUMMY) (INTEGER 4 0 0 INTEGER ()) 0 0 () () 0 () () () 0 0)
582 'halo' '' 'halo' 581 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DUMMY) (DERIVED 144 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
583 'process' '' 'process' 581 ((VARIABLE IN UNKNOWN-PROC UNKNOWN
UNKNOWN 0 0 DUMMY) (INTEGER 4 0 0 INTEGER ()) 0 0 () () 0 () () () 0 0)
591 'process' '' 'process' 589 ((VARIABLE IN UNKNOWN-PROC UNKNOWN
UNKNOWN 0 0 DUMMY) (INTEGER 4 0 0 INTEGER ()) 0 0 () () 0 () () () 0 0)
577 'halo' '' 'halo' 576 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DUMMY) (DERIVED 144 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
612 'process' '' 'process' 610 ((VARIABLE IN UNKNOWN-PROC UNKNOWN
UNKNOWN 0 0 DUMMY) (INTEGER 4 0 0 INTEGER ()) 0 0 () () 0 () () () 0 0)
656 'halo' '' 'halo' 655 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DUMMY) (DERIVED 144 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
640 'halo' '' 'halo' 639 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DUMMY) (DERIVED 144 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
611 'halo' '' 'halo' 610 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DUMMY) (DERIVED 144 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
573 'halo' '' 'halo' 572 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DUMMY) (DERIVED 144 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
664 'halo' '' 'halo' 663 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DUMMY) (DERIVED 144 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
665 'min_node' '' 'min_node' 663 ((VARIABLE UNKNOWN-INTENT UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0 RESULT) (INTEGER 4 0 0 INTEGER ()) 0 0 () () 0 () ()
() 0 0)
661 'max_node' '' 'max_node' 659 ((VARIABLE UNKNOWN-INTENT UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0 RESULT) (INTEGER 4 0 0 INTEGER ()) 0 0 () () 0 () ()
() 0 0)
700 'receives' '' 'receives' 698 ((VARIABLE IN UNKNOWN-PROC UNKNOWN
UNKNOWN 0 0 DIMENSION DUMMY) (INTEGER 4 0 0 INTEGER ()) 0 0 () (1 0
ASSUMED_SHAPE (CONSTANT (INTEGER 4 0 0 INTEGER ()) 0 '1') ()) 0 () () ()
0 0)
704 'sends' '' 'sends' 702 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0
0 DIMENSION DUMMY) (INTEGER 4 0 0 INTEGER ()) 0 0 () (1 0 ASSUMED_SHAPE
(CONSTANT (INTEGER 4 0 0 INTEGER ()) 0 '1') ()) 0 () () () 0 0)
703 'halo' '' 'halo' 702 ((VARIABLE INOUT UNKNOWN-PROC UNKNOWN UNKNOWN 0
0 DUMMY) (DERIVED 144 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
708 'communicator' '' 'communicator' 706 ((VARIABLE IN UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0 DUMMY) (INTEGER 4 0 0 INTEGER ()) 0 0 () () 0 () ()
() 0 0)
707 'halo' '' 'halo' 706 ((VARIABLE INOUT UNKNOWN-PROC UNKNOWN UNKNOWN 0
0 DUMMY) (DERIVED 144 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
728 'process' '' 'process' 726 ((VARIABLE IN UNKNOWN-PROC UNKNOWN
UNKNOWN 0 0 DUMMY) (INTEGER 4 0 0 INTEGER ()) 0 0 () () 0 () () () 0 0)
727 'halo' '' 'halo' 726 ((VARIABLE INOUT UNKNOWN-PROC UNKNOWN UNKNOWN 0
0 DUMMY) (DERIVED 144 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
724 'ordering_scheme' '' 'ordering_scheme' 722 ((VARIABLE IN
UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 DUMMY) (INTEGER 4 0 0 INTEGER ()) 0 0 ()
() 0 () () () 0 0)
723 'halo' '' 'halo' 722 ((VARIABLE INOUT UNKNOWN-PROC UNKNOWN UNKNOWN 0
0 DUMMY) (DERIVED 144 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
733 'halo' '' 'halo' 732 ((VARIABLE INOUT UNKNOWN-PROC UNKNOWN UNKNOWN 0
0 DUMMY) (DERIVED 144 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
735 'receives' '' 'receives' 732 ((VARIABLE IN UNKNOWN-PROC UNKNOWN
UNKNOWN 0 0 DIMENSION DUMMY) (INTEGER 4 0 0 INTEGER ()) 0 0 () (1 0
EXPLICIT (CONSTANT (INTEGER 4 0 0 INTEGER ()) 0 '1') (FUNCTION (INTEGER
4 0 0 INTEGER ()) 0 563 (('' (VARIABLE (DERIVED 144 0 0 DERIVED ()) 0
733 ())) ('' (VARIABLE (INTEGER 4 0 0 INTEGER ()) 0 734 ())))
'halo_receive_count' 1 563)) 0 () () () 0 0)
734 'process' '' 'process' 732 ((VARIABLE IN UNKNOWN-PROC UNKNOWN
UNKNOWN 0 0 DUMMY) (INTEGER 4 0 0 INTEGER ()) 0 0 () () 0 () () () 0 0)
745 'process' '' 'process' 743 ((VARIABLE IN UNKNOWN-PROC UNKNOWN
UNKNOWN 0 0 DUMMY) (INTEGER 4 0 0 INTEGER ()) 0 0 () () 0 () () () 0 0)
419 'halo' '' 'halo' 418 ((VARIABLE INOUT UNKNOWN-PROC UNKNOWN UNKNOWN 0
0 DUMMY) (DERIVED 144 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
746 'sends' '' 'sends' 743 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0
0 DIMENSION DUMMY) (INTEGER 4 0 0 INTEGER ()) 0 0 () (1 0 EXPLICIT (
CONSTANT (INTEGER 4 0 0 INTEGER ()) 0 '1') (FUNCTION (INTEGER 4 0 0
INTEGER ()) 0 580 (('' (VARIABLE (DERIVED 144 0 0 DERIVED ()) 0 744 ()))
('' (VARIABLE (INTEGER 4 0 0 INTEGER ()) 0 745 ()))) 'halo_send_count' 1
580)) 0 () () () 0 0)
411 'communicator' '' 'communicator' 409 ((VARIABLE UNKNOWN-INTENT
UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 RESULT) (INTEGER 4 0 0 INTEGER ()) 0 0
() () 0 () () () 0 0)
410 'halo' '' 'halo' 409 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DUMMY) (DERIVED 144 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
744 'halo' '' 'halo' 743 ((VARIABLE INOUT UNKNOWN-PROC UNKNOWN UNKNOWN 0
0 DUMMY) (DERIVED 144 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
730 'node' '' 'node' 726 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DUMMY) (INTEGER 4 0 0 INTEGER ()) 0 0 () () 0 () () () 0 0)
729 'index' '' 'index' 726 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0
0 DUMMY) (INTEGER 4 0 0 INTEGER ()) 0 0 () () 0 () () () 0 0)
382 'node_count' '' 'node_count' 380 ((VARIABLE UNKNOWN-INTENT
UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 RESULT) (INTEGER 4 0 0 INTEGER ()) 0 0
() () 0 () () () 0 0)
396 'serial' '' 'serial' 394 ((VARIABLE UNKNOWN-INTENT UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0 RESULT) (LOGICAL 4 0 0 LOGICAL ()) 0 0 () () 0 () ()
() 0 0)
395 'halo' '' 'halo' 394 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DUMMY) (DERIVED 144 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
406 'halos' '' 'halos' 405 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0
0 DIMENSION DUMMY) (DERIVED 144 0 0 DERIVED ()) 0 0 () (1 0
ASSUMED_SHAPE (CONSTANT (INTEGER 4 0 0 INTEGER ()) 0 '1') ()) 0 () () ()
0 0)
407 'size' '(intrinsic)' 'size' 405 ((PROCEDURE UNKNOWN-INTENT
UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 FUNCTION) (REAL 8 0 0 REAL ()) 0 0 () ()
407 () () () 0 0)
408 'serial' '' 'serial' 405 ((VARIABLE UNKNOWN-INTENT UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0 DIMENSION RESULT ALWAYS_EXPLICIT) (LOGICAL 4 0 0
LOGICAL ()) 0 0 () (1 0 EXPLICIT (CONSTANT (INTEGER 4 0 0 INTEGER ()) 0
'1') (FUNCTION (INTEGER 4 0 0 INTEGER ()) 0 407 (('' (VARIABLE (DERIVED
144 0 0 DERIVED ()) 1 406 ((ARRAY (FULL 0))))) ('' ()) ('' ())) '' 0
'size')) 0 () () () 0 0)
381 'halo' '' 'halo' 380 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DUMMY) (DERIVED 144 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
657 'max_node' '' 'max_node' 655 ((VARIABLE UNKNOWN-INTENT UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0 RESULT) (INTEGER 4 0 0 INTEGER ()) 0 0 () () 0 () ()
() 0 0)
532 'halo' '' 'halo' 534 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DUMMY) (DERIVED 144 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
468 'halo' '' 'halo' 467 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DUMMY) (DERIVED 144 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
469 'receives' '' 'receives' 467 ((VARIABLE OUT UNKNOWN-PROC UNKNOWN
UNKNOWN 0 0 DIMENSION DUMMY) (INTEGER 4 0 0 INTEGER ()) 0 0 () (1 0
EXPLICIT (CONSTANT (INTEGER 4 0 0 INTEGER ()) 0 '1') (FUNCTION (INTEGER
4 0 0 INTEGER ()) 0 518 (('' (VARIABLE (DERIVED 144 0 0 DERIVED ()) 0
468 ()))) 'halo_all_receives_count' 1 518)) 0 () () () 0 0)
520 'halo' '' 'halo' 519 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DUMMY) (DERIVED 144 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
529 'halo' '' 'halo' 528 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DUMMY) (DERIVED 144 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
569 'halo' '' 'halo' 568 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DUMMY) (DERIVED 144 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
570 'nreceives' '' 'nreceives' 568 ((VARIABLE OUT UNKNOWN-PROC UNKNOWN
UNKNOWN 0 0 DIMENSION DUMMY) (INTEGER 4 0 0 INTEGER ()) 0 0 () (1 0
EXPLICIT (CONSTANT (INTEGER 4 0 0 INTEGER ()) 0 '1') (FUNCTION (INTEGER
4 0 0 INTEGER ()) 0 391 (('' (VARIABLE (DERIVED 144 0 0 DERIVED ()) 0
569 ()))) 'halo_proc_count' 1 391)) 0 () () () 0 0)
586 'halo' '' 'halo' 585 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DUMMY) (DERIVED 144 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
587 'nsends' '' 'nsends' 585 ((VARIABLE OUT UNKNOWN-PROC UNKNOWN UNKNOWN
0 0 DIMENSION DUMMY) (INTEGER 4 0 0 INTEGER ()) 0 0 () (1 0 EXPLICIT (
CONSTANT (INTEGER 4 0 0 INTEGER ()) 0 '1') (FUNCTION (INTEGER 4 0 0
INTEGER ()) 0 391 (('' (VARIABLE (DERIVED 144 0 0 DERIVED ()) 0 586 ())))
'halo_proc_count' 1 391)) 0 () () () 0 0)
590 'halo' '' 'halo' 589 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DUMMY) (DERIVED 144 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
652 'halo' '' 'halo' 651 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DUMMY) (DERIVED 144 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
653 'max_node' '' 'max_node' 651 ((VARIABLE UNKNOWN-INTENT UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0 RESULT) (INTEGER 4 0 0 INTEGER ()) 0 0 () () 0 () ()
() 0 0)
660 'halo' '' 'halo' 659 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DUMMY) (DERIVED 144 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
669 'min_node' '' 'min_node' 667 ((VARIABLE UNKNOWN-INTENT UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0 RESULT) (INTEGER 4 0 0 INTEGER ()) 0 0 () () 0 () ()
() 0 0)
668 'halo' '' 'halo' 667 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DUMMY) (DERIVED 144 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
673 'min_node' '' 'min_node' 671 ((VARIABLE UNKNOWN-INTENT UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0 RESULT) (INTEGER 4 0 0 INTEGER ()) 0 0 () () 0 () ()
() 0 0)
672 'halo' '' 'halo' 671 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DUMMY) (DERIVED 144 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
699 'halo' '' 'halo' 698 ((VARIABLE INOUT UNKNOWN-PROC UNKNOWN UNKNOWN 0
0 DUMMY) (DERIVED 144 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
712 'data_type' '' 'data_type' 710 ((VARIABLE IN UNKNOWN-PROC UNKNOWN
UNKNOWN 0 0 DUMMY) (INTEGER 4 0 0 INTEGER ()) 0 0 () () 0 () () () 0 0)
716 'name' '' 'name' 714 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DUMMY) (CHARACTER 1 0 0 CHARACTER (())) 0 0 () () 0 () () () 0 0)
715 'halo' '' 'halo' 714 ((VARIABLE INOUT UNKNOWN-PROC UNKNOWN UNKNOWN 0
0 DUMMY) (DERIVED 144 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
711 'halo' '' 'halo' 710 ((VARIABLE INOUT UNKNOWN-PROC UNKNOWN UNKNOWN 0
0 DUMMY) (DERIVED 144 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
719 'halo' '' 'halo' 718 ((VARIABLE INOUT UNKNOWN-PROC UNKNOWN UNKNOWN 0
0 DUMMY) (DERIVED 144 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
720 'nowned_nodes' '' 'nowned_nodes' 718 ((VARIABLE IN UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0 DUMMY) (INTEGER 4 0 0 INTEGER ()) 0 0 () () 0 () ()
() 0 0)
738 'halo' '' 'halo' 737 ((VARIABLE INOUT UNKNOWN-PROC UNKNOWN UNKNOWN 0
0 DUMMY) (DERIVED 144 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
740 'index' '' 'index' 737 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0
0 DUMMY) (INTEGER 4 0 0 INTEGER ()) 0 0 () () 0 () () () 0 0)
741 'node' '' 'node' 737 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DUMMY) (INTEGER 4 0 0 INTEGER ()) 0 0 () () 0 () () () 0 0)
739 'process' '' 'process' 737 ((VARIABLE IN UNKNOWN-PROC UNKNOWN
UNKNOWN 0 0 DUMMY) (INTEGER 4 0 0 INTEGER ()) 0 0 () () 0 () () () 0 0)
789 'name' '' 'name' 785 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DUMMY) (CHARACTER 1 0 0 CHARACTER (())) 0 0 () () 0 () () () 0 0)
788 'position' '' 'position' 785 ((VARIABLE IN UNKNOWN-PROC UNKNOWN
UNKNOWN 0 0 DUMMY) (DERIVED 108 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
787 'mesh' '' 'mesh' 785 ((VARIABLE INOUT UNKNOWN-PROC UNKNOWN UNKNOWN 0
0 DUMMY) (DERIVED 118 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
786 'halo' '' 'halo' 785 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DUMMY) (DERIVED 144 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
439 'create_caches' '' 'create_caches' 436 ((VARIABLE IN UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0 OPTIONAL DUMMY) (LOGICAL 4 0 0 LOGICAL ()) 0 0 () ()
0 () () () 0 0)
438 'ordering_scheme' '' 'ordering_scheme' 436 ((VARIABLE IN
UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 OPTIONAL DUMMY) (INTEGER 4 0 0 INTEGER
()) 0 0 () () 0 () () () 0 0)
442 'mesh' '' 'mesh' 441 ((VARIABLE INOUT UNKNOWN-PROC UNKNOWN UNKNOWN 0
0 DUMMY) (DERIVED 118 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
443 'element_halo' '' 'element_halo' 441 ((VARIABLE IN UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0 DUMMY) (DERIVED 144 0 0 DERIVED ()) 0 0 () () 0 () ()
() 0 0)
445 'create_caches' '' 'create_caches' 441 ((VARIABLE IN UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0 OPTIONAL DUMMY) (LOGICAL 4 0 0 LOGICAL ()) 0 0 () ()
0 () () () 0 0)
449 'new_positions' '' 'new_positions' 448 ((VARIABLE INOUT UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0 DUMMY) (DERIVED 108 0 0 DERIVED ()) 0 0 () () 0 () ()
() 0 0)
446 'selement_halo' '' 'selement_halo' 441 ((VARIABLE UNKNOWN-INTENT
UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 RESULT ALWAYS_EXPLICIT) (DERIVED 144 0
0 DERIVED ()) 0 0 () () 0 () () () 0 0)
444 'ordering_scheme' '' 'ordering_scheme' 441 ((VARIABLE IN
UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 OPTIONAL DUMMY) (INTEGER 4 0 0 INTEGER
()) 0 0 () () 0 () () () 0 0)
451 'aliased_to_new_node_number' '' 'aliased_to_new_node_number' 448 ((
VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 DUMMY) (DERIVED 288 0 0
DERIVED ()) 0 0 () () 0 () () () 0 0)
461 'node_halo' '' 'node_halo' 458 ((VARIABLE IN UNKNOWN-PROC UNKNOWN
UNKNOWN 0 0 DUMMY) (DERIVED 144 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
460 'mesh' '' 'mesh' 458 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DUMMY) (DERIVED 118 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
649 'unknowns_sizes' '' 'unknowns_sizes' 645 ((VARIABLE UNKNOWN-INTENT
UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 DIMENSION RESULT ALWAYS_EXPLICIT) (
INTEGER 4 0 0 INTEGER ()) 0 0 () (1 0 EXPLICIT (CONSTANT (INTEGER 4 0 0
INTEGER ()) 0 '1') (FUNCTION (INTEGER 4 0 0 INTEGER ()) 0 648 (('' (
VARIABLE (INTEGER 4 0 0 INTEGER ()) 1 646 ((ARRAY (FULL 0))))) ('' ()) (
'' ())) '' 0 'size')) 0 () () () 0 0)
648 'size' '(intrinsic)' 'size' 645 ((PROCEDURE UNKNOWN-INTENT
UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 FUNCTION) (REAL 8 0 0 REAL ()) 0 0 () ()
648 () () () 0 0)
647 'communicator' '' 'communicator' 645 ((VARIABLE IN UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0 OPTIONAL DUMMY) (INTEGER 4 0 0 INTEGER ()) 0 0 () ()
0 () () () 0 0)
646 'knowns_sizes' '' 'knowns_sizes' 645 ((VARIABLE IN UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0 DIMENSION DUMMY) (INTEGER 4 0 0 INTEGER ()) 0 0 () (
1 0 ASSUMED_SHAPE (CONSTANT (INTEGER 4 0 0 INTEGER ()) 0 '1') ()) 0 () ()
() 0 0)
450 'model' '' 'model' 448 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0
0 DUMMY) (DERIVED 108 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
285 'mesh' '' 'mesh' 284 ((VARIABLE INOUT UNKNOWN-PROC UNKNOWN UNKNOWN 0
0 DUMMY) (DERIVED 118 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
287 'create_caches' '' 'create_caches' 284 ((VARIABLE IN UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0 OPTIONAL DUMMY) (LOGICAL 4 0 0 LOGICAL ()) 0 0 () ()
0 () () () 0 0)
286 'ordering_scheme' '' 'ordering_scheme' 284 ((VARIABLE IN
UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 OPTIONAL DUMMY) (INTEGER 4 0 0 INTEGER
()) 0 0 () () 0 () () () 0 0)
437 'mesh' '' 'mesh' 436 ((VARIABLE INOUT UNKNOWN-PROC UNKNOWN UNKNOWN 0
0 DUMMY) (DERIVED 118 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
454 'halo' '' 'halo' 453 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DUMMY) (DERIVED 144 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
455 'sub_nodes' '' 'sub_nodes' 453 ((VARIABLE IN UNKNOWN-PROC UNKNOWN
UNKNOWN 0 0 DIMENSION DUMMY) (INTEGER 4 0 0 INTEGER ()) 0 0 () (1 0
ASSUMED_SHAPE (CONSTANT (INTEGER 4 0 0 INTEGER ()) 0 '1') ()) 0 () () ()
0 0)
459 'ele' '' 'ele' 458 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DUMMY) (INTEGER 4 0 0 INTEGER ()) 0 0 () () 0 () () () 0 0)
456 'sub_halo' '' 'sub_halo' 453 ((VARIABLE UNKNOWN-INTENT UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0 RESULT ALWAYS_EXPLICIT) (DERIVED 144 0 0 DERIVED ())
0 0 () () 0 () () () 0 0)
354 'halo' '' 'halo' 353 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DUMMY) (DERIVED 144 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
355 'integer_data' '' 'integer_data' 353 ((VARIABLE INOUT UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0 DIMENSION DUMMY) (INTEGER 4 0 0 INTEGER ()) 0 0 () (
1 0 ASSUMED_SHAPE (CONSTANT (INTEGER 4 0 0 INTEGER ()) 0 '1') ()) 0 () ()
() 0 0)
312 'v_field' '' 'v_field' 310 ((VARIABLE INOUT UNKNOWN-PROC UNKNOWN
UNKNOWN 0 0 DUMMY) (DERIVED 108 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
308 'halo' '' 'halo' 307 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DUMMY) (DERIVED 144 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
296 'v_field' '' 'v_field' 295 ((VARIABLE INOUT UNKNOWN-PROC UNKNOWN
UNKNOWN 0 0 DUMMY) (DERIVED 108 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
300 'level' '' 'level' 298 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0
0 OPTIONAL DUMMY) (INTEGER 4 0 0 INTEGER ()) 0 0 () () 0 () () () 0 0)
299 't_field' '' 't_field' 298 ((VARIABLE INOUT UNKNOWN-PROC UNKNOWN
UNKNOWN 0 0 DUMMY) (DERIVED 359 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
303 'real_data' '' 'real_data' 301 ((VARIABLE INOUT UNKNOWN-PROC UNKNOWN
UNKNOWN 0 0 DIMENSION DUMMY) (REAL 8 0 0 REAL ()) 0 0 () (1 0
ASSUMED_SHAPE (CONSTANT (INTEGER 4 0 0 INTEGER ()) 0 '1') ()) 0 () () ()
0 0)
302 'halo' '' 'halo' 301 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DUMMY) (DERIVED 144 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
297 'level' '' 'level' 295 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0
0 OPTIONAL DUMMY) (INTEGER 4 0 0 INTEGER ()) 0 0 () () 0 () () () 0 0)
306 's_field' '' 's_field' 304 ((VARIABLE INOUT UNKNOWN-PROC UNKNOWN
UNKNOWN 0 0 DUMMY) (DERIVED 256 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
294 'level' '' 'level' 292 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0
0 OPTIONAL DUMMY) (INTEGER 4 0 0 INTEGER ()) 0 0 () () 0 () () () 0 0)
293 's_field' '' 's_field' 292 ((VARIABLE INOUT UNKNOWN-PROC UNKNOWN
UNKNOWN 0 0 DUMMY) (DERIVED 256 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
305 'halo' '' 'halo' 304 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DUMMY) (DERIVED 144 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
309 't_field' '' 't_field' 307 ((VARIABLE INOUT UNKNOWN-PROC UNKNOWN
UNKNOWN 0 0 DUMMY) (DERIVED 359 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
351 'integer_array' '' 'integer_array' 349 ((VARIABLE IN UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0 DIMENSION DUMMY) (INTEGER 4 0 0 INTEGER ()) 0 0 () (
1 0 ASSUMED_SHAPE (CONSTANT (INTEGER 4 0 0 INTEGER ()) 0 '1') ()) 0 () ()
() 0 0)
340 'verifies' '' 'verifies' 337 ((VARIABLE UNKNOWN-INTENT UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0 RESULT ALWAYS_EXPLICIT) (LOGICAL 4 0 0 LOGICAL ()) 0
0 () () 0 () () () 0 0)
339 'real_array' '' 'real_array' 337 ((VARIABLE IN UNKNOWN-PROC UNKNOWN
UNKNOWN 0 0 DIMENSION DUMMY) (REAL 8 0 0 REAL ()) 0 0 () (1 0
ASSUMED_SHAPE (CONSTANT (INTEGER 4 0 0 INTEGER ()) 0 '1') ()) 0 () () ()
0 0)
338 'halo' '' 'halo' 337 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DUMMY) (DERIVED 144 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
352 'verifies' '' 'verifies' 349 ((VARIABLE UNKNOWN-INTENT UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0 RESULT ALWAYS_EXPLICIT) (LOGICAL 4 0 0 LOGICAL ()) 0
0 () () 0 () () () 0 0)
343 'sfield' '' 'sfield' 341 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN
0 0 DUMMY) (DERIVED 256 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
344 'verifies' '' 'verifies' 341 ((VARIABLE UNKNOWN-INTENT UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0 RESULT) (LOGICAL 4 0 0 LOGICAL ()) 0 0 () () 0 () ()
() 0 0)
342 'halo' '' 'halo' 341 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DUMMY) (DERIVED 144 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
350 'halo' '' 'halo' 349 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DUMMY) (DERIVED 144 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
335 'dim' '' 'dim' 332 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DUMMY) (INTEGER 4 0 0 INTEGER ()) 0 0 () () 0 () () () 0 0)
334 'vfield' '' 'vfield' 332 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN
0 0 DUMMY) (DERIVED 108 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
346 'halo' '' 'halo' 345 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DUMMY) (DERIVED 144 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
348 'verifies' '' 'verifies' 345 ((VARIABLE UNKNOWN-INTENT UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0 RESULT) (LOGICAL 4 0 0 LOGICAL ()) 0 0 () () 0 () ()
() 0 0)
347 'vfield' '' 'vfield' 345 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN
0 0 DUMMY) (DERIVED 108 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
336 'verifies' '' 'verifies' 332 ((VARIABLE UNKNOWN-INTENT UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0 RESULT) (LOGICAL 4 0 0 LOGICAL ()) 0 0 () () 0 () ()
() 0 0)
333 'halo' '' 'halo' 332 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DUMMY) (DERIVED 144 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
311 'halo' '' 'halo' 310 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DUMMY) (DERIVED 144 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
330 'halo' '' 'halo' 329 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DUMMY) (DERIVED 144 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
331 'integer_data' '' 'integer_data' 329 ((VARIABLE INOUT UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0 DIMENSION DUMMY) (INTEGER 4 0 0 INTEGER ()) 0 0 () (
2 0 ASSUMED_SHAPE (CONSTANT (INTEGER 4 0 0 INTEGER ()) 0 '1') () (
CONSTANT (INTEGER 4 0 0 INTEGER ()) 0 '1') ()) 0 () () () 0 0)
358 'integer_data' '' 'integer_data' 356 ((VARIABLE INOUT UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0 DIMENSION DUMMY) (INTEGER 4 0 0 INTEGER ()) 0 0 () (
3 0 ASSUMED_SHAPE (CONSTANT (INTEGER 4 0 0 INTEGER ()) 0 '1') () (
CONSTANT (INTEGER 4 0 0 INTEGER ()) 0 '1') () (CONSTANT (INTEGER 4 0 0
INTEGER ()) 0 '1') ()) 0 () () () 0 0)
326 'halo' '' 'halo' 325 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DUMMY) (DERIVED 144 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
328 'block_size' '' 'block_size' 325 ((VARIABLE IN UNKNOWN-PROC UNKNOWN
UNKNOWN 0 0 DUMMY) (INTEGER 4 0 0 INTEGER ()) 0 0 () () 0 () () () 0 0)
327 'integer_data' '' 'integer_data' 325 ((VARIABLE INOUT UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0 DIMENSION DUMMY) (INTEGER 4 0 0 INTEGER ()) 0 0 () (
1 0 ASSUMED_SIZE (CONSTANT (INTEGER 4 0 0 INTEGER ()) 0 '1') ()) 0 () ()
() 0 0)
318 'real_data' '' 'real_data' 316 ((VARIABLE INOUT UNKNOWN-PROC UNKNOWN
UNKNOWN 0 0 DIMENSION DUMMY) (REAL 8 0 0 REAL ()) 0 0 () (1 0
ASSUMED_SHAPE (CONSTANT (INTEGER 4 0 0 INTEGER ()) 0 '1') ()) 0 () () ()
0 0)
320 'halo' '' 'halo' 319 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DUMMY) (DERIVED 144 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
315 'real_data' '' 'real_data' 313 ((VARIABLE INOUT UNKNOWN-PROC UNKNOWN
UNKNOWN 0 0 DIMENSION DUMMY) (REAL 8 0 0 REAL ()) 0 0 () (2 0
ASSUMED_SHAPE (CONSTANT (INTEGER 4 0 0 INTEGER ()) 0 '1') () (CONSTANT (
INTEGER 4 0 0 INTEGER ()) 0 '1') ()) 0 () () () 0 0)
314 'halo' '' 'halo' 313 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DUMMY) (DERIVED 144 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
317 'halo' '' 'halo' 316 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DUMMY) (DERIVED 144 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
357 'halo' '' 'halo' 356 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DUMMY) (DERIVED 144 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
324 's_field' '' 's_field' 322 ((VARIABLE INOUT UNKNOWN-PROC UNKNOWN
UNKNOWN 0 0 DUMMY) (DERIVED 256 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
323 'halo' '' 'halo' 322 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DUMMY) (DERIVED 144 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
321 'real_data' '' 'real_data' 319 ((VARIABLE INOUT UNKNOWN-PROC UNKNOWN
UNKNOWN 0 0 DIMENSION DUMMY) (REAL 8 0 0 REAL ()) 0 0 () (3 0
ASSUMED_SHAPE (CONSTANT (INTEGER 4 0 0 INTEGER ()) 0 '1') () (CONSTANT (
INTEGER 4 0 0 INTEGER ()) 0 '1') () (CONSTANT (INTEGER 4 0 0 INTEGER ())
0 '1') ()) 0 () () () 0 0)
371 'level' '' 'level' 369 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0
0 OPTIONAL DUMMY) (INTEGER 4 0 0 INTEGER ()) 0 0 () () 0 () () () 0 0)
370 's_field' '' 's_field' 369 ((VARIABLE INOUT UNKNOWN-PROC UNKNOWN
UNKNOWN 0 0 DUMMY) (DERIVED 256 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
465 'debug_level' '' 'debug_level' 463 ((VARIABLE IN UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0 DUMMY) (INTEGER 4 0 0 INTEGER ()) 0 0 () () 0 () ()
() 0 0)
512 'halo' '' 'halo' 511 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DUMMY) (DERIVED 144 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
513 'unns' '' 'unns' 511 ((VARIABLE OUT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DIMENSION DUMMY) (INTEGER 4 0 0 INTEGER ()) 0 0 () (1 0 EXPLICIT (
CONSTANT (INTEGER 4 0 0 INTEGER ()) 0 '1') (FUNCTION (INTEGER 4 0 0
INTEGER ()) 0 790 (('' (VARIABLE (DERIVED 144 0 0 DERIVED ()) 0 512 ())))
'node_count_halo' 1 34)) 0 () () () 0 0)
620 'halo' '' 'halo' 619 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DUMMY) (DERIVED 144 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
625 'halo' '' 'halo' 624 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DUMMY) (DERIVED 144 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
622 'unn' '' 'unn' 619 ((VARIABLE UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN
UNKNOWN 0 0 RESULT) (INTEGER 4 0 0 INTEGER ()) 0 0 () () 0 () () () 0 0)
628 'unns' '' 'unns' 624 ((VARIABLE UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN
UNKNOWN 0 0 DIMENSION RESULT ALWAYS_EXPLICIT) (INTEGER 4 0 0 INTEGER ())
0 0 () (1 0 EXPLICIT (CONSTANT (INTEGER 4 0 0 INTEGER ()) 0 '1') (
FUNCTION (INTEGER 4 0 0 INTEGER ()) 0 627 (('' (VARIABLE (INTEGER 4 0 0
INTEGER ()) 1 626 ((ARRAY (FULL 0))))) ('' ()) ('' ())) '' 0 'size')) 0
() () () 0 0)
627 'size' '(intrinsic)' 'size' 624 ((PROCEDURE UNKNOWN-INTENT
UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 FUNCTION) (REAL 8 0 0 REAL ()) 0 0 () ()
627 () () () 0 0)
626 'global_numbers' '' 'global_numbers' 624 ((VARIABLE IN UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0 DIMENSION DUMMY) (INTEGER 4 0 0 INTEGER ()) 0 0 () (
1 0 ASSUMED_SHAPE (CONSTANT (INTEGER 4 0 0 INTEGER ()) 0 '1') ()) 0 () ()
() 0 0)
759 'halo' '' 'halo' 758 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DUMMY) (DERIVED 144 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
752 'stat' '' 'stat' 748 ((VARIABLE OUT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
OPTIONAL DUMMY) (INTEGER 4 0 0 INTEGER ()) 0 0 () () 0 () () () 0 0)
760 'unn_count' '' 'unn_count' 758 ((VARIABLE UNKNOWN-INTENT
UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 RESULT) (INTEGER 4 0 0 INTEGER ()) 0 0
() () 0 () () () 0 0)
751 'universal_number' '' 'universal_number' 748 ((VARIABLE IN
UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 DUMMY) (INTEGER 4 0 0 INTEGER ()) 0 0 ()
() 0 () () () 0 0)
764 'valid' '' 'valid' 762 ((VARIABLE UNKNOWN-INTENT UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0 RESULT) (LOGICAL 4 0 0 LOGICAL ()) 0 0 () () 0 () ()
() 0 0)
373 'halo' '' 'halo' 372 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DUMMY) (DERIVED 144 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
374 'global_number' '' 'global_number' 372 ((VARIABLE IN UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0 DIMENSION DUMMY) (INTEGER 4 0 0 INTEGER ()) 0 0 () (
1 0 ASSUMED_SHAPE (CONSTANT (INTEGER 4 0 0 INTEGER ()) 0 '1') ()) 0 () ()
() 0 0)
376 'unn' '' 'unn' 372 ((VARIABLE UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN
UNKNOWN 0 0 DIMENSION RESULT ALWAYS_EXPLICIT) (INTEGER 4 0 0 INTEGER ())
0 0 () (1 0 EXPLICIT (CONSTANT (INTEGER 4 0 0 INTEGER ()) 0 '1') (
FUNCTION (INTEGER 4 0 0 INTEGER ()) 0 375 (('' (VARIABLE (INTEGER 4 0 0
INTEGER ()) 1 374 ((ARRAY (FULL 0))))) ('' ()) ('' ())) '' 0 'size')) 0
() () () 0 0)
375 'size' '(intrinsic)' 'size' 372 ((PROCEDURE UNKNOWN-INTENT
UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 FUNCTION) (REAL 8 0 0 REAL ()) 0 0 () ()
375 () () () 0 0)
763 'halo' '' 'halo' 762 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DUMMY) (DERIVED 144 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
379 'unns' '' 'unns' 377 ((VARIABLE OUT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DIMENSION DUMMY) (INTEGER 4 0 0 INTEGER ()) 0 0 () (2 0 ASSUMED_SHAPE (
CONSTANT (INTEGER 4 0 0 INTEGER ()) 0 '1') () (CONSTANT (INTEGER 4 0 0
INTEGER ()) 0 '1') ()) 0 () () () 0 0)
378 'halo' '' 'halo' 377 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DUMMY) (DERIVED 144 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
621 'global_number' '' 'global_number' 619 ((VARIABLE IN UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0 DUMMY) (INTEGER 4 0 0 INTEGER ()) 0 0 () () 0 () ()
() 0 0)
424 'halo' '' 'halo' 423 ((VARIABLE INOUT UNKNOWN-PROC UNKNOWN UNKNOWN 0
0 DUMMY) (DERIVED 144 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
425 'local_only' '' 'local_only' 423 ((VARIABLE IN UNKNOWN-PROC UNKNOWN
UNKNOWN 0 0 OPTIONAL DUMMY) (LOGICAL 4 0 0 LOGICAL ()) 0 0 () () 0 () ()
() 0 0)
464 'halo' '' 'halo' 463 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DUMMY) (DERIVED 144 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
517 'gnns' '' 'gnns' 515 ((VARIABLE OUT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DUMMY) (DERIVED 288 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
516 'halo' '' 'halo' 515 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DUMMY) (DERIVED 144 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
637 'has_gnn_to_unn' '' 'has_gnn_to_unn' 635 ((VARIABLE UNKNOWN-INTENT
UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 RESULT) (LOGICAL 4 0 0 LOGICAL ()) 0 0
() () 0 () () () 0 0)
636 'halo' '' 'halo' 635 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DUMMY) (DERIVED 144 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
750 'node' '' 'node' 748 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DUMMY) (INTEGER 4 0 0 INTEGER ()) 0 0 () () 0 () () () 0 0)
749 'halo' '' 'halo' 748 ((VARIABLE INOUT UNKNOWN-PROC UNKNOWN UNKNOWN 0
0 DUMMY) (DERIVED 144 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
428 'halo' '' 'halo' 427 ((VARIABLE INOUT UNKNOWN-PROC UNKNOWN UNKNOWN 0
0 DUMMY) (DERIVED 144 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
506 'owners' '' 'owners' 504 ((VARIABLE OUT UNKNOWN-PROC UNKNOWN UNKNOWN
0 0 DIMENSION DUMMY) (INTEGER 4 0 0 INTEGER ()) 0 0 () (1 0
ASSUMED_SHAPE (CONSTANT (INTEGER 4 0 0 INTEGER ()) 0 '1') ()) 0 () () ()
0 0)
505 'halo' '' 'halo' 504 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DUMMY) (DERIVED 144 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
540 'node_owner' '' 'node_owner' 536 ((VARIABLE UNKNOWN-INTENT
UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 RESULT ALWAYS_EXPLICIT) (INTEGER 4 0 0
INTEGER ()) 0 0 () () 0 () () () 0 0)
539 'permit_extended' '' 'permit_extended' 536 ((VARIABLE IN
UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 OPTIONAL DUMMY) (LOGICAL 4 0 0 LOGICAL
()) 0 0 () () 0 () () () 0 0)
538 'node' '' 'node' 536 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DUMMY) (INTEGER 4 0 0 INTEGER ()) 0 0 () () 0 () () () 0 0)
545 'permit_extended' '' 'permit_extended' 542 ((VARIABLE IN
UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 OPTIONAL DUMMY) (LOGICAL 4 0 0 LOGICAL
()) 0 0 () () 0 () () () 0 0)
544 'nodes' '' 'nodes' 542 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0
0 DIMENSION DUMMY) (INTEGER 4 0 0 INTEGER ()) 0 0 () (1 0 ASSUMED_SHAPE
(CONSTANT (INTEGER 4 0 0 INTEGER ()) 0 '1') ()) 0 () () () 0 0)
547 'node_owners' '' 'node_owners' 542 ((VARIABLE UNKNOWN-INTENT
UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 DIMENSION RESULT ALWAYS_EXPLICIT) (
INTEGER 4 0 0 INTEGER ()) 0 0 () (1 0 EXPLICIT (CONSTANT (INTEGER 4 0 0
INTEGER ()) 0 '1') (FUNCTION (INTEGER 4 0 0 INTEGER ()) 0 546 (('' (
VARIABLE (INTEGER 4 0 0 INTEGER ()) 1 544 ((ARRAY (FULL 0))))) ('' ()) (
'' ())) '' 0 'size')) 0 () () () 0 0)
546 'size' '(intrinsic)' 'size' 542 ((PROCEDURE UNKNOWN-INTENT
UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 FUNCTION) (REAL 8 0 0 REAL ()) 0 0 () ()
546 () () () 0 0)
616 'unns' '' 'unns' 614 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DIMENSION DUMMY) (INTEGER 4 0 0 INTEGER ()) 0 0 () (1 0 ASSUMED_SHAPE (
CONSTANT (INTEGER 4 0 0 INTEGER ()) 0 '1') ()) 0 () () () 0 0)
618 'node_owners' '' 'node_owners' 614 ((VARIABLE UNKNOWN-INTENT
UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 DIMENSION RESULT ALWAYS_EXPLICIT) (
INTEGER 4 0 0 INTEGER ()) 0 0 () (1 0 EXPLICIT (CONSTANT (INTEGER 4 0 0
INTEGER ()) 0 '1') (FUNCTION (INTEGER 4 0 0 INTEGER ()) 0 617 (('' (
VARIABLE (INTEGER 4 0 0 INTEGER ()) 1 616 ((ARRAY (FULL 0))))) ('' ()) (
'' ())) '' 0 'size')) 0 () () () 0 0)
617 'size' '(intrinsic)' 'size' 614 ((PROCEDURE UNKNOWN-INTENT
UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 FUNCTION) (REAL 8 0 0 REAL ()) 0 0 () ()
617 () () () 0 0)
615 'halo' '' 'halo' 614 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DUMMY) (DERIVED 144 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
384 'halo' '' 'halo' 383 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DUMMY) (DERIVED 144 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
385 'node' '' 'node' 383 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DUMMY) (INTEGER 4 0 0 INTEGER ()) 0 0 () () 0 () () () 0 0)
388 'nodes' '' 'nodes' 386 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0
0 DIMENSION DUMMY) (INTEGER 4 0 0 INTEGER ()) 0 0 () (1 0 ASSUMED_SHAPE
(CONSTANT (INTEGER 4 0 0 INTEGER ()) 0 '1') ()) 0 () () () 0 0)
387 'halo' '' 'halo' 386 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DUMMY) (DERIVED 144 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
543 'halo' '' 'halo' 542 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DUMMY) (DERIVED 144 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
390 'owned' '' 'owned' 386 ((VARIABLE UNKNOWN-INTENT UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0 DIMENSION RESULT ALWAYS_EXPLICIT) (LOGICAL 4 0 0
LOGICAL ()) 0 0 () (1 0 EXPLICIT (CONSTANT (INTEGER 4 0 0 INTEGER ()) 0
'1') (FUNCTION (INTEGER 4 0 0 INTEGER ()) 0 389 (('' (VARIABLE (INTEGER
4 0 0 INTEGER ()) 1 388 ((ARRAY (FULL 0))))) ('' ()) ('' ())) '' 0 'size'))
0 () () () 0 0)
389 'size' '(intrinsic)' 'size' 386 ((PROCEDURE UNKNOWN-INTENT
UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 FUNCTION) (REAL 8 0 0 REAL ()) 0 0 () ()
389 () () () 0 0)
537 'halo' '' 'halo' 536 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DUMMY) (DERIVED 144 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
509 'halo' '' 'halo' 508 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DUMMY) (DERIVED 144 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
510 'owned_nodes' '' 'owned_nodes' 508 ((VARIABLE OUT UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0 DIMENSION DUMMY) (INTEGER 4 0 0 INTEGER ()) 0 0 () (
1 0 ASSUMED_SHAPE (CONSTANT (INTEGER 4 0 0 INTEGER ()) 0 '1') ()) 0 () ()
() 0 0)
643 'halo' '' 'halo' 642 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DUMMY) (DERIVED 144 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
481 'sends' '' 'sends' 479 ((VARIABLE OUT UNKNOWN-PROC UNKNOWN UNKNOWN 0
0 DIMENSION DUMMY) (INTEGER 4 0 0 INTEGER ()) 0 0 () (1 0 ASSUMED_SHAPE
(CONSTANT (INTEGER 4 0 0 INTEGER ()) 0 '1') ()) 0 () () () 0 0)
483 'receives' '' 'receives' 479 ((VARIABLE OUT UNKNOWN-PROC UNKNOWN
UNKNOWN 0 0 DIMENSION DUMMY) (INTEGER 4 0 0 INTEGER ()) 0 0 () (1 0
ASSUMED_SHAPE (CONSTANT (INTEGER 4 0 0 INTEGER ()) 0 '1') ()) 0 () () ()
0 0)
484 'receive_starts' '' 'receive_starts' 479 ((VARIABLE OUT UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0 DIMENSION DUMMY) (INTEGER 4 0 0 INTEGER ()) 0 0 () (
1 0 ASSUMED_SHAPE (CONSTANT (INTEGER 4 0 0 INTEGER ()) 0 '1') ()) 0 () ()
() 0 0)
482 'send_starts' '' 'send_starts' 479 ((VARIABLE OUT UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0 DIMENSION DUMMY) (INTEGER 4 0 0 INTEGER ()) 0 0 () (
1 0 ASSUMED_SHAPE (CONSTANT (INTEGER 4 0 0 INTEGER ()) 0 '1') ()) 0 () ()
() 0 0)
485 'nowned_nodes' '' 'nowned_nodes' 479 ((VARIABLE OUT UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0 OPTIONAL DUMMY) (INTEGER 4 0 0 INTEGER ()) 0 0 () ()
0 () () () 0 0)
488 'halo' '' 'halo' 487 ((VARIABLE INOUT UNKNOWN-PROC UNKNOWN UNKNOWN 0
0 DUMMY) (DERIVED 144 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
491 'send_starts' '' 'send_starts' 487 ((VARIABLE IN UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0 DIMENSION DUMMY) (INTEGER 4 0 0 INTEGER ()) 0 0 () (
1 0 ASSUMED_SHAPE (CONSTANT (INTEGER 4 0 0 INTEGER ()) 0 '1') ()) 0 () ()
() 0 0)
490 'sends' '' 'sends' 487 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0
0 DIMENSION DUMMY) (INTEGER 4 0 0 INTEGER ()) 0 0 () (1 0 ASSUMED_SHAPE
(CONSTANT (INTEGER 4 0 0 INTEGER ()) 0 '1') ()) 0 () () () 0 0)
489 'nprocs' '' 'nprocs' 487 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN
0 0 DUMMY) (INTEGER 4 0 0 INTEGER ()) 0 0 () () 0 () () () 0 0)
480 'halo' '' 'halo' 479 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DUMMY) (DERIVED 144 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
494 'nowned_nodes' '' 'nowned_nodes' 487 ((VARIABLE IN UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0 OPTIONAL DUMMY) (INTEGER 4 0 0 INTEGER ()) 0 0 () ()
0 () () () 0 0)
495 'ordering_scheme' '' 'ordering_scheme' 487 ((VARIABLE IN
UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 OPTIONAL DUMMY) (INTEGER 4 0 0 INTEGER
()) 0 0 () () 0 () () () 0 0)
496 'create_caches' '' 'create_caches' 487 ((VARIABLE IN UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0 OPTIONAL DUMMY) (LOGICAL 4 0 0 LOGICAL ()) 0 0 () ()
0 () () () 0 0)
493 'receive_starts' '' 'receive_starts' 487 ((VARIABLE IN UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0 DIMENSION DUMMY) (INTEGER 4 0 0 INTEGER ()) 0 0 () (
1 0 ASSUMED_SHAPE (CONSTANT (INTEGER 4 0 0 INTEGER ()) 0 '1') ()) 0 () ()
() 0 0)
502 'inverse_node_list' '' 'inverse_node_list' 498 ((VARIABLE
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 DIMENSION DUMMY) (
INTEGER 4 0 0 INTEGER ()) 0 0 () (1 0 ASSUMED_SHAPE (CONSTANT (INTEGER 4
0 0 INTEGER ()) 0 '1') ()) 0 () () () 0 0)
783 'mesh' '' 'mesh' 781 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DUMMY) (DERIVED 118 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
782 'filename' '' 'filename' 781 ((VARIABLE IN UNKNOWN-PROC UNKNOWN
UNKNOWN 0 0 DUMMY) (CHARACTER 1 0 0 CHARACTER (())) 0 0 () () 0 () () ()
0 0)
492 'receives' '' 'receives' 487 ((VARIABLE IN UNKNOWN-PROC UNKNOWN
UNKNOWN 0 0 DIMENSION DUMMY) (INTEGER 4 0 0 INTEGER ()) 0 0 () (1 0
ASSUMED_SHAPE (CONSTANT (INTEGER 4 0 0 INTEGER ()) 0 '1') ()) 0 () () ()
0 0)
399 'mesh' '' 'mesh' 397 ((VARIABLE INOUT UNKNOWN-PROC UNKNOWN UNKNOWN 0
0 DUMMY) (DERIVED 118 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
402 'filename' '' 'filename' 401 ((VARIABLE IN UNKNOWN-PROC UNKNOWN
UNKNOWN 0 0 DUMMY) (CHARACTER 1 0 0 CHARACTER (())) 0 0 () () 0 () () ()
0 0)
400 'communicator' '' 'communicator' 397 ((VARIABLE IN UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0 OPTIONAL DUMMY) (INTEGER 4 0 0 INTEGER ()) 0 0 () ()
0 () () () 0 0)
404 'communicator' '' 'communicator' 401 ((VARIABLE IN UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0 OPTIONAL DUMMY) (INTEGER 4 0 0 INTEGER ()) 0 0 () ()
0 () () () 0 0)
403 'positions' '' 'positions' 401 ((VARIABLE INOUT UNKNOWN-PROC UNKNOWN
UNKNOWN 0 0 DUMMY) (DERIVED 108 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
398 'filename' '' 'filename' 397 ((VARIABLE IN UNKNOWN-PROC UNKNOWN
UNKNOWN 0 0 DUMMY) (CHARACTER 1 0 0 CHARACTER (())) 0 0 () () 0 () () ()
0 0)
500 'subdomain_mesh' '' 'subdomain_mesh' 498 ((VARIABLE INOUT
UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 DUMMY) (DERIVED 118 0 0 DERIVED ()) 0 0
() () 0 () () () 0 0)
499 'external_mesh' '' 'external_mesh' 498 ((VARIABLE IN UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0 DUMMY) (DERIVED 118 0 0 DERIVED ()) 0 0 () () 0 () ()
() 0 0)
501 'node_list' '' 'node_list' 498 ((VARIABLE UNKNOWN-INTENT
UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 DIMENSION DUMMY) (INTEGER 4 0 0 INTEGER
()) 0 0 () (1 0 ASSUMED_SHAPE (CONSTANT (INTEGER 4 0 0 INTEGER ()) 0 '1')
()) 0 () () () 0 0)
779 'positions' '' 'positions' 778 ((VARIABLE INOUT UNKNOWN-PROC UNKNOWN
UNKNOWN 0 0 DUMMY) (DERIVED 108 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
687 'mesh' '' 'mesh' 684 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DUMMY) (DERIVED 118 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
690 'halo' '' 'halo' 689 ((VARIABLE INOUT UNKNOWN-PROC UNKNOWN UNKNOWN 0
0 DUMMY) (DERIVED 144 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
695 'l2_halo' '' 'l2_halo' 693 ((VARIABLE IN UNKNOWN-PROC UNKNOWN
UNKNOWN 0 0 DUMMY) (DERIVED 144 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
694 'l1_halo' '' 'l1_halo' 693 ((VARIABLE INOUT UNKNOWN-PROC UNKNOWN
UNKNOWN 0 0 DUMMY) (DERIVED 144 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
691 'repair_field' '' 'repair_field' 689 ((VARIABLE IN UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0 DUMMY) (DERIVED 108 0 0 DERIVED ()) 0 0 () () 0 () ()
() 0 0)
413 'halo' '' 'halo' 412 ((VARIABLE INOUT UNKNOWN-PROC UNKNOWN UNKNOWN 0
0 DUMMY) (DERIVED 144 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
417 'repair_halo' '' 'repair_halo' 415 ((VARIABLE IN UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0 DUMMY) (DERIVED 144 0 0 DERIVED ()) 0 0 () () 0 () ()
() 0 0)
416 'halo' '' 'halo' 415 ((VARIABLE INOUT UNKNOWN-PROC UNKNOWN UNKNOWN 0
0 DUMMY) (DERIVED 144 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
414 'repair_field' '' 'repair_field' 412 ((VARIABLE IN UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0 DUMMY) (DERIVED 108 0 0 DERIVED ()) 0 0 () () 0 () ()
() 0 0)
696 'sorted_l1_halo' '' 'sorted_l1_halo' 693 ((VARIABLE IN UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0 OPTIONAL DUMMY) (LOGICAL 4 0 0 LOGICAL ()) 0 0 () ()
0 () () () 0 0)
680 'element_halo' '' 'element_halo' 679 ((VARIABLE INOUT UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0 DUMMY) (DERIVED 144 0 0 DERIVED ()) 0 0 () () 0 () ()
() 0 0)
681 'node_halo' '' 'node_halo' 679 ((VARIABLE IN UNKNOWN-PROC UNKNOWN
UNKNOWN 0 0 DUMMY) (DERIVED 144 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
682 'mesh' '' 'mesh' 679 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DUMMY) (DERIVED 118 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
685 'node_halo' '' 'node_halo' 684 ((VARIABLE INOUT UNKNOWN-PROC UNKNOWN
UNKNOWN 0 0 DUMMY) (DERIVED 144 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
686 'element_halo' '' 'element_halo' 684 ((VARIABLE IN UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0 DUMMY) (DERIVED 144 0 0 DERIVED ()) 0 0 () () 0 () ()
() 0 0)
790 'node_count' 'halos_base' 'node_count' 1 ((PROCEDURE UNKNOWN-INTENT
UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 GENERIC) (UNKNOWN 0 0 0 UNKNOWN ()) 0 0
() () 0 () () () 0 0)
)

('create_global_to_universal_numbering' 0 422 'create_ownership' 0 426
'deallocate_ownership_cache' 0 429 'deallocate_universal_numbering_cache'
0 432 'derive_element_halo_from_node_halo' 0 435
'derive_maximal_surface_element_halo' 0 440
'derive_nonperiodic_halos_from_periodic_halos' 0 447 'derive_sub_halo' 0
452 'ele_owner' 0 457 'ewrite_universal_numbers' 0 462
'extract_all_halo_receives' 0 466 'extract_all_halo_sends' 0 472
'extract_raw_halo_data' 0 478 'form_halo_from_raw_data' 0 486
'generate_substate_halos' 0 497 'get_node_owners' 0 503 'get_owned_nodes'
0 507 'get_universal_numbering' 0 7 'get_universal_numbering_inverse' 0
514 'halo_all_receives_count' 0 518 'halo_all_sends_count' 0 521
'halo_all_unique_receives_count' 0 524 'halo_data_type' 0 527 'halo_name'
0 530 'halo_node_owner' 0 535 'halo_node_owners' 0 541 'halo_nowned_nodes'
0 548 'halo_order_general' 0 551 'halo_order_trailing_receives' 0 552
'halo_ordering_scheme' 0 553 'halo_pointer' 0 556 'halo_proc_count' 0
391 'halo_receive' 0 558 'halo_receive_count' 0 563 'halo_receive_counts'
0 567 'halo_receives' 0 571 'halo_send' 0 575 'halo_send_count' 0 580
'halo_send_counts' 0 584 'halo_sends' 0 588 'halo_type' 0 144
'halo_type_cg_node' 0 606 'halo_type_dg_node' 0 607 'halo_type_element'
0 608 'halo_unique_receive_count' 0 609 'halo_universal_node_owners' 0
613 'halo_universal_number' 0 13 'halo_universal_numbers' 0 623
'halo_valid_for_communication' 0 629 'halos' 0 633
'has_global_to_universal_numbering' 0 634 'has_nowned_nodes' 0 638
'has_ownership' 0 641 'invert_comms_sizes' 0 644 'max_halo_node' 0 650
'max_halo_receive_node' 0 654 'max_halo_send_node' 0 658 'min_halo_node'
0 662 'min_halo_receive_node' 0 666 'min_halo_send_node' 0 670
'print_halo' 0 674 'reorder_element_halo' 0 678
'reorder_halo_from_element_halo' 0 683 'reorder_halo_receives' 0 688
'reorder_l1_from_l2_halo' 0 692 'set_all_halo_receives' 0 697
'set_all_halo_sends' 0 701 'set_halo_communicator' 0 705
'set_halo_data_type' 0 709 'set_halo_name' 0 713 'set_halo_nowned_nodes'
0 717 'set_halo_ordering_scheme' 0 721 'set_halo_receive' 0 725
'set_halo_receives' 0 731 'set_halo_send' 0 736 'set_halo_sends' 0 742
'set_halo_universal_number' 0 747 'trailing_receives_consistent' 0 753
'universal_numbering_count' 0 757 'valid_global_to_universal_numbering'
0 761 'valid_halo_communicator' 0 765 'valid_halo_node_counts' 0 769
'valid_serial_halo' 0 773 'verify_halos' 0 777 'write_halos' 0 780
'write_universal_numbering' 0 784)
