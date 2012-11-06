GFORTRAN module version '6' created from Adjoint_Functionals.F90 on Tue Oct  2 17:41:57 2012
MD5:c8e959aa679e027baf58531fe12dbcb0 -- If you edit this, you'll get what you deserve.

(() () () () () () () () () () () () () () () () () () () () () () ()
() () () ())

()

()

(('mpi_fortran_argv_null' 2 0 0 'mpi_fortran_argv_null') (
'mpi_fortran_bottom' 3 0 0 'mpi_fortran_bottom') ('petscfortran4' 4 0 0
'petscfortran4') ('petscfortran5' 5 0 0 'petscfortran5') ('petscfortran6'
6 0 0 'petscfortran6') ('petscfortran7' 7 0 0 'petscfortran7') (
'petscfortran8' 8 0 0 'petscfortran8') ('petscfortran9' 9 0 0
'petscfortran9') ('mpi_fortran_argvs_null' 10 0 0 'mpi_fortran_argvs_null')
('mpi_fortran_errcodes_ignore' 11 0 0 'mpi_fortran_errcodes_ignore') (
'mpi_fortran_in_place' 12 0 0 'mpi_fortran_in_place') (
'mpi_fortran_status_ignore' 13 0 0 'mpi_fortran_status_ignore') (
'mpi_fortran_statuses_ignore' 14 0 0 'mpi_fortran_statuses_ignore') (
'petscfortran1' 15 0 0 'petscfortran1') ('petscfortran10' 16 0 0
'petscfortran10') ('petscfortran2' 17 0 0 'petscfortran2') (
'petscfortran3' 18 0 0 'petscfortran3'))

()

()

(19 'functional' 'adjoint_functionals' 'functional' 1 ((PROCEDURE
UNKNOWN-INTENT MODULE-PROC DECL UNKNOWN 0 0 FUNCTION ALWAYS_EXPLICIT) (
REAL 8 0 0 REAL ()) 20 0 (21 22 23 24) () 25 () () () 0 0)
26 'functional_derivative' 'adjoint_functionals' 'functional_derivative'
1 ((PROCEDURE UNKNOWN-INTENT MODULE-PROC DECL UNKNOWN 0 0 SUBROUTINE
ALWAYS_EXPLICIT) (UNKNOWN 0 0 0 UNKNOWN ()) 27 0 (28 29 30 31 32) () 0 ()
() () 0 0)
21 'states' '' 'states' 20 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0
0 DIMENSION DUMMY) (DERIVED 33 0 0 DERIVED ()) 0 0 () (2 0 ASSUMED_SHAPE
(CONSTANT (INTEGER 4 0 0 INTEGER ()) 0 '1') () (CONSTANT (INTEGER 4 0 0
INTEGER ()) 0 '1') ()) 0 () () () 0 0)
22 'current_time' '' 'current_time' 20 ((VARIABLE IN UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0 DUMMY) (REAL 8 0 0 REAL ()) 0 0 () () 0 () () () 0 0)
23 'dt' '' 'dt' 20 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 DUMMY)
(REAL 8 0 0 REAL ()) 0 0 () () 0 () () () 0 0)
24 'n' '' 'n' 20 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 DUMMY) (
INTEGER 4 0 0 INTEGER ()) 0 0 () () 0 () () () 0 0)
25 'j' '' 'j' 20 ((VARIABLE UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN
0 0 RESULT ALWAYS_EXPLICIT) (REAL 8 0 0 REAL ()) 0 0 () () 0 () () () 0
0)
28 'states' '' 'states' 27 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0
0 DIMENSION DUMMY) (DERIVED 33 0 0 DERIVED ()) 0 0 () (2 0 ASSUMED_SHAPE
(CONSTANT (INTEGER 4 0 0 INTEGER ()) 0 '1') () (CONSTANT (INTEGER 4 0 0
INTEGER ()) 0 '1') ()) 0 () () () 0 0)
29 'current_time' '' 'current_time' 27 ((VARIABLE IN UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0 DUMMY) (REAL 8 0 0 REAL ()) 0 0 () () 0 () () () 0 0)
30 'dt' '' 'dt' 27 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 DUMMY)
(REAL 8 0 0 REAL ()) 0 0 () () 0 () () () 0 0)
31 'n' '' 'n' 27 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 DUMMY) (
INTEGER 4 0 0 INTEGER ()) 0 0 () () 0 () () () 0 0)
32 'derivatives' '' 'derivatives' 27 ((VARIABLE OUT UNKNOWN-PROC UNKNOWN
UNKNOWN 0 0 DIMENSION DUMMY) (DERIVED 33 0 0 DERIVED ()) 0 0 () (1 0
EXPLICIT (CONSTANT (INTEGER 4 0 0 INTEGER ()) 0 '1') (FUNCTION (INTEGER
4 0 0 INTEGER ()) 0 34 (('' (VARIABLE (DERIVED 33 0 0 DERIVED ()) 2 28 (
(ARRAY (FULL 0))))) ('' (CONSTANT (INTEGER 4 0 0 INTEGER ()) 0 '1')) (''
())) '' 0 'size')) 0 () () () 0 0)
12 'mpi_in_place' 'mpi_interfaces' 'mpi_in_place' 1 ((VARIABLE
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 IN_COMMON) (INTEGER 4 0
0 INTEGER ()) 0 0 () () 0 () () () 0 0)
11 'mpi_errcodes_ignore' 'mpi_interfaces' 'mpi_errcodes_ignore' 1 ((
VARIABLE UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 DIMENSION
IN_COMMON) (INTEGER 4 0 0 INTEGER ()) 0 0 () (1 0 EXPLICIT (CONSTANT (
INTEGER 4 0 0 INTEGER ()) 0 '1') (CONSTANT (INTEGER 4 0 0 INTEGER ()) 0
'1')) 0 () () () 0 0)
2 'mpi_argv_null' 'mpi_interfaces' 'mpi_argv_null' 1 ((VARIABLE
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 DIMENSION IN_COMMON) (
CHARACTER 1 0 0 CHARACTER ((CONSTANT (INTEGER 4 0 0 INTEGER ()) 0 '1')))
0 0 () (1 0 EXPLICIT (CONSTANT (INTEGER 4 0 0 INTEGER ()) 0 '1') (
CONSTANT (INTEGER 4 0 0 INTEGER ()) 0 '1')) 0 () () () 0 0)
15 'petsc_null_character' 'petscsysdef' 'petsc_null_character' 1 ((
VARIABLE UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN IMPLICIT-SAVE 0 0 IN_COMMON)
(CHARACTER 1 0 0 CHARACTER ((CONSTANT (INTEGER 4 0 0 INTEGER ()) 0 '80')))
0 0 () () 0 () () () 0 0)
18 'petsc_null' 'petscsysdef' 'petsc_null' 1 ((VARIABLE UNKNOWN-INTENT
UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 IN_COMMON) (INTEGER 4 0 0 INTEGER ()) 0
0 () () 0 () () () 0 0)
17 'petsc_null_integer' 'petscsysdef' 'petsc_null_integer' 1 ((VARIABLE
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 IN_COMMON) (INTEGER 4 0
0 INTEGER ()) 0 0 () () 0 () () () 0 0)
10 'mpi_argvs_null' 'mpi_interfaces' 'mpi_argvs_null' 1 ((VARIABLE
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 IN_COMMON) (REAL 8 0 0
REAL ()) 0 0 () () 0 () () () 0 0)
16 'petsc_comm_self' 'petscsysdef' 'petsc_comm_self' 1 ((VARIABLE
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 IN_COMMON) (INTEGER 4 0
0 INTEGER ()) 0 0 () () 0 () () () 0 0)
14 'mpi_statuses_ignore' 'mpi_interfaces' 'mpi_statuses_ignore' 1 ((
VARIABLE UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 IN_COMMON) (
REAL 8 0 0 REAL ()) 0 0 () () 0 () () () 0 0)
13 'mpi_status_ignore' 'mpi_interfaces' 'mpi_status_ignore' 1 ((
VARIABLE UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 DIMENSION
IN_COMMON) (INTEGER 4 0 0 INTEGER ()) 0 0 () (1 0 EXPLICIT (CONSTANT (
INTEGER 4 0 0 INTEGER ()) 0 '1') (CONSTANT (INTEGER 4 0 0 INTEGER ()) 0
'5')) 0 () () () 0 0)
5 'petsc_null_double' 'petscsysdef' 'petsc_null_double' 1 ((VARIABLE
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 IN_COMMON) (REAL 8 0 0
REAL ()) 0 0 () () 0 () () () 0 0)
6 'petsc_null_real' 'petscsysdef' 'petsc_null_real' 1 ((VARIABLE
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 IN_COMMON) (REAL 8 0 0
REAL ()) 0 0 () () 0 () () () 0 0)
4 'petsc_null_scalar' 'petscsysdef' 'petsc_null_scalar' 1 ((VARIABLE
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 IN_COMMON) (REAL 8 0 0
REAL ()) 0 0 () () 0 () () () 0 0)
3 'mpi_bottom' 'mpi_interfaces' 'mpi_bottom' 1 ((VARIABLE UNKNOWN-INTENT
UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 IN_COMMON) (INTEGER 4 0 0 INTEGER ()) 0
0 () () 0 () () () 0 0)
9 'petsc_comm_world' 'petscsysdef' 'petsc_comm_world' 1 ((VARIABLE
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 IN_COMMON) (INTEGER 4 0
0 INTEGER ()) 0 0 () () 0 () () () 0 0)
8 'petsc_null_object' 'petscsysdef' 'petsc_null_object' 1 ((VARIABLE
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 IN_COMMON) (INTEGER 8 0
0 INTEGER ()) 0 0 () () 0 () () () 0 0)
7 'petsc_null_truth' 'petscsysdef' 'petsc_null_truth' 1 ((VARIABLE
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 IN_COMMON) (LOGICAL 4 0
0 LOGICAL ()) 0 0 () () 0 () () () 0 0)
33 'state_type' 'state_module' 'state_type' 1 ((DERIVED UNKNOWN-INTENT
UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 POINTER_COMP) (UNKNOWN 0 0 0 UNKNOWN ())
0 0 () () 0 ((35 'name' (CHARACTER 1 0 0 CHARACTER ((CONSTANT (INTEGER 4
0 0 INTEGER ()) 0 '101'))) () (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0) UNKNOWN-ACCESS (CONSTANT (CHARACTER 1 0 0 CHARACTER
(())) 0 101
'                                                                                                     '))
(36 'option_path' (CHARACTER 1 0 0 CHARACTER ((CONSTANT (INTEGER 4 0 0
INTEGER ()) 0 '8192'))) () (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0) UNKNOWN-ACCESS (CONSTANT (CHARACTER 1 0 0 CHARACTER
(())) 0 8192
'/uninitialised_path/                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            '))
(37 'vector_names' (CHARACTER 1 0 0 CHARACTER ((CONSTANT (INTEGER 4 0 0
INTEGER ()) 0 '101'))) (1 0 DEFERRED () ()) (UNKNOWN-FL UNKNOWN-INTENT
UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 DIMENSION POINTER) UNKNOWN-ACCESS (
NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0)) (38 'scalar_names' (CHARACTER 1 0 0
CHARACTER ((CONSTANT (INTEGER 4 0 0 INTEGER ()) 0 '101'))) (1 0 DEFERRED
() ()) (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DIMENSION POINTER) UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0)) (
39 'mesh_names' (CHARACTER 1 0 0 CHARACTER ((CONSTANT (INTEGER 4 0 0
INTEGER ()) 0 '101'))) (1 0 DEFERRED () ()) (UNKNOWN-FL UNKNOWN-INTENT
UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 DIMENSION POINTER) UNKNOWN-ACCESS (
NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0)) (40 'halo_names' (CHARACTER 1 0 0
CHARACTER ((CONSTANT (INTEGER 4 0 0 INTEGER ()) 0 '101'))) (1 0 DEFERRED
() ()) (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DIMENSION POINTER) UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0)) (
41 'tensor_names' (CHARACTER 1 0 0 CHARACTER ((CONSTANT (INTEGER 4 0 0
INTEGER ()) 0 '101'))) (1 0 DEFERRED () ()) (UNKNOWN-FL UNKNOWN-INTENT
UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 DIMENSION POINTER) UNKNOWN-ACCESS (
NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0)) (42 'csr_sparsity_names' (CHARACTER
1 0 0 CHARACTER ((CONSTANT (INTEGER 4 0 0 INTEGER ()) 0 '101'))) (1 0
DEFERRED () ()) (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN
0 0 DIMENSION POINTER) UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0))
(43 'csr_matrix_names' (CHARACTER 1 0 0 CHARACTER ((CONSTANT (INTEGER 4
0 0 INTEGER ()) 0 '101'))) (1 0 DEFERRED () ()) (UNKNOWN-FL
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 DIMENSION POINTER)
UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0)) (44
'block_csr_matrix_names' (CHARACTER 1 0 0 CHARACTER ((CONSTANT (INTEGER
4 0 0 INTEGER ()) 0 '101'))) (1 0 DEFERRED () ()) (UNKNOWN-FL
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 DIMENSION POINTER)
UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0)) (45
'petsc_csr_matrix_names' (CHARACTER 1 0 0 CHARACTER ((CONSTANT (INTEGER
4 0 0 INTEGER ()) 0 '101'))) (1 0 DEFERRED () ()) (UNKNOWN-FL
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 DIMENSION POINTER)
UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0)) (46 'vector_fields'
(DERIVED 47 0 0 DERIVED ()) (1 0 DEFERRED () ()) (UNKNOWN-FL
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 DIMENSION POINTER)
UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0)) (48 'tensor_fields'
(DERIVED 49 0 0 DERIVED ()) (1 0 DEFERRED () ()) (UNKNOWN-FL
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 DIMENSION POINTER)
UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0)) (50 'scalar_fields'
(DERIVED 51 0 0 DERIVED ()) (1 0 DEFERRED () ()) (UNKNOWN-FL
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 DIMENSION POINTER)
UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0)) (52 'meshes' (
DERIVED 53 0 0 DERIVED ()) (1 0 DEFERRED () ()) (UNKNOWN-FL
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 DIMENSION POINTER)
UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0)) (54 'halos' (
DERIVED 55 0 0 DERIVED ()) (1 0 DEFERRED () ()) (UNKNOWN-FL
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 DIMENSION POINTER)
UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0)) (56 'csr_sparsities'
(DERIVED 57 0 0 DERIVED ()) (1 0 DEFERRED () ()) (UNKNOWN-FL
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 DIMENSION POINTER)
UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0)) (58 'csr_matrices' (
DERIVED 59 0 0 DERIVED ()) (1 0 DEFERRED () ()) (UNKNOWN-FL
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 DIMENSION POINTER)
UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0)) (60
'block_csr_matrices' (DERIVED 61 0 0 DERIVED ()) (1 0 DEFERRED () ()) (
UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 DIMENSION
POINTER) UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0)) (62
'petsc_csr_matrices' (DERIVED 63 0 0 DERIVED ()) (1 0 DEFERRED () ()) (
UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 DIMENSION
POINTER) UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0))) PUBLIC (()
() () ()) () 0 0 7575341)
61 'block_csr_matrix_pointer' 'sparse_tools' 'block_csr_matrix_pointer'
1 ((DERIVED UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 POINTER_COMP)
(UNKNOWN 0 0 0 UNKNOWN ()) 0 0 () () 0 ((64 'ptr' (DERIVED 65 0 0
DERIVED ()) () (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0
0 POINTER) UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0))) PUBLIC (
() () () ()) () 0 0 83606449)
55 'halo_pointer' 'halo_data_types' 'halo_pointer' 1 ((DERIVED
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 POINTER_COMP) (UNKNOWN 0
0 0 UNKNOWN ()) 0 0 () () 0 ((66 'ptr' (DERIVED 67 0 0 DERIVED ()) () (
UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 POINTER)
UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0))) PUBLIC (() () () ())
() 0 0 72023282)
65 'block_csr_matrix' 'sparse_tools' 'block_csr_matrix' 1 ((DERIVED
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 POINTER_COMP) (UNKNOWN 0
0 0 UNKNOWN ()) 0 0 () () 0 ((68 'sparsity' (DERIVED 69 0 0 DERIVED ())
() (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0)
UNKNOWN-ACCESS (STRUCTURE (DERIVED 69 0 0 DERIVED ()) 0 (((NULL (
UNKNOWN 0 0 0 UNKNOWN ()) 0) ()) ((NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0) ())
((NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0) ()) (() ()) ((NULL (UNKNOWN 0 0 0
UNKNOWN ()) 0) ()) ((NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0) ()) ((NULL (
UNKNOWN 0 0 0 UNKNOWN ()) 0) ()) ((CONSTANT (CHARACTER 1 0 0 CHARACTER (
())) 0 101
'                                                                                                     ')
()) ((CONSTANT (LOGICAL 4 0 0 LOGICAL ()) 0 0) ()) ((CONSTANT (LOGICAL 4
0 0 LOGICAL ()) 0 0) ())) ())) (70 'val' (DERIVED 71 0 0 DERIVED ()) (2
0 DEFERRED () () () ()) (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN
UNKNOWN 0 0 DIMENSION POINTER) UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0
UNKNOWN ()) 0)) (72 'contiguous_val' (REAL 8 0 0 REAL ()) (1 0 DEFERRED
() ()) (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DIMENSION POINTER) UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0)) (
73 'ival' (DERIVED 74 0 0 DERIVED ()) (2 0 DEFERRED () () () ()) (
UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 DIMENSION
POINTER) UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0)) (75 'blocks'
(INTEGER 4 0 0 INTEGER ()) (1 0 EXPLICIT (CONSTANT (INTEGER 4 0 0
INTEGER ()) 0 '1') (CONSTANT (INTEGER 4 0 0 INTEGER ()) 0 '2')) (
UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 DIMENSION)
UNKNOWN-ACCESS (ARRAY (INTEGER 4 0 0 INTEGER ()) 1 (((CONSTANT (INTEGER
4 0 0 INTEGER ()) 0 '0') ()) ((CONSTANT (INTEGER 4 0 0 INTEGER ()) 0 '0')
())) ('2'))) (76 'clone' (LOGICAL 4 0 0 LOGICAL ()) () (UNKNOWN-FL
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0) UNKNOWN-ACCESS (
CONSTANT (LOGICAL 4 0 0 LOGICAL ()) 0 0)) (77 'external_val' (LOGICAL 4
0 0 LOGICAL ()) () (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN
UNKNOWN 0 0) UNKNOWN-ACCESS (CONSTANT (LOGICAL 4 0 0 LOGICAL ()) 0 0)) (
78 'columns' (INTEGER 4 0 0 INTEGER ()) () (UNKNOWN-FL UNKNOWN-INTENT
UNKNOWN-PROC UNKNOWN UNKNOWN 0 0) UNKNOWN-ACCESS ()) (79 'refcount' (
DERIVED 80 0 0 DERIVED ()) () (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0 POINTER) UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN
()) 0)) (81 'name' (CHARACTER 1 0 0 CHARACTER ((CONSTANT (INTEGER 4 0 0
INTEGER ()) 0 '101'))) () (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0) UNKNOWN-ACCESS (CONSTANT (CHARACTER 1 0 0 CHARACTER
(())) 0 101
'                                                                                                     '))
(82 'diagonal' (LOGICAL 4 0 0 LOGICAL ()) () (UNKNOWN-FL UNKNOWN-INTENT
UNKNOWN-PROC UNKNOWN UNKNOWN 0 0) UNKNOWN-ACCESS (CONSTANT (LOGICAL 4 0
0 LOGICAL ()) 0 0)) (83 'equal_diagonal_blocks' (LOGICAL 4 0 0 LOGICAL ())
() (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0)
UNKNOWN-ACCESS (CONSTANT (LOGICAL 4 0 0 LOGICAL ()) 0 0)) (84 'ksp' (
INTEGER 8 0 0 INTEGER ()) () (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0 POINTER) UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN
()) 0))) PUBLIC (() () () ()) () 0 0 61921235)
59 'csr_matrix_pointer' 'sparse_tools' 'csr_matrix_pointer' 1 ((DERIVED
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 POINTER_COMP) (UNKNOWN 0
0 0 UNKNOWN ()) 0 0 () () 0 ((85 'ptr' (DERIVED 86 0 0 DERIVED ()) () (
UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 POINTER)
UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0))) PUBLIC (() () () ())
() 0 0 82770239)
69 'csr_sparsity' 'sparse_tools' 'csr_sparsity' 1 ((DERIVED
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 POINTER_COMP) (UNKNOWN 0
0 0 UNKNOWN ()) 0 0 () () 0 ((87 'findrm' (INTEGER 4 0 0 INTEGER ()) (1
0 DEFERRED () ()) (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN
UNKNOWN 0 0 DIMENSION POINTER) UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0
UNKNOWN ()) 0)) (88 'centrm' (INTEGER 4 0 0 INTEGER ()) (1 0 DEFERRED ()
()) (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DIMENSION POINTER) UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0)) (
89 'colm' (INTEGER 4 0 0 INTEGER ()) (1 0 DEFERRED () ()) (UNKNOWN-FL
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 DIMENSION POINTER)
UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0)) (90 'columns' (
INTEGER 4 0 0 INTEGER ()) () (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0) UNKNOWN-ACCESS ()) (91 'row_halo' (DERIVED 67 0 0
DERIVED ()) () (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0
0 POINTER) UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0)) (92
'column_halo' (DERIVED 67 0 0 DERIVED ()) () (UNKNOWN-FL UNKNOWN-INTENT
UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 POINTER) UNKNOWN-ACCESS (NULL (UNKNOWN
0 0 0 UNKNOWN ()) 0)) (93 'refcount' (DERIVED 80 0 0 DERIVED ()) () (
UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 POINTER)
UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0)) (94 'name' (
CHARACTER 1 0 0 CHARACTER ((CONSTANT (INTEGER 4 0 0 INTEGER ()) 0 '101')))
() (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0)
UNKNOWN-ACCESS (CONSTANT (CHARACTER 1 0 0 CHARACTER (())) 0 101
'                                                                                                     '))
(95 'wrapped' (LOGICAL 4 0 0 LOGICAL ()) () (UNKNOWN-FL UNKNOWN-INTENT
UNKNOWN-PROC UNKNOWN UNKNOWN 0 0) UNKNOWN-ACCESS (CONSTANT (LOGICAL 4 0
0 LOGICAL ()) 0 0)) (96 'sorted_rows' (LOGICAL 4 0 0 LOGICAL ()) () (
UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0)
UNKNOWN-ACCESS (CONSTANT (LOGICAL 4 0 0 LOGICAL ()) 0 0))) PUBLIC (() ()
() ()) () 0 0 78882345)
57 'csr_sparsity_pointer' 'sparse_tools' 'csr_sparsity_pointer' 1 ((
DERIVED UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 POINTER_COMP) (
UNKNOWN 0 0 0 UNKNOWN ()) 0 0 () () 0 ((97 'ptr' (DERIVED 69 0 0 DERIVED
()) () (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
POINTER) UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0))) PUBLIC (()
() () ()) () 0 0 9687815)
74 'integer_vector' 'futils' 'integer_vector' 1 ((DERIVED UNKNOWN-INTENT
UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 POINTER_COMP) (UNKNOWN 0 0 0 UNKNOWN ())
0 0 () () 0 ((98 'ptr' (INTEGER 4 0 0 INTEGER ()) (1 0 DEFERRED () ()) (
UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 DIMENSION
POINTER) UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0))) PUBLIC (()
() () ()) () 0 0 9661976)
53 'mesh_pointer' 'fields_data_types' 'mesh_pointer' 1 ((DERIVED
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 POINTER_COMP) (UNKNOWN 0
0 0 UNKNOWN ()) 0 0 () () 0 ((99 'ptr' (DERIVED 100 0 0 DERIVED ()) () (
UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 POINTER)
UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0))) PUBLIC (() () () ())
() 0 0 65432384)
100 'mesh_type' 'fields_data_types' 'mesh_type' 1 ((DERIVED
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 POINTER_COMP) (UNKNOWN 0
0 0 UNKNOWN ()) 0 0 () () 0 ((101 'ndglno' (INTEGER 4 0 0 INTEGER ()) (
1 0 DEFERRED () ()) (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN
UNKNOWN 0 0 DIMENSION POINTER) UNKNOWN-ACCESS ()) (102 'wrapped' (
LOGICAL 4 0 0 LOGICAL ()) () (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0) UNKNOWN-ACCESS (CONSTANT (LOGICAL 4 0 0 LOGICAL ())
0 1)) (103 'shape' (DERIVED 104 0 0 DERIVED ()) () (UNKNOWN-FL
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0) UNKNOWN-ACCESS (
STRUCTURE (DERIVED 104 0 0 DERIVED ()) 0 ((() ()) (() ()) (() ()) (() ())
((NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0) ()) ((NULL (UNKNOWN 0 0 0 UNKNOWN ())
0) ()) ((NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0) ()) ((NULL (UNKNOWN 0 0 0
UNKNOWN ()) 0) ()) ((NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0) ()) ((NULL (
UNKNOWN 0 0 0 UNKNOWN ()) 0) ()) ((NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0) ())
((STRUCTURE (DERIVED 105 0 0 DERIVED ()) 0 ((() ()) (() ()) (() ()) (()
()) ((NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0) ()) ((NULL (UNKNOWN 0 0 0
UNKNOWN ()) 0) ()) (() ()) ((NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0) ()) (()
())) ()) ()) ((NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0) ()) ((NULL (UNKNOWN 0
0 0 UNKNOWN ()) 0) ()) ((NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0) ()) ((NULL (
UNKNOWN 0 0 0 UNKNOWN ()) 0) ()) (() ())) ())) (106 'elements' (INTEGER
4 0 0 INTEGER ()) () (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN
UNKNOWN 0 0) UNKNOWN-ACCESS ()) (107 'nodes' (INTEGER 4 0 0 INTEGER ())
() (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0)
UNKNOWN-ACCESS ()) (108 'name' (CHARACTER 1 0 0 CHARACTER ((CONSTANT (
INTEGER 4 0 0 INTEGER ()) 0 '101'))) () (UNKNOWN-FL UNKNOWN-INTENT
UNKNOWN-PROC UNKNOWN UNKNOWN 0 0) UNKNOWN-ACCESS ()) (109 'option_path'
(CHARACTER 1 0 0 CHARACTER ((CONSTANT (INTEGER 4 0 0 INTEGER ()) 0 '8192')))
() (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0)
UNKNOWN-ACCESS (CONSTANT (CHARACTER 1 0 0 CHARACTER (())) 0 8192
'/uninitialised_path/                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            '))
(110 'continuity' (INTEGER 4 0 0 INTEGER ()) () (UNKNOWN-FL
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0) UNKNOWN-ACCESS (
CONSTANT (INTEGER 4 0 0 INTEGER ()) 0 '0')) (111 'refcount' (DERIVED 80
0 0 DERIVED ()) () (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN
UNKNOWN 0 0 POINTER) UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0))
(112 'faces' (DERIVED 113 0 0 DERIVED ()) () (UNKNOWN-FL UNKNOWN-INTENT
UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 POINTER) UNKNOWN-ACCESS (NULL (UNKNOWN
0 0 0 UNKNOWN ()) 0)) (114 'subdomain_mesh' (DERIVED 115 0 0 DERIVED ())
() (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 POINTER)
UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0)) (116 'adj_lists' (
DERIVED 117 0 0 DERIVED ()) () (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0 POINTER) UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN
()) 0)) (118 'columns' (INTEGER 4 0 0 INTEGER ()) (1 0 DEFERRED () ()) (
UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 DIMENSION
POINTER) UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0)) (119
'element_columns' (INTEGER 4 0 0 INTEGER ()) (1 0 DEFERRED () ()) (
UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 DIMENSION
POINTER) UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0)) (120
'region_ids' (INTEGER 4 0 0 INTEGER ()) (1 0 DEFERRED () ()) (
UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 DIMENSION
POINTER) UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0)) (121 'halos'
(DERIVED 67 0 0 DERIVED ()) (1 0 DEFERRED () ()) (UNKNOWN-FL
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 DIMENSION POINTER)
UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0)) (122 'element_halos'
(DERIVED 67 0 0 DERIVED ()) (1 0 DEFERRED () ()) (UNKNOWN-FL
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 DIMENSION POINTER)
UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0)) (123 'periodic' (
LOGICAL 4 0 0 LOGICAL ()) () (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0) UNKNOWN-ACCESS (CONSTANT (LOGICAL 4 0 0 LOGICAL ())
0 0))) PUBLIC (() () () ()) () 0 0 2572791)
105 'quadrature_type' 'quadrature' 'quadrature_type' 1 ((DERIVED
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 POINTER_COMP) (UNKNOWN 0
0 0 UNKNOWN ()) 0 0 () () 0 ((124 'dim' (INTEGER 4 0 0 INTEGER ()) () (
UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0)
UNKNOWN-ACCESS ()) (125 'degree' (INTEGER 4 0 0 INTEGER ()) () (
UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0)
UNKNOWN-ACCESS ()) (126 'vertices' (INTEGER 4 0 0 INTEGER ()) () (
UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0)
UNKNOWN-ACCESS ()) (127 'ngi' (INTEGER 4 0 0 INTEGER ()) () (UNKNOWN-FL
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0) UNKNOWN-ACCESS ()) (
128 'weight' (REAL 8 0 0 REAL ()) (1 0 DEFERRED () ()) (UNKNOWN-FL
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 DIMENSION POINTER)
UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0)) (129 'l' (REAL 8 0 0
REAL ()) (2 0 DEFERRED () () () ()) (UNKNOWN-FL UNKNOWN-INTENT
UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 DIMENSION POINTER) UNKNOWN-ACCESS (
NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0)) (130 'name' (CHARACTER 1 0 0
CHARACTER ((CONSTANT (INTEGER 4 0 0 INTEGER ()) 0 '0'))) () (UNKNOWN-FL
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0) UNKNOWN-ACCESS ()) (
131 'refcount' (DERIVED 80 0 0 DERIVED ()) () (UNKNOWN-FL UNKNOWN-INTENT
UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 POINTER) UNKNOWN-ACCESS (NULL (UNKNOWN
0 0 0 UNKNOWN ()) 0)) (132 'family' (INTEGER 4 0 0 INTEGER ()) () (
UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0)
UNKNOWN-ACCESS ())) PUBLIC (() () () ()) () 0 0 59837722)
67 'halo_type' 'halo_data_types' 'halo_type' 1 ((DERIVED UNKNOWN-INTENT
UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 POINTER_COMP) (UNKNOWN 0 0 0 UNKNOWN ())
0 0 () () 0 ((133 'name' (CHARACTER 1 0 0 CHARACTER ((CONSTANT (INTEGER
4 0 0 INTEGER ()) 0 '101'))) () (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0) UNKNOWN-ACCESS ()) (134 'refcount' (DERIVED 80 0 0
DERIVED ()) () (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0
0 POINTER) UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0)) (135
'data_type' (INTEGER 4 0 0 INTEGER ()) () (UNKNOWN-FL UNKNOWN-INTENT
UNKNOWN-PROC UNKNOWN UNKNOWN 0 0) UNKNOWN-ACCESS (CONSTANT (INTEGER 4 0
0 INTEGER ()) 0 '0')) (136 'ordering_scheme' (INTEGER 4 0 0 INTEGER ())
() (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0)
UNKNOWN-ACCESS (CONSTANT (INTEGER 4 0 0 INTEGER ()) 0 '0')) (137
'communicator' (INTEGER 4 0 0 INTEGER ()) () (UNKNOWN-FL UNKNOWN-INTENT
UNKNOWN-PROC UNKNOWN UNKNOWN 0 0) UNKNOWN-ACCESS ()) (138 'nprocs' (
INTEGER 4 0 0 INTEGER ()) () (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0) UNKNOWN-ACCESS (CONSTANT (INTEGER 4 0 0 INTEGER ())
0 '0')) (139 'sends' (DERIVED 74 0 0 DERIVED ()) (1 0 DEFERRED () ()) (
UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 DIMENSION
POINTER) UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0)) (140
'receives' (DERIVED 74 0 0 DERIVED ()) (1 0 DEFERRED () ()) (UNKNOWN-FL
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 DIMENSION POINTER)
UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0)) (141 'nowned_nodes'
(INTEGER 4 0 0 INTEGER ()) () (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0) UNKNOWN-ACCESS (CONSTANT (INTEGER 4 0 0 INTEGER ())
0 '-1')) (142 'owners' (INTEGER 4 0 0 INTEGER ()) (1 0 DEFERRED () ()) (
UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 DIMENSION
POINTER) UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0)) (143
'unn_count' (INTEGER 4 0 0 INTEGER ()) () (UNKNOWN-FL UNKNOWN-INTENT
UNKNOWN-PROC UNKNOWN UNKNOWN 0 0) UNKNOWN-ACCESS (CONSTANT (INTEGER 4 0
0 INTEGER ()) 0 '-1')) (144 'owned_nodes_unn_base' (INTEGER 4 0 0
INTEGER ()) (1 0 DEFERRED () ()) (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0 DIMENSION POINTER) UNKNOWN-ACCESS (NULL (UNKNOWN 0 0
0 UNKNOWN ()) 0)) (145 'my_owned_nodes_unn_base' (INTEGER 4 0 0 INTEGER
()) () (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0)
UNKNOWN-ACCESS (CONSTANT (INTEGER 4 0 0 INTEGER ()) 0 '-1')) (146
'receives_gnn_to_unn' (INTEGER 4 0 0 INTEGER ()) (1 0 DEFERRED () ()) (
UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 DIMENSION
POINTER) UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0)) (147
'gnn_to_unn' (INTEGER 4 0 0 INTEGER ()) (1 0 DEFERRED () ()) (
UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 DIMENSION
POINTER) UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0))) PUBLIC (()
() () ()) () 0 0 63994053)
63 'petsc_csr_matrix_pointer' 'sparse_tools_petsc'
'petsc_csr_matrix_pointer' 1 ((DERIVED UNKNOWN-INTENT UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0 POINTER_COMP) (UNKNOWN 0 0 0 UNKNOWN ()) 0 0 () () 0
((148 'ptr' (DERIVED 149 0 0 DERIVED ()) () (UNKNOWN-FL UNKNOWN-INTENT
UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 POINTER) UNKNOWN-ACCESS (NULL (UNKNOWN
0 0 0 UNKNOWN ()) 0))) PUBLIC (() () () ()) () 0 0 85445535)
117 'adjacency_cache' 'fields_data_types' 'adjacency_cache' 1 ((DERIVED
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 POINTER_COMP) (UNKNOWN 0
0 0 UNKNOWN ()) 0 0 () () 0 ((150 'nnlist' (DERIVED 69 0 0 DERIVED ()) ()
(UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 POINTER)
UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0)) (151 'nelist' (
DERIVED 69 0 0 DERIVED ()) () (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0 POINTER) UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN
()) 0)) (152 'eelist' (DERIVED 69 0 0 DERIVED ()) () (UNKNOWN-FL
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 POINTER) UNKNOWN-ACCESS
(NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0))) PUBLIC (() () () ()) () 0 0
37419158)
86 'csr_matrix' 'sparse_tools' 'csr_matrix' 1 ((DERIVED UNKNOWN-INTENT
UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 POINTER_COMP) (UNKNOWN 0 0 0 UNKNOWN ())
0 0 () () 0 ((153 'sparsity' (DERIVED 69 0 0 DERIVED ()) () (UNKNOWN-FL
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0) UNKNOWN-ACCESS (
STRUCTURE (DERIVED 69 0 0 DERIVED ()) 0 (((NULL (UNKNOWN 0 0 0 UNKNOWN ())
0) ()) ((NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0) ()) ((NULL (UNKNOWN 0 0 0
UNKNOWN ()) 0) ()) (() ()) ((NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0) ()) ((
NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0) ()) ((NULL (UNKNOWN 0 0 0 UNKNOWN ())
0) ()) ((CONSTANT (CHARACTER 1 0 0 CHARACTER (())) 0 101
'                                                                                                     ')
()) ((CONSTANT (LOGICAL 4 0 0 LOGICAL ()) 0 0) ()) ((CONSTANT (LOGICAL 4
0 0 LOGICAL ()) 0 0) ())) ())) (154 'val' (REAL 8 0 0 REAL ()) (1 0
DEFERRED () ()) (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN
0 0 DIMENSION POINTER) UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0))
(155 'ival' (INTEGER 4 0 0 INTEGER ()) (1 0 DEFERRED () ()) (UNKNOWN-FL
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 DIMENSION POINTER)
UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0)) (156 'clone' (
LOGICAL 4 0 0 LOGICAL ()) () (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0) UNKNOWN-ACCESS (CONSTANT (LOGICAL 4 0 0 LOGICAL ())
0 0)) (157 'external_val' (LOGICAL 4 0 0 LOGICAL ()) () (UNKNOWN-FL
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0) UNKNOWN-ACCESS (
CONSTANT (LOGICAL 4 0 0 LOGICAL ()) 0 0)) (158 'inactive' (DERIVED 159 0
0 DERIVED ()) () (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN
0 0 POINTER) UNKNOWN-ACCESS ()) (160 'ksp' (INTEGER 8 0 0 INTEGER ()) ()
(UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 POINTER)
UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0)) (161 'refcount' (
DERIVED 80 0 0 DERIVED ()) () (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0 POINTER) UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN
()) 0)) (162 'name' (CHARACTER 1 0 0 CHARACTER ((CONSTANT (INTEGER 4 0 0
INTEGER ()) 0 '101'))) () (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0) UNKNOWN-ACCESS (CONSTANT (CHARACTER 1 0 0 CHARACTER
(())) 0 101
'                                                                                                     ')))
PUBLIC (() () () ()) () 0 0 66371681)
104 'element_type' 'elements' 'element_type' 1 ((DERIVED UNKNOWN-INTENT
UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 POINTER_COMP) (UNKNOWN 0 0 0 UNKNOWN ())
0 0 () () 0 ((163 'dim' (INTEGER 4 0 0 INTEGER ()) () (UNKNOWN-FL
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0) UNKNOWN-ACCESS ()) (
164 'loc' (INTEGER 4 0 0 INTEGER ()) () (UNKNOWN-FL UNKNOWN-INTENT
UNKNOWN-PROC UNKNOWN UNKNOWN 0 0) UNKNOWN-ACCESS ()) (165 'ngi' (
INTEGER 4 0 0 INTEGER ()) () (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0) UNKNOWN-ACCESS ()) (166 'degree' (INTEGER 4 0 0
INTEGER ()) () (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0
0) UNKNOWN-ACCESS ()) (167 'n' (REAL 8 0 0 REAL ()) (2 0 DEFERRED () ()
() ()) (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DIMENSION POINTER) UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0)) (
168 'dn' (REAL 8 0 0 REAL ()) (3 0 DEFERRED () () () () () ()) (
UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 DIMENSION
POINTER) UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0)) (169 'n_s'
(REAL 8 0 0 REAL ()) (3 0 DEFERRED () () () () () ()) (UNKNOWN-FL
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 DIMENSION POINTER)
UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0)) (170 'dn_s' (REAL 8
0 0 REAL ()) (4 0 DEFERRED () () () () () () () ()) (UNKNOWN-FL
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 DIMENSION POINTER)
UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0)) (171 'spoly' (
DERIVED 172 0 0 DERIVED ()) (2 0 DEFERRED () () () ()) (UNKNOWN-FL
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 DIMENSION POINTER)
UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0)) (173 'dspoly' (
DERIVED 172 0 0 DERIVED ()) (2 0 DEFERRED () () () ()) (UNKNOWN-FL
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 DIMENSION POINTER)
UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0)) (174 'numbering' (
DERIVED 175 0 0 DERIVED ()) () (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0 POINTER) UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN
()) 0)) (176 'quadrature' (DERIVED 105 0 0 DERIVED ()) () (UNKNOWN-FL
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0) UNKNOWN-ACCESS (
STRUCTURE (DERIVED 105 0 0 DERIVED ()) 0 ((() ()) (() ()) (() ()) (() ())
((NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0) ()) ((NULL (UNKNOWN 0 0 0 UNKNOWN ())
0) ()) (() ()) ((NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0) ()) (() ())) ())) (
177 'surface_quadrature' (DERIVED 105 0 0 DERIVED ()) () (UNKNOWN-FL
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 POINTER) UNKNOWN-ACCESS
(NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0)) (178 'superconvergence' (DERIVED
179 0 0 DERIVED ()) () (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN
UNKNOWN 0 0 POINTER) UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0))
(180 'constraints' (DERIVED 181 0 0 DERIVED ()) () (UNKNOWN-FL
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 POINTER) UNKNOWN-ACCESS
(NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0)) (182 'refcount' (DERIVED 80 0 0
DERIVED ()) () (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0
0 POINTER) UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0)) (183 'name'
(CHARACTER 1 0 0 CHARACTER ((CONSTANT (INTEGER 4 0 0 INTEGER ()) 0 '0')))
() (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0)
UNKNOWN-ACCESS ())) PUBLIC (() () () ()) () 0 0 79461029)
159 'logical_array_ptr' 'sparse_tools' 'logical_array_ptr' 1 ((DERIVED
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 POINTER_COMP) (UNKNOWN 0
0 0 UNKNOWN ()) 0 0 () () 0 ((184 'ptr' (LOGICAL 4 0 0 LOGICAL ()) (1 0
DEFERRED () ()) (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN
0 0 DIMENSION POINTER) UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0)))
PUBLIC (() () () ()) () 0 0 30974511)
113 'mesh_faces' 'fields_data_types' 'mesh_faces' 1 ((DERIVED
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 POINTER_COMP) (UNKNOWN 0
0 0 UNKNOWN ()) 0 0 () () 0 ((185 'shape' (DERIVED 104 0 0 DERIVED ()) ()
(UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 POINTER)
UNKNOWN-ACCESS ()) (186 'face_list' (DERIVED 86 0 0 DERIVED ()) () (
UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0)
UNKNOWN-ACCESS (STRUCTURE (DERIVED 86 0 0 DERIVED ()) 0 (((STRUCTURE (
DERIVED 69 0 0 DERIVED ()) 0 (((NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0) ()) (
(NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0) ()) ((NULL (UNKNOWN 0 0 0 UNKNOWN ())
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
())) ())) (187 'face_lno' (INTEGER 4 0 0 INTEGER ()) (1 0 DEFERRED () ())
(UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 DIMENSION
POINTER) UNKNOWN-ACCESS ()) (188 'surface_mesh' (DERIVED 100 0 0 DERIVED
()) () (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0)
UNKNOWN-ACCESS (STRUCTURE (DERIVED 100 0 0 DERIVED ()) 0 ((() ()) ((
CONSTANT (LOGICAL 4 0 0 LOGICAL ()) 0 1) ()) ((STRUCTURE (DERIVED 104 0
0 DERIVED ()) 0 ((() ()) (() ()) (() ()) (() ()) ((NULL (UNKNOWN 0 0 0
UNKNOWN ()) 0) ()) ((NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0) ()) ((NULL (
UNKNOWN 0 0 0 UNKNOWN ()) 0) ()) ((NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0) ())
((NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0) ()) ((NULL (UNKNOWN 0 0 0 UNKNOWN ())
0) ()) ((NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0) ()) ((STRUCTURE (DERIVED 105
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
LOGICAL 4 0 0 LOGICAL ()) 0 0) ())) ())) (189 'surface_node_list' (
INTEGER 4 0 0 INTEGER ()) (1 0 DEFERRED () ()) (UNKNOWN-FL
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 DIMENSION POINTER)
UNKNOWN-ACCESS ()) (190 'face_element_list' (INTEGER 4 0 0 INTEGER ()) (
1 0 DEFERRED () ()) (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN
UNKNOWN 0 0 DIMENSION POINTER) UNKNOWN-ACCESS ()) (191 'boundary_ids' (
INTEGER 4 0 0 INTEGER ()) (1 0 DEFERRED () ()) (UNKNOWN-FL
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 DIMENSION POINTER)
UNKNOWN-ACCESS ()) (192 'coplanar_ids' (INTEGER 4 0 0 INTEGER ()) (1 0
DEFERRED () ()) (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN
0 0 DIMENSION POINTER) UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0))
(193 'dg_surface_mesh' (DERIVED 100 0 0 DERIVED ()) () (UNKNOWN-FL
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 POINTER) UNKNOWN-ACCESS
(NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0)) (194 'has_internal_boundaries' (
LOGICAL 4 0 0 LOGICAL ()) () (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0) UNKNOWN-ACCESS (CONSTANT (LOGICAL 4 0 0 LOGICAL ())
0 0))) PUBLIC (() () () ()) () 0 0 11936185)
115 'mesh_subdomain_mesh' 'fields_data_types' 'mesh_subdomain_mesh' 1 (
(DERIVED UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 POINTER_COMP) (
UNKNOWN 0 0 0 UNKNOWN ()) 0 0 () () 0 ((195 'element_list' (INTEGER 4 0
0 INTEGER ()) (1 0 DEFERRED () ()) (UNKNOWN-FL UNKNOWN-INTENT
UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 DIMENSION POINTER) UNKNOWN-ACCESS ()) (
196 'node_list' (INTEGER 4 0 0 INTEGER ()) (1 0 DEFERRED () ()) (
UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 DIMENSION
POINTER) UNKNOWN-ACCESS ())) PUBLIC (() () () ()) () 0 0 85302181)
172 'polynomial' 'polynomials' 'polynomial' 1 ((DERIVED UNKNOWN-INTENT
UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 POINTER_COMP) (UNKNOWN 0 0 0 UNKNOWN ())
0 0 () () 0 ((197 'coefs' (REAL 8 0 0 REAL ()) (1 0 DEFERRED () ()) (
UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 DIMENSION
POINTER) UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0)) (198 'degree'
(INTEGER 4 0 0 INTEGER ()) () (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0) UNKNOWN-ACCESS (CONSTANT (INTEGER 4 0 0 INTEGER ())
0 '-1'))) PUBLIC (() () () ()) () 0 0 87989236)
71 'real_vector' 'futils' 'real_vector' 1 ((DERIVED UNKNOWN-INTENT
UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 POINTER_COMP) (UNKNOWN 0 0 0 UNKNOWN ())
0 0 () () 0 ((199 'ptr' (REAL 8 0 0 REAL ()) (1 0 DEFERRED () ()) (
UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 DIMENSION
POINTER) UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0))) PUBLIC (()
() () ()) () 0 0 72870256)
80 'refcount_type' 'reference_counting' 'refcount_type' 1 ((DERIVED
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 POINTER_COMP) (UNKNOWN 0
0 0 UNKNOWN ()) 0 0 () () 0 ((200 'prev' (DERIVED 80 0 0 DERIVED ()) ()
(UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 POINTER)
UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0)) (201 'next' (
DERIVED 80 0 0 DERIVED ()) () (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0 POINTER) UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN
()) 0)) (202 'count' (INTEGER 4 0 0 INTEGER ()) () (UNKNOWN-FL
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0) UNKNOWN-ACCESS (
CONSTANT (INTEGER 4 0 0 INTEGER ()) 0 '0')) (203 'id' (INTEGER 4 0 0
INTEGER ()) () (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0
0) UNKNOWN-ACCESS ()) (204 'name' (CHARACTER 1 0 0 CHARACTER ((CONSTANT
(INTEGER 4 0 0 INTEGER ()) 0 '101'))) () (UNKNOWN-FL UNKNOWN-INTENT
UNKNOWN-PROC UNKNOWN UNKNOWN 0 0) UNKNOWN-ACCESS ()) (205 'type' (
CHARACTER 1 0 0 CHARACTER ((CONSTANT (INTEGER 4 0 0 INTEGER ()) 0 '101')))
() (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0)
UNKNOWN-ACCESS ()) (206 'tagged' (LOGICAL 4 0 0 LOGICAL ()) () (
UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0)
UNKNOWN-ACCESS (CONSTANT (LOGICAL 4 0 0 LOGICAL ()) 0 0))) PUBLIC (() ()
() ()) () 0 0 25948645)
51 'scalar_field_pointer' 'fields_data_types' 'scalar_field_pointer' 1 (
(DERIVED UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 POINTER_COMP) (
UNKNOWN 0 0 0 UNKNOWN ()) 0 0 () () 0 ((207 'ptr' (DERIVED 208 0 0
DERIVED ()) () (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0
0 POINTER) UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0))) PUBLIC (
() () () ()) () 0 0 54828570)
179 'superconvergence_type' 'elements' 'superconvergence_type' 1 ((
DERIVED UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 POINTER_COMP) (
UNKNOWN 0 0 0 UNKNOWN ()) 0 0 () () 0 ((209 'nsp' (INTEGER 4 0 0 INTEGER
()) () (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0)
UNKNOWN-ACCESS ()) (210 'l' (REAL 8 0 0 REAL ()) (2 0 DEFERRED () () ()
()) (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DIMENSION POINTER) UNKNOWN-ACCESS ()) (211 'n' (REAL 8 0 0 REAL ()) (2 0
DEFERRED () () () ()) (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN
UNKNOWN 0 0 DIMENSION POINTER) UNKNOWN-ACCESS ()) (212 'dn' (REAL 8 0 0
REAL ()) (3 0 DEFERRED () () () () () ()) (UNKNOWN-FL UNKNOWN-INTENT
UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 DIMENSION POINTER) UNKNOWN-ACCESS ()))
PUBLIC (() () () ()) () 0 0 18282395)
49 'tensor_field_pointer' 'fields_data_types' 'tensor_field_pointer' 1 (
(DERIVED UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 POINTER_COMP) (
UNKNOWN 0 0 0 UNKNOWN ()) 0 0 () () 0 ((213 'ptr' (DERIVED 214 0 0
DERIVED ()) () (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0
0 POINTER) UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0))) PUBLIC (
() () () ()) () 0 0 16722823)
47 'vector_field_pointer' 'fields_data_types' 'vector_field_pointer' 1 (
(DERIVED UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 POINTER_COMP) (
UNKNOWN 0 0 0 UNKNOWN ()) 0 0 () () 0 ((215 'ptr' (DERIVED 216 0 0
DERIVED ()) () (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0
0 POINTER) UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0))) PUBLIC (
() () () ()) () 0 0 74448721)
34 'size' 'fields' 'size' 1 ((PROCEDURE UNKNOWN-INTENT UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0 GENERIC) (UNKNOWN 0 0 0 UNKNOWN ()) 0 0 () () 0 () ()
() 0 0)
149 'petsc_csr_matrix' 'sparse_tools_petsc' 'petsc_csr_matrix' 1 ((
DERIVED UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 POINTER_COMP) (
UNKNOWN 0 0 0 UNKNOWN ()) 0 0 () () 0 ((217 'm' (INTEGER 8 0 0 INTEGER ())
() (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0)
UNKNOWN-ACCESS ()) (218 'row_numbering' (DERIVED 219 0 0 DERIVED ()) ()
(UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0)
UNKNOWN-ACCESS (STRUCTURE (DERIVED 219 0 0 DERIVED ()) 0 (((NULL (
UNKNOWN 0 0 0 UNKNOWN ()) 0) ()) (() ()) (() ()) (() ()) (() ()) ((NULL
(UNKNOWN 0 0 0 UNKNOWN ()) 0) ()) ((NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0) ())
((NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0) ()) ((CONSTANT (CHARACTER 1 0 0
CHARACTER (())) 0 101
'                                                                                                     ')
())) ())) (220 'column_numbering' (DERIVED 219 0 0 DERIVED ()) () (
UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0)
UNKNOWN-ACCESS (STRUCTURE (DERIVED 219 0 0 DERIVED ()) 0 (((NULL (
UNKNOWN 0 0 0 UNKNOWN ()) 0) ()) (() ()) (() ()) (() ()) (() ()) ((NULL
(UNKNOWN 0 0 0 UNKNOWN ()) 0) ()) ((NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0) ())
((NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0) ()) ((CONSTANT (CHARACTER 1 0 0
CHARACTER (())) 0 101
'                                                                                                     ')
())) ())) (221 'row_halo' (DERIVED 67 0 0 DERIVED ()) () (UNKNOWN-FL
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 POINTER) UNKNOWN-ACCESS
(NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0)) (222 'column_halo' (DERIVED 67 0 0
DERIVED ()) () (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0
0 POINTER) UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0)) (223
'refcount' (DERIVED 80 0 0 DERIVED ()) () (UNKNOWN-FL UNKNOWN-INTENT
UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 POINTER) UNKNOWN-ACCESS (NULL (UNKNOWN
0 0 0 UNKNOWN ()) 0)) (224 'name' (CHARACTER 1 0 0 CHARACTER ((CONSTANT
(INTEGER 4 0 0 INTEGER ()) 0 '101'))) () (UNKNOWN-FL UNKNOWN-INTENT
UNKNOWN-PROC UNKNOWN UNKNOWN 0 0) UNKNOWN-ACCESS (CONSTANT (CHARACTER 1
0 0 CHARACTER (())) 0 101
'                                                                                                     '))
(225 'is_assembled' (LOGICAL 4 0 0 LOGICAL ()) () (UNKNOWN-FL
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0) UNKNOWN-ACCESS (
CONSTANT (LOGICAL 4 0 0 LOGICAL ()) 0 0))) PUBLIC (() () () ()) () 0 0
75554753)
181 'constraints_type' 'elements' 'constraints_type' 1 ((DERIVED
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 POINTER_COMP) (UNKNOWN 0
0 0 UNKNOWN ()) 0 0 () () 0 ((226 'type' (INTEGER 4 0 0 INTEGER ()) () (
UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0)
UNKNOWN-ACCESS ()) (227 'dim' (INTEGER 4 0 0 INTEGER ()) () (UNKNOWN-FL
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0) UNKNOWN-ACCESS ()) (
228 'degree' (INTEGER 4 0 0 INTEGER ()) () (UNKNOWN-FL UNKNOWN-INTENT
UNKNOWN-PROC UNKNOWN UNKNOWN 0 0) UNKNOWN-ACCESS ()) (229 'loc' (
INTEGER 4 0 0 INTEGER ()) () (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0) UNKNOWN-ACCESS ()) (230 'n_constraints' (INTEGER 4
0 0 INTEGER ()) () (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN
UNKNOWN 0 0) UNKNOWN-ACCESS ()) (231 'orthogonal' (REAL 8 0 0 REAL ()) (
3 0 DEFERRED () () () () () ()) (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0 DIMENSION POINTER) UNKNOWN-ACCESS (NULL (UNKNOWN 0 0
0 UNKNOWN ()) 0))) PUBLIC (() () () ()) () 0 0 57548971)
175 'ele_numbering_type' 'element_numbering' 'ele_numbering_type' 1 ((
DERIVED UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 POINTER_COMP) (
UNKNOWN 0 0 0 UNKNOWN ()) 0 0 () () 0 ((232 'faces' (INTEGER 4 0 0
INTEGER ()) () (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0
0) UNKNOWN-ACCESS ()) (233 'vertices' (INTEGER 4 0 0 INTEGER ()) () (
UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0)
UNKNOWN-ACCESS ()) (234 'edges' (INTEGER 4 0 0 INTEGER ()) () (
UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0)
UNKNOWN-ACCESS ()) (235 'boundaries' (INTEGER 4 0 0 INTEGER ()) () (
UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0)
UNKNOWN-ACCESS ()) (236 'degree' (INTEGER 4 0 0 INTEGER ()) () (
UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0)
UNKNOWN-ACCESS ()) (237 'dimension' (INTEGER 4 0 0 INTEGER ()) () (
UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0)
UNKNOWN-ACCESS ()) (238 'nodes' (INTEGER 4 0 0 INTEGER ()) () (
UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0)
UNKNOWN-ACCESS ()) (239 'type' (INTEGER 4 0 0 INTEGER ()) () (
UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0)
UNKNOWN-ACCESS (CONSTANT (INTEGER 4 0 0 INTEGER ()) 0 '1')) (240 'family'
(INTEGER 4 0 0 INTEGER ()) () (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0) UNKNOWN-ACCESS ()) (241 'count2number' (INTEGER 4 0
0 INTEGER ()) (3 0 DEFERRED () () () () () ()) (UNKNOWN-FL
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 DIMENSION POINTER)
UNKNOWN-ACCESS ()) (242 'number2count' (INTEGER 4 0 0 INTEGER ()) (2 0
DEFERRED () () () ()) (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN
UNKNOWN 0 0 DIMENSION POINTER) UNKNOWN-ACCESS ()) (243 'boundary_coord'
(INTEGER 4 0 0 INTEGER ()) (1 0 DEFERRED () ()) (UNKNOWN-FL
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 DIMENSION POINTER)
UNKNOWN-ACCESS ()) (244 'boundary_val' (INTEGER 4 0 0 INTEGER ()) (1 0
DEFERRED () ()) (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN
0 0 DIMENSION POINTER) UNKNOWN-ACCESS ())) PUBLIC (() () () ()) () 0 0
96431082)
208 'scalar_field' 'fields_data_types' 'scalar_field' 1 ((DERIVED
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 POINTER_COMP) (UNKNOWN 0
0 0 UNKNOWN ()) 0 0 () () 0 ((245 'val' (REAL 8 0 0 REAL ()) (1 0
DEFERRED () ()) (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN
0 0 DIMENSION POINTER) UNKNOWN-ACCESS ()) (246 'val_stride' (INTEGER 4 0
0 INTEGER ()) () (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN
0 0) UNKNOWN-ACCESS (CONSTANT (INTEGER 4 0 0 INTEGER ()) 0 '1')) (247
'wrapped' (LOGICAL 4 0 0 LOGICAL ()) () (UNKNOWN-FL UNKNOWN-INTENT
UNKNOWN-PROC UNKNOWN UNKNOWN 0 0) UNKNOWN-ACCESS (CONSTANT (LOGICAL 4 0
0 LOGICAL ()) 0 1)) (248 'field_type' (INTEGER 4 0 0 INTEGER ()) () (
UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0)
UNKNOWN-ACCESS (CONSTANT (INTEGER 4 0 0 INTEGER ()) 0 '0')) (249 'bc' (
DERIVED 250 0 0 DERIVED ()) () (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0 POINTER) UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN
()) 0)) (251 'name' (CHARACTER 1 0 0 CHARACTER ((CONSTANT (INTEGER 4 0 0
INTEGER ()) 0 '101'))) () (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0) UNKNOWN-ACCESS ()) (252 'option_path' (CHARACTER 1
0 0 CHARACTER ((CONSTANT (INTEGER 4 0 0 INTEGER ()) 0 '8192'))) () (
UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0)
UNKNOWN-ACCESS (CONSTANT (CHARACTER 1 0 0 CHARACTER (())) 0 8192
'/uninitialised_path/                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            '))
(253 'mesh' (DERIVED 100 0 0 DERIVED ()) () (UNKNOWN-FL UNKNOWN-INTENT
UNKNOWN-PROC UNKNOWN UNKNOWN 0 0) UNKNOWN-ACCESS (STRUCTURE (DERIVED 100
0 0 DERIVED ()) 0 ((() ()) ((CONSTANT (LOGICAL 4 0 0 LOGICAL ()) 0 1) ())
((STRUCTURE (DERIVED 104 0 0 DERIVED ()) 0 ((() ()) (() ()) (() ()) (()
()) ((NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0) ()) ((NULL (UNKNOWN 0 0 0
UNKNOWN ()) 0) ()) ((NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0) ()) ((NULL (
UNKNOWN 0 0 0 UNKNOWN ()) 0) ()) ((NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0) ())
((NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0) ()) ((NULL (UNKNOWN 0 0 0 UNKNOWN ())
0) ()) ((STRUCTURE (DERIVED 105 0 0 DERIVED ()) 0 ((() ()) (() ()) (() ())
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
LOGICAL 4 0 0 LOGICAL ()) 0 0) ())) ())) (254 'refcount' (DERIVED 80 0 0
DERIVED ()) () (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0
0 POINTER) UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0)) (255
'aliased' (LOGICAL 4 0 0 LOGICAL ()) () (UNKNOWN-FL UNKNOWN-INTENT
UNKNOWN-PROC UNKNOWN UNKNOWN 0 0) UNKNOWN-ACCESS (CONSTANT (LOGICAL 4 0
0 LOGICAL ()) 0 0)) (256 'py_locweight' (REAL 8 0 0 REAL ()) (2 0
DEFERRED () () () ()) (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN
UNKNOWN 0 0 DIMENSION POINTER) UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0
UNKNOWN ()) 0)) (257 'py_func' (CHARACTER 1 0 0 CHARACTER ((CONSTANT (
INTEGER 4 0 0 INTEGER ()) 0 '8192'))) () (UNKNOWN-FL UNKNOWN-INTENT
UNKNOWN-PROC UNKNOWN UNKNOWN 0 0) UNKNOWN-ACCESS ()) (258 'py_positions'
(DERIVED 216 0 0 DERIVED ()) () (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0 POINTER) UNKNOWN-ACCESS ()) (259
'py_positions_same_mesh' (LOGICAL 4 0 0 LOGICAL ()) () (UNKNOWN-FL
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0) UNKNOWN-ACCESS ()) (
260 'py_dim' (INTEGER 4 0 0 INTEGER ()) () (UNKNOWN-FL UNKNOWN-INTENT
UNKNOWN-PROC UNKNOWN UNKNOWN 0 0) UNKNOWN-ACCESS ()) (261
'py_positions_shape' (DERIVED 104 0 0 DERIVED ()) () (UNKNOWN-FL
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 POINTER) UNKNOWN-ACCESS
(NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0))) PUBLIC (() () () ()) () 0 0
64912956)
214 'tensor_field' 'fields_data_types' 'tensor_field' 1 ((DERIVED
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 POINTER_COMP) (UNKNOWN 0
0 0 UNKNOWN ()) 0 0 () () 0 ((262 'val' (REAL 8 0 0 REAL ()) (3 0
DEFERRED () () () () () ()) (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0 DIMENSION POINTER) UNKNOWN-ACCESS ()) (263 'wrapped'
(LOGICAL 4 0 0 LOGICAL ()) () (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0) UNKNOWN-ACCESS (CONSTANT (LOGICAL 4 0 0 LOGICAL ())
0 1)) (264 'field_type' (INTEGER 4 0 0 INTEGER ()) () (UNKNOWN-FL
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0) UNKNOWN-ACCESS (
CONSTANT (INTEGER 4 0 0 INTEGER ()) 0 '0')) (265 'name' (CHARACTER 1 0 0
CHARACTER ((CONSTANT (INTEGER 4 0 0 INTEGER ()) 0 '101'))) () (
UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0)
UNKNOWN-ACCESS ()) (266 'dim' (INTEGER 4 0 0 INTEGER ()) (1 0 EXPLICIT (
CONSTANT (INTEGER 4 0 0 INTEGER ()) 0 '1') (CONSTANT (INTEGER 4 0 0
INTEGER ()) 0 '2')) (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN
UNKNOWN 0 0 DIMENSION) UNKNOWN-ACCESS ()) (267 'option_path' (CHARACTER
1 0 0 CHARACTER ((CONSTANT (INTEGER 4 0 0 INTEGER ()) 0 '8192'))) () (
UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0)
UNKNOWN-ACCESS (CONSTANT (CHARACTER 1 0 0 CHARACTER (())) 0 8192
'/uninitialised_path/                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            '))
(268 'mesh' (DERIVED 100 0 0 DERIVED ()) () (UNKNOWN-FL UNKNOWN-INTENT
UNKNOWN-PROC UNKNOWN UNKNOWN 0 0) UNKNOWN-ACCESS (STRUCTURE (DERIVED 100
0 0 DERIVED ()) 0 ((() ()) ((CONSTANT (LOGICAL 4 0 0 LOGICAL ()) 0 1) ())
((STRUCTURE (DERIVED 104 0 0 DERIVED ()) 0 ((() ()) (() ()) (() ()) (()
()) ((NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0) ()) ((NULL (UNKNOWN 0 0 0
UNKNOWN ()) 0) ()) ((NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0) ()) ((NULL (
UNKNOWN 0 0 0 UNKNOWN ()) 0) ()) ((NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0) ())
((NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0) ()) ((NULL (UNKNOWN 0 0 0 UNKNOWN ())
0) ()) ((STRUCTURE (DERIVED 105 0 0 DERIVED ()) 0 ((() ()) (() ()) (() ())
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
LOGICAL 4 0 0 LOGICAL ()) 0 0) ())) ())) (269 'refcount' (DERIVED 80 0 0
DERIVED ()) () (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0
0 POINTER) UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0)) (270
'aliased' (LOGICAL 4 0 0 LOGICAL ()) () (UNKNOWN-FL UNKNOWN-INTENT
UNKNOWN-PROC UNKNOWN UNKNOWN 0 0) UNKNOWN-ACCESS (CONSTANT (LOGICAL 4 0
0 LOGICAL ()) 0 0))) PUBLIC (() () () ()) () 0 0 25110185)
216 'vector_field' 'fields_data_types' 'vector_field' 1 ((DERIVED
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 POINTER_COMP) (UNKNOWN 0
0 0 UNKNOWN ()) 0 0 () () 0 ((271 'val' (REAL 8 0 0 REAL ()) (2 0
DEFERRED () () () ()) (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN
UNKNOWN 0 0 DIMENSION POINTER) UNKNOWN-ACCESS ()) (272 'wrapped' (
LOGICAL 4 0 0 LOGICAL ()) () (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0) UNKNOWN-ACCESS (CONSTANT (LOGICAL 4 0 0 LOGICAL ())
0 1)) (273 'field_type' (INTEGER 4 0 0 INTEGER ()) () (UNKNOWN-FL
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0) UNKNOWN-ACCESS (
CONSTANT (INTEGER 4 0 0 INTEGER ()) 0 '0')) (274 'bc' (DERIVED 275 0 0
DERIVED ()) () (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0
0 POINTER) UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0)) (276 'name'
(CHARACTER 1 0 0 CHARACTER ((CONSTANT (INTEGER 4 0 0 INTEGER ()) 0 '101')))
() (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0)
UNKNOWN-ACCESS ()) (277 'dim' (INTEGER 4 0 0 INTEGER ()) () (UNKNOWN-FL
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0) UNKNOWN-ACCESS ()) (
278 'option_path' (CHARACTER 1 0 0 CHARACTER ((CONSTANT (INTEGER 4 0 0
INTEGER ()) 0 '8192'))) () (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0) UNKNOWN-ACCESS (CONSTANT (CHARACTER 1 0 0 CHARACTER
(())) 0 8192
'/uninitialised_path/                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            '))
(279 'mesh' (DERIVED 100 0 0 DERIVED ()) () (UNKNOWN-FL UNKNOWN-INTENT
UNKNOWN-PROC UNKNOWN UNKNOWN 0 0) UNKNOWN-ACCESS (STRUCTURE (DERIVED 100
0 0 DERIVED ()) 0 ((() ()) ((CONSTANT (LOGICAL 4 0 0 LOGICAL ()) 0 1) ())
((STRUCTURE (DERIVED 104 0 0 DERIVED ()) 0 ((() ()) (() ()) (() ()) (()
()) ((NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0) ()) ((NULL (UNKNOWN 0 0 0
UNKNOWN ()) 0) ()) ((NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0) ()) ((NULL (
UNKNOWN 0 0 0 UNKNOWN ()) 0) ()) ((NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0) ())
((NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0) ()) ((NULL (UNKNOWN 0 0 0 UNKNOWN ())
0) ()) ((STRUCTURE (DERIVED 105 0 0 DERIVED ()) 0 ((() ()) (() ()) (() ())
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
LOGICAL 4 0 0 LOGICAL ()) 0 0) ())) ())) (280 'refcount' (DERIVED 80 0 0
DERIVED ()) () (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0
0 POINTER) UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0)) (281
'aliased' (LOGICAL 4 0 0 LOGICAL ()) () (UNKNOWN-FL UNKNOWN-INTENT
UNKNOWN-PROC UNKNOWN UNKNOWN 0 0) UNKNOWN-ACCESS (CONSTANT (LOGICAL 4 0
0 LOGICAL ()) 0 0)) (282 'picker' (DERIVED 283 0 0 DERIVED ()) () (
UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 POINTER)
UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0))) PUBLIC (() () () ())
() 0 0 43598963)
219 'petsc_numbering_type' 'petsc_tools' 'petsc_numbering_type' 1 ((
DERIVED UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 POINTER_COMP) (
UNKNOWN 0 0 0 UNKNOWN ()) 0 0 () () 0 ((284 'halo' (DERIVED 67 0 0
DERIVED ()) () (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0
0 POINTER) UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0)) (285
'nprivatenodes' (INTEGER 4 0 0 INTEGER ()) () (UNKNOWN-FL UNKNOWN-INTENT
UNKNOWN-PROC UNKNOWN UNKNOWN 0 0) UNKNOWN-ACCESS ()) (286
'universal_length' (INTEGER 4 0 0 INTEGER ()) () (UNKNOWN-FL
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0) UNKNOWN-ACCESS ()) (
287 'offset' (INTEGER 4 0 0 INTEGER ()) () (UNKNOWN-FL UNKNOWN-INTENT
UNKNOWN-PROC UNKNOWN UNKNOWN 0 0) UNKNOWN-ACCESS ()) (288 'gnn2unn' (
INTEGER 4 0 0 INTEGER ()) (2 0 DEFERRED () () () ()) (UNKNOWN-FL
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 DIMENSION POINTER)
UNKNOWN-ACCESS ()) (289 'ghost_nodes' (INTEGER 4 0 0 INTEGER ()) (1 0
DEFERRED () ()) (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN
0 0 DIMENSION POINTER) UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0))
(290 'ghost2unn' (INTEGER 4 0 0 INTEGER ()) (2 0 DEFERRED () () () ()) (
UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 DIMENSION
POINTER) UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0)) (291
'refcount' (DERIVED 80 0 0 DERIVED ()) () (UNKNOWN-FL UNKNOWN-INTENT
UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 POINTER) UNKNOWN-ACCESS (NULL (UNKNOWN
0 0 0 UNKNOWN ()) 0)) (292 'name' (CHARACTER 1 0 0 CHARACTER ((CONSTANT
(INTEGER 4 0 0 INTEGER ()) 0 '101'))) () (UNKNOWN-FL UNKNOWN-INTENT
UNKNOWN-PROC UNKNOWN UNKNOWN 0 0) UNKNOWN-ACCESS (CONSTANT (CHARACTER 1
0 0 CHARACTER (())) 0 101
'                                                                                                     ')))
PUBLIC (() () () ()) () 0 0 16529124)
250 'scalar_boundary_conditions_ptr' 'fields_data_types'
'scalar_boundary_conditions_ptr' 1 ((DERIVED UNKNOWN-INTENT UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0 POINTER_COMP) (UNKNOWN 0 0 0 UNKNOWN ()) 0 0 () () 0
((293 'boundary_condition' (DERIVED 294 0 0 DERIVED ()) (1 0 DEFERRED ()
()) (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DIMENSION POINTER) UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0)))
PUBLIC (() () () ()) () 0 0 47502942)
275 'vector_boundary_conditions_ptr' 'fields_data_types'
'vector_boundary_conditions_ptr' 1 ((DERIVED UNKNOWN-INTENT UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0 POINTER_COMP) (UNKNOWN 0 0 0 UNKNOWN ()) 0 0 () () 0
((295 'boundary_condition' (DERIVED 296 0 0 DERIVED ()) (1 0 DEFERRED ()
()) (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DIMENSION POINTER) UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0)))
PUBLIC (() () () ()) () 0 0 84476181)
283 'picker_ptr' 'picker_data_types' 'picker_ptr' 1 ((DERIVED
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 POINTER_COMP) (UNKNOWN 0
0 0 UNKNOWN ()) 0 0 () () 0 ((297 'ptr' (DERIVED 298 0 0 DERIVED ()) ()
(UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 POINTER)
UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0))) PUBLIC (() () () ())
() 0 0 71120007)
294 'scalar_boundary_condition' 'fields_data_types'
'scalar_boundary_condition' 1 ((DERIVED UNKNOWN-INTENT UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0 POINTER_COMP) (UNKNOWN 0 0 0 UNKNOWN ()) 0 0 () () 0
((299 'name' (CHARACTER 1 0 0 CHARACTER ((CONSTANT (INTEGER 4 0 0
INTEGER ()) 0 '101'))) () (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0) UNKNOWN-ACCESS ()) (300 'type' (CHARACTER 1 0 0
CHARACTER ((CONSTANT (INTEGER 4 0 0 INTEGER ()) 0 '101'))) () (
UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0)
UNKNOWN-ACCESS (CONSTANT (CHARACTER 1 0 0 CHARACTER (())) 0 101
'                                                                                                     '))
(301 'surface_element_list' (INTEGER 4 0 0 INTEGER ()) (1 0 DEFERRED ()
()) (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DIMENSION POINTER) UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0)) (
302 'surface_node_list' (INTEGER 4 0 0 INTEGER ()) (1 0 DEFERRED () ())
(UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 DIMENSION
POINTER) UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0)) (303
'surface_mesh' (DERIVED 100 0 0 DERIVED ()) () (UNKNOWN-FL
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 POINTER) UNKNOWN-ACCESS
()) (304 'surface_fields' (DERIVED 208 0 0 DERIVED ()) (1 0 DEFERRED ()
()) (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DIMENSION POINTER) UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0)) (
305 'option_path' (CHARACTER 1 0 0 CHARACTER ((CONSTANT (INTEGER 4 0 0
INTEGER ()) 0 '8192'))) () (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0) UNKNOWN-ACCESS (CONSTANT (CHARACTER 1 0 0 CHARACTER
(())) 0 8192
'/uninitialised_path/                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            ')))
PUBLIC (() () () ()) () 0 0 50740068)
296 'vector_boundary_condition' 'fields_data_types'
'vector_boundary_condition' 1 ((DERIVED UNKNOWN-INTENT UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0 POINTER_COMP) (UNKNOWN 0 0 0 UNKNOWN ()) 0 0 () () 0
((306 'name' (CHARACTER 1 0 0 CHARACTER ((CONSTANT (INTEGER 4 0 0
INTEGER ()) 0 '101'))) () (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0) UNKNOWN-ACCESS ()) (307 'type' (CHARACTER 1 0 0
CHARACTER ((CONSTANT (INTEGER 4 0 0 INTEGER ()) 0 '101'))) () (
UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0)
UNKNOWN-ACCESS (CONSTANT (CHARACTER 1 0 0 CHARACTER (())) 0 101
'                                                                                                     '))
(308 'applies' (LOGICAL 4 0 0 LOGICAL ()) (1 0 EXPLICIT (CONSTANT (
INTEGER 4 0 0 INTEGER ()) 0 '1') (CONSTANT (INTEGER 4 0 0 INTEGER ()) 0
'3')) (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DIMENSION) UNKNOWN-ACCESS ()) (309 'surface_element_list' (INTEGER 4 0 0
INTEGER ()) (1 0 DEFERRED () ()) (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0 DIMENSION POINTER) UNKNOWN-ACCESS (NULL (UNKNOWN 0 0
0 UNKNOWN ()) 0)) (310 'surface_node_list' (INTEGER 4 0 0 INTEGER ()) (
1 0 DEFERRED () ()) (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN
UNKNOWN 0 0 DIMENSION POINTER) UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0
UNKNOWN ()) 0)) (311 'surface_mesh' (DERIVED 100 0 0 DERIVED ()) () (
UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 POINTER)
UNKNOWN-ACCESS ()) (312 'surface_fields' (DERIVED 216 0 0 DERIVED ()) (
1 0 DEFERRED () ()) (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN
UNKNOWN 0 0 DIMENSION POINTER) UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0
UNKNOWN ()) 0)) (313 'scalar_surface_fields' (DERIVED 208 0 0 DERIVED ())
(1 0 DEFERRED () ()) (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN
UNKNOWN 0 0 DIMENSION POINTER) UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0
UNKNOWN ()) 0)) (314 'option_path' (CHARACTER 1 0 0 CHARACTER ((
CONSTANT (INTEGER 4 0 0 INTEGER ()) 0 '8192'))) () (UNKNOWN-FL
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0) UNKNOWN-ACCESS (
CONSTANT (CHARACTER 1 0 0 CHARACTER (())) 0 8192
'/uninitialised_path/                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            ')))
PUBLIC (() () () ()) () 0 0 42205165)
298 'picker_type' 'picker_data_types' 'picker_type' 1 ((DERIVED
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 POINTER_COMP) (UNKNOWN 0
0 0 UNKNOWN ()) 0 0 () () 0 ((315 'name' (CHARACTER 1 0 0 CHARACTER ((
CONSTANT (INTEGER 4 0 0 INTEGER ()) 0 '101'))) () (UNKNOWN-FL
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0) UNKNOWN-ACCESS ()) (
316 'refcount' (DERIVED 80 0 0 DERIVED ()) () (UNKNOWN-FL UNKNOWN-INTENT
UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 POINTER) UNKNOWN-ACCESS (NULL (UNKNOWN
0 0 0 UNKNOWN ()) 0)) (317 'picker_id' (INTEGER 4 0 0 INTEGER ()) () (
UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0)
UNKNOWN-ACCESS (CONSTANT (INTEGER 4 0 0 INTEGER ()) 0 '0')) (318
'last_mesh_movement' (INTEGER 4 0 0 INTEGER ()) () (UNKNOWN-FL
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0) UNKNOWN-ACCESS (
CONSTANT (INTEGER 4 0 0 INTEGER ()) 0 '0'))) PUBLIC (() () () ()) () 0 0
8821665)
)

('functional' 0 19 'functional_derivative' 0 26)
