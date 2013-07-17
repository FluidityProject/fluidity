GFORTRAN module version '6' created from Slope_limiters_DG.F90 on Mon Jun 10 19:48:04 2013
MD5:1b282f78c707e10d0d09b09dbe42d72b -- If you edit this, you'll get what you deserve.

(() () () () () () () () () () () () () () () () () () () () () () ()
() () () ())

()

()

(('mpi_fortran_argv_null' 2 0 0 'mpi_fortran_argv_null') (
'mpi_fortran_errcodes_ignore' 3 0 0 'mpi_fortran_errcodes_ignore') (
'mpi_fortran_in_place' 4 0 0 'mpi_fortran_in_place') (
'mpi_fortran_bottom' 5 0 0 'mpi_fortran_bottom') ('petscfortran4' 6 0 0
'petscfortran4') ('petscfortran5' 7 0 0 'petscfortran5') ('petscfortran6'
8 0 0 'petscfortran6') ('petscfortran7' 9 0 0 'petscfortran7') (
'petscfortran8' 10 0 0 'petscfortran8') ('petscfortran9' 11 0 0
'petscfortran9') ('mpi_fortran_argvs_null' 12 0 0 'mpi_fortran_argvs_null')
('mpi_fortran_status_ignore' 13 0 0 'mpi_fortran_status_ignore') (
'mpi_fortran_statuses_ignore' 14 0 0 'mpi_fortran_statuses_ignore') (
'petscfortran1' 15 0 0 'petscfortran1') ('petscfortran10' 16 0 0
'petscfortran10') ('petscfortran2' 17 0 0 'petscfortran2') (
'petscfortran3' 18 0 0 'petscfortran3'))

()

()

(19 'limit_fpn' 'slope_limiters_dg' 'limit_fpn' 1 ((PROCEDURE
UNKNOWN-INTENT MODULE-PROC DECL UNKNOWN 0 0 SUBROUTINE) (UNKNOWN 0 0 0
UNKNOWN ()) 20 0 (21 22) () 0 () () () 0 0)
23 'limit_slope_dg' 'slope_limiters_dg' 'limit_slope_dg' 1 ((PROCEDURE
UNKNOWN-INTENT MODULE-PROC DECL UNKNOWN 0 0 SUBROUTINE) (UNKNOWN 0 0 0
UNKNOWN ()) 24 0 (25 26 27 28 29) () 0 () () () 0 0)
30 'limit_vb' 'slope_limiters_dg' 'limit_vb' 1 ((PROCEDURE
UNKNOWN-INTENT MODULE-PROC DECL UNKNOWN 0 0 SUBROUTINE) (UNKNOWN 0 0 0
UNKNOWN ()) 31 0 (32 33) () 0 () () () 0 0)
34 'limiter_cockburn' 'slope_limiters_dg' 'limiter_cockburn' 1 ((
PARAMETER UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN IMPLICIT-SAVE 0 0) (
INTEGER 4 0 0 INTEGER ()) 0 0 () (CONSTANT (INTEGER 4 0 0 INTEGER ()) 0
'2') () 0 () () () 0 0)
35 'limiter_fpn' 'slope_limiters_dg' 'limiter_fpn' 1 ((PARAMETER
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN IMPLICIT-SAVE 0 0) (INTEGER 4 0 0
INTEGER ()) 0 0 () (CONSTANT (INTEGER 4 0 0 INTEGER ()) 0 '4') () 0 () ()
() 0 0)
36 'limiter_hermite_weno' 'slope_limiters_dg' 'limiter_hermite_weno' 1 (
(PARAMETER UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN IMPLICIT-SAVE 0 0) (
INTEGER 4 0 0 INTEGER ()) 0 0 () (CONSTANT (INTEGER 4 0 0 INTEGER ()) 0
'3') () 0 () () () 0 0)
37 'limiter_minimal' 'slope_limiters_dg' 'limiter_minimal' 1 ((
PARAMETER UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN IMPLICIT-SAVE 0 0) (
INTEGER 4 0 0 INTEGER ()) 0 0 () (CONSTANT (INTEGER 4 0 0 INTEGER ()) 0
'1') () 0 () () () 0 0)
38 'limiter_vb' 'slope_limiters_dg' 'limiter_vb' 1 ((PARAMETER
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN IMPLICIT-SAVE 0 0) (INTEGER 4 0 0
INTEGER ()) 0 0 () (CONSTANT (INTEGER 4 0 0 INTEGER ()) 0 '5') () 0 () ()
() 0 0)
25 't' '' 't' 24 ((VARIABLE INOUT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 DUMMY)
(DERIVED 39 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
26 'u' '' 'u' 24 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 DUMMY) (
DERIVED 40 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
27 'x' '' 'x' 24 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 DUMMY) (
DERIVED 40 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
28 'state' '' 'state' 24 ((VARIABLE INOUT UNKNOWN-PROC UNKNOWN UNKNOWN 0
0 DUMMY) (DERIVED 41 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
29 'limiter' '' 'limiter' 24 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN
0 0 DUMMY) (INTEGER 4 0 0 INTEGER ()) 0 0 () () 0 () () () 0 0)
39 'scalar_field' 'fields_data_types' 'scalar_field' 1 ((DERIVED
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 POINTER_COMP) (UNKNOWN 0
0 0 UNKNOWN ()) 0 0 () () 0 ((42 'val' (REAL 8 0 0 REAL ()) (1 0
DEFERRED () ()) (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN
0 0 DIMENSION POINTER) UNKNOWN-ACCESS ()) (43 'val_stride' (INTEGER 4 0
0 INTEGER ()) () (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN
0 0) UNKNOWN-ACCESS (CONSTANT (INTEGER 4 0 0 INTEGER ()) 0 '1')) (44
'wrapped' (LOGICAL 4 0 0 LOGICAL ()) () (UNKNOWN-FL UNKNOWN-INTENT
UNKNOWN-PROC UNKNOWN UNKNOWN 0 0) UNKNOWN-ACCESS (CONSTANT (LOGICAL 4 0
0 LOGICAL ()) 0 1)) (45 'field_type' (INTEGER 4 0 0 INTEGER ()) () (
UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0)
UNKNOWN-ACCESS (CONSTANT (INTEGER 4 0 0 INTEGER ()) 0 '0')) (46 'bc' (
DERIVED 47 0 0 DERIVED ()) () (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0 POINTER) UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN
()) 0)) (48 'name' (CHARACTER 1 0 0 CHARACTER ((CONSTANT (INTEGER 4 0 0
INTEGER ()) 0 '101'))) () (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0) UNKNOWN-ACCESS ()) (49 'option_path' (CHARACTER 1 0
0 CHARACTER ((CONSTANT (INTEGER 4 0 0 INTEGER ()) 0 '8192'))) () (
UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0)
UNKNOWN-ACCESS (CONSTANT (CHARACTER 1 0 0 CHARACTER (())) 0 8192
'/uninitialised_path/                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            '))
(50 'mesh' (DERIVED 51 0 0 DERIVED ()) () (UNKNOWN-FL UNKNOWN-INTENT
UNKNOWN-PROC UNKNOWN UNKNOWN 0 0) UNKNOWN-ACCESS (STRUCTURE (DERIVED 51
0 0 DERIVED ()) 0 ((() ()) ((CONSTANT (LOGICAL 4 0 0 LOGICAL ()) 0 1) ())
((STRUCTURE (DERIVED 52 0 0 DERIVED ()) 0 ((() ()) (() ()) (() ()) (() ())
((NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0) ()) ((NULL (UNKNOWN 0 0 0 UNKNOWN ())
0) ()) ((NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0) ()) ((NULL (UNKNOWN 0 0 0
UNKNOWN ()) 0) ()) ((NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0) ()) ((NULL (
UNKNOWN 0 0 0 UNKNOWN ()) 0) ()) ((NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0) ())
((STRUCTURE (DERIVED 53 0 0 DERIVED ()) 0 ((() ()) (() ()) (() ()) (() ())
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
UNKNOWN ()) 0) ()) ((NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0) ()) ((NULL (
UNKNOWN 0 0 0 UNKNOWN ()) 0) ()) ((CONSTANT (LOGICAL 4 0 0 LOGICAL ()) 0
0) ())) ())) (54 'refcount' (DERIVED 55 0 0 DERIVED ()) () (UNKNOWN-FL
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 POINTER) UNKNOWN-ACCESS
(NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0)) (56 'aliased' (LOGICAL 4 0 0
LOGICAL ()) () (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0
0) UNKNOWN-ACCESS (CONSTANT (LOGICAL 4 0 0 LOGICAL ()) 0 0)) (57
'py_locweight' (REAL 8 0 0 REAL ()) (2 0 DEFERRED () () () ()) (
UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 DIMENSION
POINTER) UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0)) (58 'py_func'
(CHARACTER 1 0 0 CHARACTER ((CONSTANT (INTEGER 4 0 0 INTEGER ()) 0 '8192')))
() (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0)
UNKNOWN-ACCESS ()) (59 'py_positions' (DERIVED 40 0 0 DERIVED ()) () (
UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 POINTER)
UNKNOWN-ACCESS ()) (60 'py_positions_same_mesh' (LOGICAL 4 0 0 LOGICAL ())
() (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0)
UNKNOWN-ACCESS ()) (61 'py_dim' (INTEGER 4 0 0 INTEGER ()) () (
UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0)
UNKNOWN-ACCESS ()) (62 'py_positions_shape' (DERIVED 52 0 0 DERIVED ())
() (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 POINTER)
UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0))) PUBLIC (() () () ())
() 0 0 64912956)
40 'vector_field' 'fields_data_types' 'vector_field' 1 ((DERIVED
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 POINTER_COMP) (UNKNOWN 0
0 0 UNKNOWN ()) 0 0 () () 0 ((63 'val' (REAL 8 0 0 REAL ()) (2 0
DEFERRED () () () ()) (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN
UNKNOWN 0 0 DIMENSION POINTER) UNKNOWN-ACCESS ()) (64 'wrapped' (
LOGICAL 4 0 0 LOGICAL ()) () (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0) UNKNOWN-ACCESS (CONSTANT (LOGICAL 4 0 0 LOGICAL ())
0 1)) (65 'field_type' (INTEGER 4 0 0 INTEGER ()) () (UNKNOWN-FL
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0) UNKNOWN-ACCESS (
CONSTANT (INTEGER 4 0 0 INTEGER ()) 0 '0')) (66 'bc' (DERIVED 67 0 0
DERIVED ()) () (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0
0 POINTER) UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0)) (68 'name'
(CHARACTER 1 0 0 CHARACTER ((CONSTANT (INTEGER 4 0 0 INTEGER ()) 0 '101')))
() (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0)
UNKNOWN-ACCESS ()) (69 'dim' (INTEGER 4 0 0 INTEGER ()) () (UNKNOWN-FL
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0) UNKNOWN-ACCESS ()) (70
'option_path' (CHARACTER 1 0 0 CHARACTER ((CONSTANT (INTEGER 4 0 0
INTEGER ()) 0 '8192'))) () (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0) UNKNOWN-ACCESS (CONSTANT (CHARACTER 1 0 0 CHARACTER
(())) 0 8192
'/uninitialised_path/                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            '))
(71 'mesh' (DERIVED 51 0 0 DERIVED ()) () (UNKNOWN-FL UNKNOWN-INTENT
UNKNOWN-PROC UNKNOWN UNKNOWN 0 0) UNKNOWN-ACCESS (STRUCTURE (DERIVED 51
0 0 DERIVED ()) 0 ((() ()) ((CONSTANT (LOGICAL 4 0 0 LOGICAL ()) 0 1) ())
((STRUCTURE (DERIVED 52 0 0 DERIVED ()) 0 ((() ()) (() ()) (() ()) (() ())
((NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0) ()) ((NULL (UNKNOWN 0 0 0 UNKNOWN ())
0) ()) ((NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0) ()) ((NULL (UNKNOWN 0 0 0
UNKNOWN ()) 0) ()) ((NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0) ()) ((NULL (
UNKNOWN 0 0 0 UNKNOWN ()) 0) ()) ((NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0) ())
((STRUCTURE (DERIVED 53 0 0 DERIVED ()) 0 ((() ()) (() ()) (() ()) (() ())
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
UNKNOWN ()) 0) ()) ((NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0) ()) ((NULL (
UNKNOWN 0 0 0 UNKNOWN ()) 0) ()) ((CONSTANT (LOGICAL 4 0 0 LOGICAL ()) 0
0) ())) ())) (72 'refcount' (DERIVED 55 0 0 DERIVED ()) () (UNKNOWN-FL
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 POINTER) UNKNOWN-ACCESS
(NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0)) (73 'aliased' (LOGICAL 4 0 0
LOGICAL ()) () (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0
0) UNKNOWN-ACCESS (CONSTANT (LOGICAL 4 0 0 LOGICAL ()) 0 0)) (74 'picker'
(DERIVED 75 0 0 DERIVED ()) () (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0 POINTER) UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN
()) 0))) PUBLIC (() () () ()) () 0 0 43598963)
2 'mpi_argv_null' 'mpi_interfaces' 'mpi_argv_null' 1 ((VARIABLE
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 DIMENSION IN_COMMON) (
CHARACTER 1 0 0 CHARACTER ((CONSTANT (INTEGER 4 0 0 INTEGER ()) 0 '1')))
0 0 () (1 0 EXPLICIT (CONSTANT (INTEGER 4 0 0 INTEGER ()) 0 '1') (
CONSTANT (INTEGER 4 0 0 INTEGER ()) 0 '1')) 0 () () () 0 0)
3 'mpi_errcodes_ignore' 'mpi_interfaces' 'mpi_errcodes_ignore' 1 ((
VARIABLE UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 DIMENSION
IN_COMMON) (INTEGER 4 0 0 INTEGER ()) 0 0 () (1 0 EXPLICIT (CONSTANT (
INTEGER 4 0 0 INTEGER ()) 0 '1') (CONSTANT (INTEGER 4 0 0 INTEGER ()) 0
'1')) 0 () () () 0 0)
14 'mpi_statuses_ignore' 'mpi_interfaces' 'mpi_statuses_ignore' 1 ((
VARIABLE UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 IN_COMMON) (
REAL 8 0 0 REAL ()) 0 0 () () 0 () () () 0 0)
17 'petsc_null_integer' 'petscsysdef' 'petsc_null_integer' 1 ((VARIABLE
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 IN_COMMON) (INTEGER 4 0
0 INTEGER ()) 0 0 () () 0 () () () 0 0)
16 'petsc_comm_self' 'petscsysdef' 'petsc_comm_self' 1 ((VARIABLE
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 IN_COMMON) (INTEGER 4 0
0 INTEGER ()) 0 0 () () 0 () () () 0 0)
18 'petsc_null' 'petscsysdef' 'petsc_null' 1 ((VARIABLE UNKNOWN-INTENT
UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 IN_COMMON) (INTEGER 4 0 0 INTEGER ()) 0
0 () () 0 () () () 0 0)
15 'petsc_null_character' 'petscsysdef' 'petsc_null_character' 1 ((
VARIABLE UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN IMPLICIT-SAVE 0 0 IN_COMMON)
(CHARACTER 1 0 0 CHARACTER ((CONSTANT (INTEGER 4 0 0 INTEGER ()) 0 '80')))
0 0 () () 0 () () () 0 0)
5 'mpi_bottom' 'mpi_interfaces' 'mpi_bottom' 1 ((VARIABLE UNKNOWN-INTENT
UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 IN_COMMON) (INTEGER 4 0 0 INTEGER ()) 0
0 () () 0 () () () 0 0)
12 'mpi_argvs_null' 'mpi_interfaces' 'mpi_argvs_null' 1 ((VARIABLE
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 IN_COMMON) (REAL 8 0 0
REAL ()) 0 0 () () 0 () () () 0 0)
13 'mpi_status_ignore' 'mpi_interfaces' 'mpi_status_ignore' 1 ((
VARIABLE UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 DIMENSION
IN_COMMON) (INTEGER 4 0 0 INTEGER ()) 0 0 () (1 0 EXPLICIT (CONSTANT (
INTEGER 4 0 0 INTEGER ()) 0 '1') (CONSTANT (INTEGER 4 0 0 INTEGER ()) 0
'5')) 0 () () () 0 0)
4 'mpi_in_place' 'mpi_interfaces' 'mpi_in_place' 1 ((VARIABLE
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 IN_COMMON) (INTEGER 4 0
0 INTEGER ()) 0 0 () () 0 () () () 0 0)
10 'petsc_null_object' 'petscsysdef' 'petsc_null_object' 1 ((VARIABLE
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 IN_COMMON) (INTEGER 8 0
0 INTEGER ()) 0 0 () () 0 () () () 0 0)
11 'petsc_comm_world' 'petscsysdef' 'petsc_comm_world' 1 ((VARIABLE
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 IN_COMMON) (INTEGER 4 0
0 INTEGER ()) 0 0 () () 0 () () () 0 0)
9 'petsc_null_truth' 'petscsysdef' 'petsc_null_truth' 1 ((VARIABLE
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 IN_COMMON) (LOGICAL 4 0
0 LOGICAL ()) 0 0 () () 0 () () () 0 0)
8 'petsc_null_real' 'petscsysdef' 'petsc_null_real' 1 ((VARIABLE
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 IN_COMMON) (REAL 8 0 0
REAL ()) 0 0 () () 0 () () () 0 0)
7 'petsc_null_double' 'petscsysdef' 'petsc_null_double' 1 ((VARIABLE
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 IN_COMMON) (REAL 8 0 0
REAL ()) 0 0 () () 0 () () () 0 0)
6 'petsc_null_scalar' 'petscsysdef' 'petsc_null_scalar' 1 ((VARIABLE
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 IN_COMMON) (REAL 8 0 0
REAL ()) 0 0 () () 0 () () () 0 0)
33 't' '' 't' 31 ((VARIABLE INOUT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 DUMMY)
(DERIVED 39 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
21 'state' '' 'state' 20 ((VARIABLE INOUT UNKNOWN-PROC UNKNOWN UNKNOWN 0
0 DUMMY) (DERIVED 41 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
22 't' '' 't' 20 ((VARIABLE INOUT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 DUMMY)
(DERIVED 39 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
41 'state_type' 'state_module' 'state_type' 1 ((DERIVED UNKNOWN-INTENT
UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 POINTER_COMP) (UNKNOWN 0 0 0 UNKNOWN ())
0 0 () () 0 ((76 'name' (CHARACTER 1 0 0 CHARACTER ((CONSTANT (INTEGER 4
0 0 INTEGER ()) 0 '101'))) () (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0) UNKNOWN-ACCESS (CONSTANT (CHARACTER 1 0 0 CHARACTER
(())) 0 101
'                                                                                                     '))
(77 'option_path' (CHARACTER 1 0 0 CHARACTER ((CONSTANT (INTEGER 4 0 0
INTEGER ()) 0 '8192'))) () (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0) UNKNOWN-ACCESS (CONSTANT (CHARACTER 1 0 0 CHARACTER
(())) 0 8192
'/uninitialised_path/                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            '))
(78 'vector_names' (CHARACTER 1 0 0 CHARACTER ((CONSTANT (INTEGER 4 0 0
INTEGER ()) 0 '101'))) (1 0 DEFERRED () ()) (UNKNOWN-FL UNKNOWN-INTENT
UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 DIMENSION POINTER) UNKNOWN-ACCESS (
NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0)) (79 'scalar_names' (CHARACTER 1 0 0
CHARACTER ((CONSTANT (INTEGER 4 0 0 INTEGER ()) 0 '101'))) (1 0 DEFERRED
() ()) (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DIMENSION POINTER) UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0)) (
80 'mesh_names' (CHARACTER 1 0 0 CHARACTER ((CONSTANT (INTEGER 4 0 0
INTEGER ()) 0 '101'))) (1 0 DEFERRED () ()) (UNKNOWN-FL UNKNOWN-INTENT
UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 DIMENSION POINTER) UNKNOWN-ACCESS (
NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0)) (81 'halo_names' (CHARACTER 1 0 0
CHARACTER ((CONSTANT (INTEGER 4 0 0 INTEGER ()) 0 '101'))) (1 0 DEFERRED
() ()) (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DIMENSION POINTER) UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0)) (
82 'tensor_names' (CHARACTER 1 0 0 CHARACTER ((CONSTANT (INTEGER 4 0 0
INTEGER ()) 0 '101'))) (1 0 DEFERRED () ()) (UNKNOWN-FL UNKNOWN-INTENT
UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 DIMENSION POINTER) UNKNOWN-ACCESS (
NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0)) (83 'csr_sparsity_names' (CHARACTER
1 0 0 CHARACTER ((CONSTANT (INTEGER 4 0 0 INTEGER ()) 0 '101'))) (1 0
DEFERRED () ()) (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN
0 0 DIMENSION POINTER) UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0))
(84 'csr_matrix_names' (CHARACTER 1 0 0 CHARACTER ((CONSTANT (INTEGER 4
0 0 INTEGER ()) 0 '101'))) (1 0 DEFERRED () ()) (UNKNOWN-FL
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 DIMENSION POINTER)
UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0)) (85
'block_csr_matrix_names' (CHARACTER 1 0 0 CHARACTER ((CONSTANT (INTEGER
4 0 0 INTEGER ()) 0 '101'))) (1 0 DEFERRED () ()) (UNKNOWN-FL
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 DIMENSION POINTER)
UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0)) (86
'petsc_csr_matrix_names' (CHARACTER 1 0 0 CHARACTER ((CONSTANT (INTEGER
4 0 0 INTEGER ()) 0 '101'))) (1 0 DEFERRED () ()) (UNKNOWN-FL
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 DIMENSION POINTER)
UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0)) (87 'vector_fields'
(DERIVED 88 0 0 DERIVED ()) (1 0 DEFERRED () ()) (UNKNOWN-FL
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 DIMENSION POINTER)
UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0)) (89 'tensor_fields'
(DERIVED 90 0 0 DERIVED ()) (1 0 DEFERRED () ()) (UNKNOWN-FL
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 DIMENSION POINTER)
UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0)) (91 'scalar_fields'
(DERIVED 92 0 0 DERIVED ()) (1 0 DEFERRED () ()) (UNKNOWN-FL
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 DIMENSION POINTER)
UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0)) (93 'meshes' (
DERIVED 94 0 0 DERIVED ()) (1 0 DEFERRED () ()) (UNKNOWN-FL
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 DIMENSION POINTER)
UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0)) (95 'halos' (
DERIVED 96 0 0 DERIVED ()) (1 0 DEFERRED () ()) (UNKNOWN-FL
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 DIMENSION POINTER)
UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0)) (97 'csr_sparsities'
(DERIVED 98 0 0 DERIVED ()) (1 0 DEFERRED () ()) (UNKNOWN-FL
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 DIMENSION POINTER)
UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0)) (99 'csr_matrices' (
DERIVED 100 0 0 DERIVED ()) (1 0 DEFERRED () ()) (UNKNOWN-FL
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 DIMENSION POINTER)
UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0)) (101
'block_csr_matrices' (DERIVED 102 0 0 DERIVED ()) (1 0 DEFERRED () ()) (
UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 DIMENSION
POINTER) UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0)) (103
'petsc_csr_matrices' (DERIVED 104 0 0 DERIVED ()) (1 0 DEFERRED () ()) (
UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 DIMENSION
POINTER) UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0))) PUBLIC (()
() () ()) () 0 0 7575341)
102 'block_csr_matrix_pointer' 'sparse_tools' 'block_csr_matrix_pointer'
1 ((DERIVED UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 POINTER_COMP)
(UNKNOWN 0 0 0 UNKNOWN ()) 0 0 () () 0 ((105 'ptr' (DERIVED 106 0 0
DERIVED ()) () (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0
0 POINTER) UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0))) PUBLIC (
() () () ()) () 0 0 83606449)
100 'csr_matrix_pointer' 'sparse_tools' 'csr_matrix_pointer' 1 ((
DERIVED UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 POINTER_COMP) (
UNKNOWN 0 0 0 UNKNOWN ()) 0 0 () () 0 ((107 'ptr' (DERIVED 108 0 0
DERIVED ()) () (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0
0 POINTER) UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0))) PUBLIC (
() () () ()) () 0 0 82770239)
98 'csr_sparsity_pointer' 'sparse_tools' 'csr_sparsity_pointer' 1 ((
DERIVED UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 POINTER_COMP) (
UNKNOWN 0 0 0 UNKNOWN ()) 0 0 () () 0 ((109 'ptr' (DERIVED 110 0 0
DERIVED ()) () (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0
0 POINTER) UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0))) PUBLIC (
() () () ()) () 0 0 9687815)
52 'element_type' 'elements' 'element_type' 1 ((DERIVED UNKNOWN-INTENT
UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 POINTER_COMP) (UNKNOWN 0 0 0 UNKNOWN ())
0 0 () () 0 ((111 'dim' (INTEGER 4 0 0 INTEGER ()) () (UNKNOWN-FL
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0) UNKNOWN-ACCESS ()) (
112 'loc' (INTEGER 4 0 0 INTEGER ()) () (UNKNOWN-FL UNKNOWN-INTENT
UNKNOWN-PROC UNKNOWN UNKNOWN 0 0) UNKNOWN-ACCESS ()) (113 'ngi' (
INTEGER 4 0 0 INTEGER ()) () (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0) UNKNOWN-ACCESS ()) (114 'degree' (INTEGER 4 0 0
INTEGER ()) () (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0
0) UNKNOWN-ACCESS ()) (115 'n' (REAL 8 0 0 REAL ()) (2 0 DEFERRED () ()
() ()) (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DIMENSION POINTER) UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0)) (
116 'dn' (REAL 8 0 0 REAL ()) (3 0 DEFERRED () () () () () ()) (
UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 DIMENSION
POINTER) UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0)) (117 'n_s'
(REAL 8 0 0 REAL ()) (3 0 DEFERRED () () () () () ()) (UNKNOWN-FL
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 DIMENSION POINTER)
UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0)) (118 'dn_s' (REAL 8
0 0 REAL ()) (4 0 DEFERRED () () () () () () () ()) (UNKNOWN-FL
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 DIMENSION POINTER)
UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0)) (119 'spoly' (
DERIVED 120 0 0 DERIVED ()) (2 0 DEFERRED () () () ()) (UNKNOWN-FL
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 DIMENSION POINTER)
UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0)) (121 'dspoly' (
DERIVED 120 0 0 DERIVED ()) (2 0 DEFERRED () () () ()) (UNKNOWN-FL
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 DIMENSION POINTER)
UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0)) (122 'numbering' (
DERIVED 123 0 0 DERIVED ()) () (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0 POINTER) UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN
()) 0)) (124 'quadrature' (DERIVED 53 0 0 DERIVED ()) () (UNKNOWN-FL
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0) UNKNOWN-ACCESS (
STRUCTURE (DERIVED 53 0 0 DERIVED ()) 0 ((() ()) (() ()) (() ()) (() ())
((NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0) ()) ((NULL (UNKNOWN 0 0 0 UNKNOWN ())
0) ()) (() ()) ((NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0) ()) (() ())) ())) (
125 'surface_quadrature' (DERIVED 53 0 0 DERIVED ()) () (UNKNOWN-FL
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 POINTER) UNKNOWN-ACCESS
(NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0)) (126 'superconvergence' (DERIVED
127 0 0 DERIVED ()) () (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN
UNKNOWN 0 0 POINTER) UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0))
(128 'constraints' (DERIVED 129 0 0 DERIVED ()) () (UNKNOWN-FL
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 POINTER) UNKNOWN-ACCESS
(NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0)) (130 'refcount' (DERIVED 55 0 0
DERIVED ()) () (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0
0 POINTER) UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0)) (131 'name'
(CHARACTER 1 0 0 CHARACTER ((CONSTANT (INTEGER 4 0 0 INTEGER ()) 0 '0')))
() (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0)
UNKNOWN-ACCESS ())) PUBLIC (() () () ()) () 0 0 79461029)
94 'mesh_pointer' 'fields_data_types' 'mesh_pointer' 1 ((DERIVED
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 POINTER_COMP) (UNKNOWN 0
0 0 UNKNOWN ()) 0 0 () () 0 ((132 'ptr' (DERIVED 51 0 0 DERIVED ()) () (
UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 POINTER)
UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0))) PUBLIC (() () () ())
() 0 0 65432384)
51 'mesh_type' 'fields_data_types' 'mesh_type' 1 ((DERIVED
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 POINTER_COMP) (UNKNOWN 0
0 0 UNKNOWN ()) 0 0 () () 0 ((133 'ndglno' (INTEGER 4 0 0 INTEGER ()) (
1 0 DEFERRED () ()) (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN
UNKNOWN 0 0 DIMENSION POINTER) UNKNOWN-ACCESS ()) (134 'wrapped' (
LOGICAL 4 0 0 LOGICAL ()) () (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0) UNKNOWN-ACCESS (CONSTANT (LOGICAL 4 0 0 LOGICAL ())
0 1)) (135 'shape' (DERIVED 52 0 0 DERIVED ()) () (UNKNOWN-FL
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0) UNKNOWN-ACCESS (
STRUCTURE (DERIVED 52 0 0 DERIVED ()) 0 ((() ()) (() ()) (() ()) (() ())
((NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0) ()) ((NULL (UNKNOWN 0 0 0 UNKNOWN ())
0) ()) ((NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0) ()) ((NULL (UNKNOWN 0 0 0
UNKNOWN ()) 0) ()) ((NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0) ()) ((NULL (
UNKNOWN 0 0 0 UNKNOWN ()) 0) ()) ((NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0) ())
((STRUCTURE (DERIVED 53 0 0 DERIVED ()) 0 ((() ()) (() ()) (() ()) (() ())
((NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0) ()) ((NULL (UNKNOWN 0 0 0 UNKNOWN ())
0) ()) (() ()) ((NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0) ()) (() ())) ()) ())
((NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0) ()) ((NULL (UNKNOWN 0 0 0 UNKNOWN ())
0) ()) ((NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0) ()) ((NULL (UNKNOWN 0 0 0
UNKNOWN ()) 0) ()) (() ())) ())) (136 'elements' (INTEGER 4 0 0 INTEGER
()) () (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0)
UNKNOWN-ACCESS ()) (137 'nodes' (INTEGER 4 0 0 INTEGER ()) () (
UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0)
UNKNOWN-ACCESS ()) (138 'name' (CHARACTER 1 0 0 CHARACTER ((CONSTANT (
INTEGER 4 0 0 INTEGER ()) 0 '101'))) () (UNKNOWN-FL UNKNOWN-INTENT
UNKNOWN-PROC UNKNOWN UNKNOWN 0 0) UNKNOWN-ACCESS ()) (139 'option_path'
(CHARACTER 1 0 0 CHARACTER ((CONSTANT (INTEGER 4 0 0 INTEGER ()) 0 '8192')))
() (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0)
UNKNOWN-ACCESS (CONSTANT (CHARACTER 1 0 0 CHARACTER (())) 0 8192
'/uninitialised_path/                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            '))
(140 'continuity' (INTEGER 4 0 0 INTEGER ()) () (UNKNOWN-FL
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0) UNKNOWN-ACCESS (
CONSTANT (INTEGER 4 0 0 INTEGER ()) 0 '0')) (141 'refcount' (DERIVED 55
0 0 DERIVED ()) () (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN
UNKNOWN 0 0 POINTER) UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0))
(142 'faces' (DERIVED 143 0 0 DERIVED ()) () (UNKNOWN-FL UNKNOWN-INTENT
UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 POINTER) UNKNOWN-ACCESS (NULL (UNKNOWN
0 0 0 UNKNOWN ()) 0)) (144 'subdomain_mesh' (DERIVED 145 0 0 DERIVED ())
() (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 POINTER)
UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0)) (146 'adj_lists' (
DERIVED 147 0 0 DERIVED ()) () (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0 POINTER) UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN
()) 0)) (148 'columns' (INTEGER 4 0 0 INTEGER ()) (1 0 DEFERRED () ()) (
UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 DIMENSION
POINTER) UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0)) (149
'element_columns' (INTEGER 4 0 0 INTEGER ()) (1 0 DEFERRED () ()) (
UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 DIMENSION
POINTER) UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0)) (150
'region_ids' (INTEGER 4 0 0 INTEGER ()) (1 0 DEFERRED () ()) (
UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 DIMENSION
POINTER) UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0)) (151 'halos'
(DERIVED 152 0 0 DERIVED ()) (1 0 DEFERRED () ()) (UNKNOWN-FL
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 DIMENSION POINTER)
UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0)) (153 'element_halos'
(DERIVED 152 0 0 DERIVED ()) (1 0 DEFERRED () ()) (UNKNOWN-FL
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 DIMENSION POINTER)
UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0)) (154 'colourings' (
DERIVED 155 0 0 DERIVED ()) (1 0 DEFERRED () ()) (UNKNOWN-FL
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 DIMENSION POINTER)
UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0)) (156 'periodic' (
LOGICAL 4 0 0 LOGICAL ()) () (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0) UNKNOWN-ACCESS (CONSTANT (LOGICAL 4 0 0 LOGICAL ())
0 0))) PUBLIC (() () () ()) () 0 0 2572791)
120 'polynomial' 'polynomials' 'polynomial' 1 ((DERIVED UNKNOWN-INTENT
UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 POINTER_COMP) (UNKNOWN 0 0 0 UNKNOWN ())
0 0 () () 0 ((157 'coefs' (REAL 8 0 0 REAL ()) (1 0 DEFERRED () ()) (
UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 DIMENSION
POINTER) UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0)) (158 'degree'
(INTEGER 4 0 0 INTEGER ()) () (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0) UNKNOWN-ACCESS (CONSTANT (INTEGER 4 0 0 INTEGER ())
0 '-1'))) PUBLIC (() () () ()) () 0 0 87989236)
53 'quadrature_type' 'quadrature' 'quadrature_type' 1 ((DERIVED
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
166 'refcount' (DERIVED 55 0 0 DERIVED ()) () (UNKNOWN-FL UNKNOWN-INTENT
UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 POINTER) UNKNOWN-ACCESS (NULL (UNKNOWN
0 0 0 UNKNOWN ()) 0)) (167 'family' (INTEGER 4 0 0 INTEGER ()) () (
UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0)
UNKNOWN-ACCESS ())) PUBLIC (() () () ()) () 0 0 59837722)
55 'refcount_type' 'reference_counting' 'refcount_type' 1 ((DERIVED
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 POINTER_COMP) (UNKNOWN 0
0 0 UNKNOWN ()) 0 0 () () 0 ((168 'prev' (DERIVED 55 0 0 DERIVED ()) ()
(UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 POINTER)
UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0)) (169 'next' (
DERIVED 55 0 0 DERIVED ()) () (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0 POINTER) UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN
()) 0)) (170 'count' (INTEGER 4 0 0 INTEGER ()) () (UNKNOWN-FL
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0) UNKNOWN-ACCESS (
CONSTANT (INTEGER 4 0 0 INTEGER ()) 0 '0')) (171 'id' (INTEGER 4 0 0
INTEGER ()) () (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0
0) UNKNOWN-ACCESS ()) (172 'name' (CHARACTER 1 0 0 CHARACTER ((CONSTANT
(INTEGER 4 0 0 INTEGER ()) 0 '101'))) () (UNKNOWN-FL UNKNOWN-INTENT
UNKNOWN-PROC UNKNOWN UNKNOWN 0 0) UNKNOWN-ACCESS ()) (173 'type' (
CHARACTER 1 0 0 CHARACTER ((CONSTANT (INTEGER 4 0 0 INTEGER ()) 0 '101')))
() (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0)
UNKNOWN-ACCESS ()) (174 'tagged' (LOGICAL 4 0 0 LOGICAL ()) () (
UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0)
UNKNOWN-ACCESS (CONSTANT (LOGICAL 4 0 0 LOGICAL ()) 0 0))) PUBLIC (() ()
() ()) () 0 0 25948645)
47 'scalar_boundary_conditions_ptr' 'fields_data_types'
'scalar_boundary_conditions_ptr' 1 ((DERIVED UNKNOWN-INTENT UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0 POINTER_COMP) (UNKNOWN 0 0 0 UNKNOWN ()) 0 0 () () 0
((175 'boundary_condition' (DERIVED 176 0 0 DERIVED ()) (1 0 DEFERRED ()
()) (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DIMENSION POINTER) UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0)))
PUBLIC (() () () ()) () 0 0 47502942)
92 'scalar_field_pointer' 'fields_data_types' 'scalar_field_pointer' 1 (
(DERIVED UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 POINTER_COMP) (
UNKNOWN 0 0 0 UNKNOWN ()) 0 0 () () 0 ((177 'ptr' (DERIVED 39 0 0
DERIVED ()) () (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0
0 POINTER) UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0))) PUBLIC (
() () () ()) () 0 0 54828570)
127 'superconvergence_type' 'elements' 'superconvergence_type' 1 ((
DERIVED UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 POINTER_COMP) (
UNKNOWN 0 0 0 UNKNOWN ()) 0 0 () () 0 ((178 'nsp' (INTEGER 4 0 0 INTEGER
()) () (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0)
UNKNOWN-ACCESS ()) (179 'l' (REAL 8 0 0 REAL ()) (2 0 DEFERRED () () ()
()) (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DIMENSION POINTER) UNKNOWN-ACCESS ()) (180 'n' (REAL 8 0 0 REAL ()) (2 0
DEFERRED () () () ()) (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN
UNKNOWN 0 0 DIMENSION POINTER) UNKNOWN-ACCESS ()) (181 'dn' (REAL 8 0 0
REAL ()) (3 0 DEFERRED () () () () () ()) (UNKNOWN-FL UNKNOWN-INTENT
UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 DIMENSION POINTER) UNKNOWN-ACCESS ()))
PUBLIC (() () () ()) () 0 0 18282395)
90 'tensor_field_pointer' 'fields_data_types' 'tensor_field_pointer' 1 (
(DERIVED UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 POINTER_COMP) (
UNKNOWN 0 0 0 UNKNOWN ()) 0 0 () () 0 ((182 'ptr' (DERIVED 183 0 0
DERIVED ()) () (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0
0 POINTER) UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0))) PUBLIC (
() () () ()) () 0 0 16722823)
67 'vector_boundary_conditions_ptr' 'fields_data_types'
'vector_boundary_conditions_ptr' 1 ((DERIVED UNKNOWN-INTENT UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0 POINTER_COMP) (UNKNOWN 0 0 0 UNKNOWN ()) 0 0 () () 0
((184 'boundary_condition' (DERIVED 185 0 0 DERIVED ()) (1 0 DEFERRED ()
()) (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DIMENSION POINTER) UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0)))
PUBLIC (() () () ()) () 0 0 84476181)
88 'vector_field_pointer' 'fields_data_types' 'vector_field_pointer' 1 (
(DERIVED UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 POINTER_COMP) (
UNKNOWN 0 0 0 UNKNOWN ()) 0 0 () () 0 ((186 'ptr' (DERIVED 40 0 0
DERIVED ()) () (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0
0 POINTER) UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0))) PUBLIC (
() () () ()) () 0 0 74448721)
152 'halo_type' 'halo_data_types' 'halo_type' 1 ((DERIVED UNKNOWN-INTENT
UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 POINTER_COMP) (UNKNOWN 0 0 0 UNKNOWN ())
0 0 () () 0 ((187 'name' (CHARACTER 1 0 0 CHARACTER ((CONSTANT (INTEGER
4 0 0 INTEGER ()) 0 '101'))) () (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0) UNKNOWN-ACCESS ()) (188 'refcount' (DERIVED 55 0 0
DERIVED ()) () (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0
0 POINTER) UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0)) (189
'data_type' (INTEGER 4 0 0 INTEGER ()) () (UNKNOWN-FL UNKNOWN-INTENT
UNKNOWN-PROC UNKNOWN UNKNOWN 0 0) UNKNOWN-ACCESS (CONSTANT (INTEGER 4 0
0 INTEGER ()) 0 '0')) (190 'ordering_scheme' (INTEGER 4 0 0 INTEGER ())
() (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0)
UNKNOWN-ACCESS (CONSTANT (INTEGER 4 0 0 INTEGER ()) 0 '0')) (191
'communicator' (INTEGER 4 0 0 INTEGER ()) () (UNKNOWN-FL UNKNOWN-INTENT
UNKNOWN-PROC UNKNOWN UNKNOWN 0 0) UNKNOWN-ACCESS ()) (192 'nprocs' (
INTEGER 4 0 0 INTEGER ()) () (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0) UNKNOWN-ACCESS (CONSTANT (INTEGER 4 0 0 INTEGER ())
0 '0')) (193 'sends' (DERIVED 194 0 0 DERIVED ()) (1 0 DEFERRED () ()) (
UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 DIMENSION
POINTER) UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0)) (195
'receives' (DERIVED 194 0 0 DERIVED ()) (1 0 DEFERRED () ()) (
UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 DIMENSION
POINTER) UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0)) (196
'nowned_nodes' (INTEGER 4 0 0 INTEGER ()) () (UNKNOWN-FL UNKNOWN-INTENT
UNKNOWN-PROC UNKNOWN UNKNOWN 0 0) UNKNOWN-ACCESS (CONSTANT (INTEGER 4 0
0 INTEGER ()) 0 '-1')) (197 'owners' (INTEGER 4 0 0 INTEGER ()) (1 0
DEFERRED () ()) (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN
0 0 DIMENSION POINTER) UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0))
(198 'unn_count' (INTEGER 4 0 0 INTEGER ()) () (UNKNOWN-FL
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0) UNKNOWN-ACCESS (
CONSTANT (INTEGER 4 0 0 INTEGER ()) 0 '-1')) (199 'owned_nodes_unn_base'
(INTEGER 4 0 0 INTEGER ()) (1 0 DEFERRED () ()) (UNKNOWN-FL
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 DIMENSION POINTER)
UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0)) (200
'my_owned_nodes_unn_base' (INTEGER 4 0 0 INTEGER ()) () (UNKNOWN-FL
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0) UNKNOWN-ACCESS (
CONSTANT (INTEGER 4 0 0 INTEGER ()) 0 '-1')) (201 'receives_gnn_to_unn'
(INTEGER 4 0 0 INTEGER ()) (1 0 DEFERRED () ()) (UNKNOWN-FL
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 DIMENSION POINTER)
UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0)) (202 'gnn_to_unn' (
INTEGER 4 0 0 INTEGER ()) (1 0 DEFERRED () ()) (UNKNOWN-FL
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 DIMENSION POINTER)
UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0))) PUBLIC (() () () ())
() 0 0 63994053)
96 'halo_pointer' 'halo_data_types' 'halo_pointer' 1 ((DERIVED
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 POINTER_COMP) (UNKNOWN 0
0 0 UNKNOWN ()) 0 0 () () 0 ((203 'ptr' (DERIVED 152 0 0 DERIVED ()) ()
(UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 POINTER)
UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0))) PUBLIC (() () () ())
() 0 0 72023282)
104 'petsc_csr_matrix_pointer' 'sparse_tools_petsc'
'petsc_csr_matrix_pointer' 1 ((DERIVED UNKNOWN-INTENT UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0 POINTER_COMP) (UNKNOWN 0 0 0 UNKNOWN ()) 0 0 () () 0
((204 'ptr' (DERIVED 205 0 0 DERIVED ()) () (UNKNOWN-FL UNKNOWN-INTENT
UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 POINTER) UNKNOWN-ACCESS (NULL (UNKNOWN
0 0 0 UNKNOWN ()) 0))) PUBLIC (() () () ()) () 0 0 85445535)
205 'petsc_csr_matrix' 'sparse_tools_petsc' 'petsc_csr_matrix' 1 ((
DERIVED UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 POINTER_COMP) (
UNKNOWN 0 0 0 UNKNOWN ()) 0 0 () () 0 ((206 'm' (INTEGER 8 0 0 INTEGER ())
() (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0)
UNKNOWN-ACCESS ()) (207 'row_numbering' (DERIVED 208 0 0 DERIVED ()) ()
(UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0)
UNKNOWN-ACCESS (STRUCTURE (DERIVED 208 0 0 DERIVED ()) 0 (((NULL (
UNKNOWN 0 0 0 UNKNOWN ()) 0) ()) (() ()) (() ()) (() ()) (() ()) ((NULL
(UNKNOWN 0 0 0 UNKNOWN ()) 0) ()) ((NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0) ())
((NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0) ()) ((CONSTANT (CHARACTER 1 0 0
CHARACTER (())) 0 101
'                                                                                                     ')
())) ())) (209 'column_numbering' (DERIVED 208 0 0 DERIVED ()) () (
UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0)
UNKNOWN-ACCESS (STRUCTURE (DERIVED 208 0 0 DERIVED ()) 0 (((NULL (
UNKNOWN 0 0 0 UNKNOWN ()) 0) ()) (() ()) (() ()) (() ()) (() ()) ((NULL
(UNKNOWN 0 0 0 UNKNOWN ()) 0) ()) ((NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0) ())
((NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0) ()) ((CONSTANT (CHARACTER 1 0 0
CHARACTER (())) 0 101
'                                                                                                     ')
())) ())) (210 'row_halo' (DERIVED 152 0 0 DERIVED ()) () (UNKNOWN-FL
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 POINTER) UNKNOWN-ACCESS
(NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0)) (211 'column_halo' (DERIVED 152 0 0
DERIVED ()) () (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0
0 POINTER) UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0)) (212
'refcount' (DERIVED 55 0 0 DERIVED ()) () (UNKNOWN-FL UNKNOWN-INTENT
UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 POINTER) UNKNOWN-ACCESS (NULL (UNKNOWN
0 0 0 UNKNOWN ()) 0)) (213 'name' (CHARACTER 1 0 0 CHARACTER ((CONSTANT
(INTEGER 4 0 0 INTEGER ()) 0 '101'))) () (UNKNOWN-FL UNKNOWN-INTENT
UNKNOWN-PROC UNKNOWN UNKNOWN 0 0) UNKNOWN-ACCESS (CONSTANT (CHARACTER 1
0 0 CHARACTER (())) 0 101
'                                                                                                     '))
(214 'is_assembled' (LOGICAL 4 0 0 LOGICAL ()) () (UNKNOWN-FL
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0) UNKNOWN-ACCESS (
CONSTANT (LOGICAL 4 0 0 LOGICAL ()) 0 0))) PUBLIC (() () () ()) () 0 0
75554753)
147 'adjacency_cache' 'fields_data_types' 'adjacency_cache' 1 ((DERIVED
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 POINTER_COMP) (UNKNOWN 0
0 0 UNKNOWN ()) 0 0 () () 0 ((215 'nnlist' (DERIVED 110 0 0 DERIVED ())
() (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 POINTER)
UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0)) (216 'nelist' (
DERIVED 110 0 0 DERIVED ()) () (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0 POINTER) UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN
()) 0)) (217 'eelist' (DERIVED 110 0 0 DERIVED ()) () (UNKNOWN-FL
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 POINTER) UNKNOWN-ACCESS
(NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0))) PUBLIC (() () () ()) () 0 0
37419158)
106 'block_csr_matrix' 'sparse_tools' 'block_csr_matrix' 1 ((DERIVED
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 POINTER_COMP) (UNKNOWN 0
0 0 UNKNOWN ()) 0 0 () () 0 ((218 'sparsity' (DERIVED 110 0 0 DERIVED ())
() (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0)
UNKNOWN-ACCESS (STRUCTURE (DERIVED 110 0 0 DERIVED ()) 0 (((NULL (
UNKNOWN 0 0 0 UNKNOWN ()) 0) ()) ((NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0) ())
((NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0) ()) (() ()) ((NULL (UNKNOWN 0 0 0
UNKNOWN ()) 0) ()) ((NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0) ()) ((NULL (
UNKNOWN 0 0 0 UNKNOWN ()) 0) ()) ((CONSTANT (CHARACTER 1 0 0 CHARACTER (
())) 0 101
'                                                                                                     ')
()) ((CONSTANT (LOGICAL 4 0 0 LOGICAL ()) 0 0) ()) ((CONSTANT (LOGICAL 4
0 0 LOGICAL ()) 0 0) ())) ())) (219 'val' (DERIVED 220 0 0 DERIVED ()) (
2 0 DEFERRED () () () ()) (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0 DIMENSION POINTER) UNKNOWN-ACCESS (NULL (UNKNOWN 0 0
0 UNKNOWN ()) 0)) (221 'contiguous_val' (REAL 8 0 0 REAL ()) (1 0
DEFERRED () ()) (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN
0 0 DIMENSION POINTER) UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0))
(222 'ival' (DERIVED 194 0 0 DERIVED ()) (2 0 DEFERRED () () () ()) (
UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 DIMENSION
POINTER) UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0)) (223 'blocks'
(INTEGER 4 0 0 INTEGER ()) (1 0 EXPLICIT (CONSTANT (INTEGER 4 0 0
INTEGER ()) 0 '1') (CONSTANT (INTEGER 4 0 0 INTEGER ()) 0 '2')) (
UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 DIMENSION)
UNKNOWN-ACCESS (ARRAY (INTEGER 4 0 0 INTEGER ()) 1 (((CONSTANT (INTEGER
4 0 0 INTEGER ()) 0 '0') ()) ((CONSTANT (INTEGER 4 0 0 INTEGER ()) 0 '0')
())) ('2'))) (224 'clone' (LOGICAL 4 0 0 LOGICAL ()) () (UNKNOWN-FL
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0) UNKNOWN-ACCESS (
CONSTANT (LOGICAL 4 0 0 LOGICAL ()) 0 0)) (225 'external_val' (LOGICAL 4
0 0 LOGICAL ()) () (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN
UNKNOWN 0 0) UNKNOWN-ACCESS (CONSTANT (LOGICAL 4 0 0 LOGICAL ()) 0 0)) (
226 'columns' (INTEGER 4 0 0 INTEGER ()) () (UNKNOWN-FL UNKNOWN-INTENT
UNKNOWN-PROC UNKNOWN UNKNOWN 0 0) UNKNOWN-ACCESS ()) (227 'refcount' (
DERIVED 55 0 0 DERIVED ()) () (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0 POINTER) UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN
()) 0)) (228 'name' (CHARACTER 1 0 0 CHARACTER ((CONSTANT (INTEGER 4 0 0
INTEGER ()) 0 '101'))) () (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0) UNKNOWN-ACCESS (CONSTANT (CHARACTER 1 0 0 CHARACTER
(())) 0 101
'                                                                                                     '))
(229 'diagonal' (LOGICAL 4 0 0 LOGICAL ()) () (UNKNOWN-FL UNKNOWN-INTENT
UNKNOWN-PROC UNKNOWN UNKNOWN 0 0) UNKNOWN-ACCESS (CONSTANT (LOGICAL 4 0
0 LOGICAL ()) 0 0)) (230 'equal_diagonal_blocks' (LOGICAL 4 0 0 LOGICAL
()) () (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0)
UNKNOWN-ACCESS (CONSTANT (LOGICAL 4 0 0 LOGICAL ()) 0 0)) (231 'ksp' (
INTEGER 8 0 0 INTEGER ()) () (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0 POINTER) UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN
()) 0))) PUBLIC (() () () ()) () 0 0 61921235)
129 'constraints_type' 'elements' 'constraints_type' 1 ((DERIVED
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 POINTER_COMP) (UNKNOWN 0
0 0 UNKNOWN ()) 0 0 () () 0 ((232 'type' (INTEGER 4 0 0 INTEGER ()) () (
UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0)
UNKNOWN-ACCESS ()) (233 'dim' (INTEGER 4 0 0 INTEGER ()) () (UNKNOWN-FL
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0) UNKNOWN-ACCESS ()) (
234 'degree' (INTEGER 4 0 0 INTEGER ()) () (UNKNOWN-FL UNKNOWN-INTENT
UNKNOWN-PROC UNKNOWN UNKNOWN 0 0) UNKNOWN-ACCESS ()) (235 'loc' (
INTEGER 4 0 0 INTEGER ()) () (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0) UNKNOWN-ACCESS ()) (236 'n_constraints' (INTEGER 4
0 0 INTEGER ()) () (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN
UNKNOWN 0 0) UNKNOWN-ACCESS ()) (237 'orthogonal' (REAL 8 0 0 REAL ()) (
3 0 DEFERRED () () () () () ()) (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0 DIMENSION POINTER) UNKNOWN-ACCESS (NULL (UNKNOWN 0 0
0 UNKNOWN ()) 0))) PUBLIC (() () () ()) () 0 0 57548971)
108 'csr_matrix' 'sparse_tools' 'csr_matrix' 1 ((DERIVED UNKNOWN-INTENT
UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 POINTER_COMP) (UNKNOWN 0 0 0 UNKNOWN ())
0 0 () () 0 ((238 'sparsity' (DERIVED 110 0 0 DERIVED ()) () (
UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0)
UNKNOWN-ACCESS (STRUCTURE (DERIVED 110 0 0 DERIVED ()) 0 (((NULL (
UNKNOWN 0 0 0 UNKNOWN ()) 0) ()) ((NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0) ())
((NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0) ()) (() ()) ((NULL (UNKNOWN 0 0 0
UNKNOWN ()) 0) ()) ((NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0) ()) ((NULL (
UNKNOWN 0 0 0 UNKNOWN ()) 0) ()) ((CONSTANT (CHARACTER 1 0 0 CHARACTER (
())) 0 101
'                                                                                                     ')
()) ((CONSTANT (LOGICAL 4 0 0 LOGICAL ()) 0 0) ()) ((CONSTANT (LOGICAL 4
0 0 LOGICAL ()) 0 0) ())) ())) (239 'val' (REAL 8 0 0 REAL ()) (1 0
DEFERRED () ()) (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN
0 0 DIMENSION POINTER) UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0))
(240 'ival' (INTEGER 4 0 0 INTEGER ()) (1 0 DEFERRED () ()) (UNKNOWN-FL
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 DIMENSION POINTER)
UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0)) (241 'clone' (
LOGICAL 4 0 0 LOGICAL ()) () (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0) UNKNOWN-ACCESS (CONSTANT (LOGICAL 4 0 0 LOGICAL ())
0 0)) (242 'external_val' (LOGICAL 4 0 0 LOGICAL ()) () (UNKNOWN-FL
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0) UNKNOWN-ACCESS (
CONSTANT (LOGICAL 4 0 0 LOGICAL ()) 0 0)) (243 'inactive' (DERIVED 244 0
0 DERIVED ()) () (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN
0 0 POINTER) UNKNOWN-ACCESS ()) (245 'ksp' (INTEGER 8 0 0 INTEGER ()) ()
(UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 POINTER)
UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0)) (246 'refcount' (
DERIVED 55 0 0 DERIVED ()) () (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0 POINTER) UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN
()) 0)) (247 'name' (CHARACTER 1 0 0 CHARACTER ((CONSTANT (INTEGER 4 0 0
INTEGER ()) 0 '101'))) () (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0) UNKNOWN-ACCESS (CONSTANT (CHARACTER 1 0 0 CHARACTER
(())) 0 101
'                                                                                                     ')))
PUBLIC (() () () ()) () 0 0 66371681)
110 'csr_sparsity' 'sparse_tools' 'csr_sparsity' 1 ((DERIVED
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 POINTER_COMP) (UNKNOWN 0
0 0 UNKNOWN ()) 0 0 () () 0 ((248 'findrm' (INTEGER 4 0 0 INTEGER ()) (
1 0 DEFERRED () ()) (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN
UNKNOWN 0 0 DIMENSION POINTER) UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0
UNKNOWN ()) 0)) (249 'centrm' (INTEGER 4 0 0 INTEGER ()) (1 0 DEFERRED ()
()) (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DIMENSION POINTER) UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0)) (
250 'colm' (INTEGER 4 0 0 INTEGER ()) (1 0 DEFERRED () ()) (UNKNOWN-FL
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 DIMENSION POINTER)
UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0)) (251 'columns' (
INTEGER 4 0 0 INTEGER ()) () (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0) UNKNOWN-ACCESS ()) (252 'row_halo' (DERIVED 152 0 0
DERIVED ()) () (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0
0 POINTER) UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0)) (253
'column_halo' (DERIVED 152 0 0 DERIVED ()) () (UNKNOWN-FL UNKNOWN-INTENT
UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 POINTER) UNKNOWN-ACCESS (NULL (UNKNOWN
0 0 0 UNKNOWN ()) 0)) (254 'refcount' (DERIVED 55 0 0 DERIVED ()) () (
UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 POINTER)
UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0)) (255 'name' (
CHARACTER 1 0 0 CHARACTER ((CONSTANT (INTEGER 4 0 0 INTEGER ()) 0 '101')))
() (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0)
UNKNOWN-ACCESS (CONSTANT (CHARACTER 1 0 0 CHARACTER (())) 0 101
'                                                                                                     '))
(256 'wrapped' (LOGICAL 4 0 0 LOGICAL ()) () (UNKNOWN-FL UNKNOWN-INTENT
UNKNOWN-PROC UNKNOWN UNKNOWN 0 0) UNKNOWN-ACCESS (CONSTANT (LOGICAL 4 0
0 LOGICAL ()) 0 0)) (257 'sorted_rows' (LOGICAL 4 0 0 LOGICAL ()) () (
UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0)
UNKNOWN-ACCESS (CONSTANT (LOGICAL 4 0 0 LOGICAL ()) 0 0))) PUBLIC (() ()
() ()) () 0 0 78882345)
123 'ele_numbering_type' 'element_numbering' 'ele_numbering_type' 1 ((
DERIVED UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 POINTER_COMP) (
UNKNOWN 0 0 0 UNKNOWN ()) 0 0 () () 0 ((258 'faces' (INTEGER 4 0 0
INTEGER ()) () (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0
0) UNKNOWN-ACCESS ()) (259 'vertices' (INTEGER 4 0 0 INTEGER ()) () (
UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0)
UNKNOWN-ACCESS ()) (260 'edges' (INTEGER 4 0 0 INTEGER ()) () (
UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0)
UNKNOWN-ACCESS ()) (261 'boundaries' (INTEGER 4 0 0 INTEGER ()) () (
UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0)
UNKNOWN-ACCESS ()) (262 'degree' (INTEGER 4 0 0 INTEGER ()) () (
UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0)
UNKNOWN-ACCESS ()) (263 'dimension' (INTEGER 4 0 0 INTEGER ()) () (
UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0)
UNKNOWN-ACCESS ()) (264 'nodes' (INTEGER 4 0 0 INTEGER ()) () (
UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0)
UNKNOWN-ACCESS ()) (265 'type' (INTEGER 4 0 0 INTEGER ()) () (
UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0)
UNKNOWN-ACCESS (CONSTANT (INTEGER 4 0 0 INTEGER ()) 0 '1')) (266 'family'
(INTEGER 4 0 0 INTEGER ()) () (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0) UNKNOWN-ACCESS ()) (267 'count2number' (INTEGER 4 0
0 INTEGER ()) (3 0 DEFERRED () () () () () ()) (UNKNOWN-FL
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 DIMENSION POINTER)
UNKNOWN-ACCESS ()) (268 'number2count' (INTEGER 4 0 0 INTEGER ()) (2 0
DEFERRED () () () ()) (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN
UNKNOWN 0 0 DIMENSION POINTER) UNKNOWN-ACCESS ()) (269 'boundary_coord'
(INTEGER 4 0 0 INTEGER ()) (1 0 DEFERRED () ()) (UNKNOWN-FL
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 DIMENSION POINTER)
UNKNOWN-ACCESS ()) (270 'boundary_val' (INTEGER 4 0 0 INTEGER ()) (1 0
DEFERRED () ()) (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN
0 0 DIMENSION POINTER) UNKNOWN-ACCESS ())) PUBLIC (() () () ()) () 0 0
96431082)
194 'integer_vector' 'futils' 'integer_vector' 1 ((DERIVED
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 POINTER_COMP) (UNKNOWN 0
0 0 UNKNOWN ()) 0 0 () () 0 ((271 'ptr' (INTEGER 4 0 0 INTEGER ()) (1 0
DEFERRED () ()) (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN
0 0 DIMENSION POINTER) UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0)))
PUBLIC (() () () ()) () 0 0 9661976)
244 'logical_array_ptr' 'sparse_tools' 'logical_array_ptr' 1 ((DERIVED
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 POINTER_COMP) (UNKNOWN 0
0 0 UNKNOWN ()) 0 0 () () 0 ((272 'ptr' (LOGICAL 4 0 0 LOGICAL ()) (1 0
DEFERRED () ()) (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN
0 0 DIMENSION POINTER) UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0)))
PUBLIC (() () () ()) () 0 0 30974511)
143 'mesh_faces' 'fields_data_types' 'mesh_faces' 1 ((DERIVED
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 POINTER_COMP) (UNKNOWN 0
0 0 UNKNOWN ()) 0 0 () () 0 ((273 'shape' (DERIVED 52 0 0 DERIVED ()) ()
(UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 POINTER)
UNKNOWN-ACCESS ()) (274 'face_list' (DERIVED 108 0 0 DERIVED ()) () (
UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0)
UNKNOWN-ACCESS (STRUCTURE (DERIVED 108 0 0 DERIVED ()) 0 (((STRUCTURE (
DERIVED 110 0 0 DERIVED ()) 0 (((NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0) ())
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
())) ())) (275 'face_lno' (INTEGER 4 0 0 INTEGER ()) (1 0 DEFERRED () ())
(UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 DIMENSION
POINTER) UNKNOWN-ACCESS ()) (276 'surface_mesh' (DERIVED 51 0 0 DERIVED
()) () (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0)
UNKNOWN-ACCESS (STRUCTURE (DERIVED 51 0 0 DERIVED ()) 0 ((() ()) ((
CONSTANT (LOGICAL 4 0 0 LOGICAL ()) 0 1) ()) ((STRUCTURE (DERIVED 52 0 0
DERIVED ()) 0 ((() ()) (() ()) (() ()) (() ()) ((NULL (UNKNOWN 0 0 0
UNKNOWN ()) 0) ()) ((NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0) ()) ((NULL (
UNKNOWN 0 0 0 UNKNOWN ()) 0) ()) ((NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0) ())
((NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0) ()) ((NULL (UNKNOWN 0 0 0 UNKNOWN ())
0) ()) ((NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0) ()) ((STRUCTURE (DERIVED 53
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
0) ())) ())) (277 'surface_node_list' (INTEGER 4 0 0 INTEGER ()) (1 0
DEFERRED () ()) (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN
0 0 DIMENSION POINTER) UNKNOWN-ACCESS ()) (278 'face_element_list' (
INTEGER 4 0 0 INTEGER ()) (1 0 DEFERRED () ()) (UNKNOWN-FL
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 DIMENSION POINTER)
UNKNOWN-ACCESS ()) (279 'boundary_ids' (INTEGER 4 0 0 INTEGER ()) (1 0
DEFERRED () ()) (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN
0 0 DIMENSION POINTER) UNKNOWN-ACCESS ()) (280 'coplanar_ids' (INTEGER 4
0 0 INTEGER ()) (1 0 DEFERRED () ()) (UNKNOWN-FL UNKNOWN-INTENT
UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 DIMENSION POINTER) UNKNOWN-ACCESS (
NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0)) (281 'dg_surface_mesh' (DERIVED 51 0
0 DERIVED ()) () (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN
0 0 POINTER) UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0)) (282
'has_internal_boundaries' (LOGICAL 4 0 0 LOGICAL ()) () (UNKNOWN-FL
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0) UNKNOWN-ACCESS (
CONSTANT (LOGICAL 4 0 0 LOGICAL ()) 0 0))) PUBLIC (() () () ()) () 0 0
11936185)
145 'mesh_subdomain_mesh' 'fields_data_types' 'mesh_subdomain_mesh' 1 (
(DERIVED UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 POINTER_COMP) (
UNKNOWN 0 0 0 UNKNOWN ()) 0 0 () () 0 ((283 'element_list' (INTEGER 4 0
0 INTEGER ()) (1 0 DEFERRED () ()) (UNKNOWN-FL UNKNOWN-INTENT
UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 DIMENSION POINTER) UNKNOWN-ACCESS ()) (
284 'node_list' (INTEGER 4 0 0 INTEGER ()) (1 0 DEFERRED () ()) (
UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 DIMENSION
POINTER) UNKNOWN-ACCESS ())) PUBLIC (() () () ()) () 0 0 85302181)
220 'real_vector' 'futils' 'real_vector' 1 ((DERIVED UNKNOWN-INTENT
UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 POINTER_COMP) (UNKNOWN 0 0 0 UNKNOWN ())
0 0 () () 0 ((285 'ptr' (REAL 8 0 0 REAL ()) (1 0 DEFERRED () ()) (
UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 DIMENSION
POINTER) UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0))) PUBLIC (()
() () ()) () 0 0 72870256)
176 'scalar_boundary_condition' 'fields_data_types'
'scalar_boundary_condition' 1 ((DERIVED UNKNOWN-INTENT UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0 POINTER_COMP) (UNKNOWN 0 0 0 UNKNOWN ()) 0 0 () () 0
((286 'name' (CHARACTER 1 0 0 CHARACTER ((CONSTANT (INTEGER 4 0 0
INTEGER ()) 0 '101'))) () (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0) UNKNOWN-ACCESS ()) (287 'type' (CHARACTER 1 0 0
CHARACTER ((CONSTANT (INTEGER 4 0 0 INTEGER ()) 0 '101'))) () (
UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0)
UNKNOWN-ACCESS (CONSTANT (CHARACTER 1 0 0 CHARACTER (())) 0 101
'                                                                                                     '))
(288 'surface_element_list' (INTEGER 4 0 0 INTEGER ()) (1 0 DEFERRED ()
()) (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DIMENSION POINTER) UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0)) (
289 'surface_node_list' (INTEGER 4 0 0 INTEGER ()) (1 0 DEFERRED () ())
(UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 DIMENSION
POINTER) UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0)) (290
'surface_mesh' (DERIVED 51 0 0 DERIVED ()) () (UNKNOWN-FL UNKNOWN-INTENT
UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 POINTER) UNKNOWN-ACCESS ()) (291
'surface_fields' (DERIVED 39 0 0 DERIVED ()) (1 0 DEFERRED () ()) (
UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 DIMENSION
POINTER) UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0)) (292
'option_path' (CHARACTER 1 0 0 CHARACTER ((CONSTANT (INTEGER 4 0 0
INTEGER ()) 0 '8192'))) () (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0) UNKNOWN-ACCESS (CONSTANT (CHARACTER 1 0 0 CHARACTER
(())) 0 8192
'/uninitialised_path/                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            ')))
PUBLIC (() () () ()) () 0 0 50740068)
183 'tensor_field' 'fields_data_types' 'tensor_field' 1 ((DERIVED
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 POINTER_COMP) (UNKNOWN 0
0 0 UNKNOWN ()) 0 0 () () 0 ((293 'val' (REAL 8 0 0 REAL ()) (3 0
DEFERRED () () () () () ()) (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0 DIMENSION POINTER) UNKNOWN-ACCESS ()) (294 'wrapped'
(LOGICAL 4 0 0 LOGICAL ()) () (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0) UNKNOWN-ACCESS (CONSTANT (LOGICAL 4 0 0 LOGICAL ())
0 1)) (295 'field_type' (INTEGER 4 0 0 INTEGER ()) () (UNKNOWN-FL
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0) UNKNOWN-ACCESS (
CONSTANT (INTEGER 4 0 0 INTEGER ()) 0 '0')) (296 'name' (CHARACTER 1 0 0
CHARACTER ((CONSTANT (INTEGER 4 0 0 INTEGER ()) 0 '101'))) () (
UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0)
UNKNOWN-ACCESS ()) (297 'dim' (INTEGER 4 0 0 INTEGER ()) (1 0 EXPLICIT (
CONSTANT (INTEGER 4 0 0 INTEGER ()) 0 '1') (CONSTANT (INTEGER 4 0 0
INTEGER ()) 0 '2')) (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN
UNKNOWN 0 0 DIMENSION) UNKNOWN-ACCESS ()) (298 'option_path' (CHARACTER
1 0 0 CHARACTER ((CONSTANT (INTEGER 4 0 0 INTEGER ()) 0 '8192'))) () (
UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0)
UNKNOWN-ACCESS (CONSTANT (CHARACTER 1 0 0 CHARACTER (())) 0 8192
'/uninitialised_path/                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            '))
(299 'mesh' (DERIVED 51 0 0 DERIVED ()) () (UNKNOWN-FL UNKNOWN-INTENT
UNKNOWN-PROC UNKNOWN UNKNOWN 0 0) UNKNOWN-ACCESS (STRUCTURE (DERIVED 51
0 0 DERIVED ()) 0 ((() ()) ((CONSTANT (LOGICAL 4 0 0 LOGICAL ()) 0 1) ())
((STRUCTURE (DERIVED 52 0 0 DERIVED ()) 0 ((() ()) (() ()) (() ()) (() ())
((NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0) ()) ((NULL (UNKNOWN 0 0 0 UNKNOWN ())
0) ()) ((NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0) ()) ((NULL (UNKNOWN 0 0 0
UNKNOWN ()) 0) ()) ((NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0) ()) ((NULL (
UNKNOWN 0 0 0 UNKNOWN ()) 0) ()) ((NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0) ())
((STRUCTURE (DERIVED 53 0 0 DERIVED ()) 0 ((() ()) (() ()) (() ()) (() ())
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
UNKNOWN ()) 0) ()) ((NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0) ()) ((NULL (
UNKNOWN 0 0 0 UNKNOWN ()) 0) ()) ((CONSTANT (LOGICAL 4 0 0 LOGICAL ()) 0
0) ())) ())) (300 'refcount' (DERIVED 55 0 0 DERIVED ()) () (UNKNOWN-FL
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 POINTER) UNKNOWN-ACCESS
(NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0)) (301 'aliased' (LOGICAL 4 0 0
LOGICAL ()) () (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0
0) UNKNOWN-ACCESS (CONSTANT (LOGICAL 4 0 0 LOGICAL ()) 0 0))) PUBLIC (()
() () ()) () 0 0 25110185)
185 'vector_boundary_condition' 'fields_data_types'
'vector_boundary_condition' 1 ((DERIVED UNKNOWN-INTENT UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0 POINTER_COMP) (UNKNOWN 0 0 0 UNKNOWN ()) 0 0 () () 0
((302 'name' (CHARACTER 1 0 0 CHARACTER ((CONSTANT (INTEGER 4 0 0
INTEGER ()) 0 '101'))) () (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0) UNKNOWN-ACCESS ()) (303 'type' (CHARACTER 1 0 0
CHARACTER ((CONSTANT (INTEGER 4 0 0 INTEGER ()) 0 '101'))) () (
UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0)
UNKNOWN-ACCESS (CONSTANT (CHARACTER 1 0 0 CHARACTER (())) 0 101
'                                                                                                     '))
(304 'applies' (LOGICAL 4 0 0 LOGICAL ()) (1 0 EXPLICIT (CONSTANT (
INTEGER 4 0 0 INTEGER ()) 0 '1') (CONSTANT (INTEGER 4 0 0 INTEGER ()) 0
'3')) (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DIMENSION) UNKNOWN-ACCESS ()) (305 'surface_element_list' (INTEGER 4 0 0
INTEGER ()) (1 0 DEFERRED () ()) (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0 DIMENSION POINTER) UNKNOWN-ACCESS (NULL (UNKNOWN 0 0
0 UNKNOWN ()) 0)) (306 'surface_node_list' (INTEGER 4 0 0 INTEGER ()) (
1 0 DEFERRED () ()) (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN
UNKNOWN 0 0 DIMENSION POINTER) UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0
UNKNOWN ()) 0)) (307 'surface_mesh' (DERIVED 51 0 0 DERIVED ()) () (
UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 POINTER)
UNKNOWN-ACCESS ()) (308 'surface_fields' (DERIVED 40 0 0 DERIVED ()) (1
0 DEFERRED () ()) (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN
UNKNOWN 0 0 DIMENSION POINTER) UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0
UNKNOWN ()) 0)) (309 'scalar_surface_fields' (DERIVED 39 0 0 DERIVED ())
(1 0 DEFERRED () ()) (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN
UNKNOWN 0 0 DIMENSION POINTER) UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0
UNKNOWN ()) 0)) (310 'option_path' (CHARACTER 1 0 0 CHARACTER ((
CONSTANT (INTEGER 4 0 0 INTEGER ()) 0 '8192'))) () (UNKNOWN-FL
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0) UNKNOWN-ACCESS (
CONSTANT (CHARACTER 1 0 0 CHARACTER (())) 0 8192
'/uninitialised_path/                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            ')))
PUBLIC (() () () ()) () 0 0 42205165)
75 'picker_ptr' 'picker_data_types' 'picker_ptr' 1 ((DERIVED
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 POINTER_COMP) (UNKNOWN 0
0 0 UNKNOWN ()) 0 0 () () 0 ((311 'ptr' (DERIVED 312 0 0 DERIVED ()) ()
(UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 POINTER)
UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0))) PUBLIC (() () () ())
() 0 0 71120007)
155 'integer_set_vector' 'integer_set_module' 'integer_set_vector' 1 ((
DERIVED UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 POINTER_COMP) (
UNKNOWN 0 0 0 UNKNOWN ()) 0 0 () () 0 ((313 'sets' (DERIVED 314 0 0
DERIVED ()) (1 0 DEFERRED () ()) (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0 DIMENSION POINTER) UNKNOWN-ACCESS ())) PUBLIC (() ()
() ()) () 0 0 40688950)
32 'state' '' 'state' 31 ((VARIABLE INOUT UNKNOWN-PROC UNKNOWN UNKNOWN 0
0 DUMMY) (DERIVED 41 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
208 'petsc_numbering_type' 'petsc_tools' 'petsc_numbering_type' 1 ((
DERIVED UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 POINTER_COMP) (
UNKNOWN 0 0 0 UNKNOWN ()) 0 0 () () 0 ((315 'halo' (DERIVED 152 0 0
DERIVED ()) () (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0
0 POINTER) UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0)) (316
'nprivatenodes' (INTEGER 4 0 0 INTEGER ()) () (UNKNOWN-FL UNKNOWN-INTENT
UNKNOWN-PROC UNKNOWN UNKNOWN 0 0) UNKNOWN-ACCESS ()) (317
'universal_length' (INTEGER 4 0 0 INTEGER ()) () (UNKNOWN-FL
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0) UNKNOWN-ACCESS ()) (
318 'offset' (INTEGER 4 0 0 INTEGER ()) () (UNKNOWN-FL UNKNOWN-INTENT
UNKNOWN-PROC UNKNOWN UNKNOWN 0 0) UNKNOWN-ACCESS ()) (319 'gnn2unn' (
INTEGER 4 0 0 INTEGER ()) (2 0 DEFERRED () () () ()) (UNKNOWN-FL
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 DIMENSION POINTER)
UNKNOWN-ACCESS ()) (320 'ghost_nodes' (INTEGER 4 0 0 INTEGER ()) (1 0
DEFERRED () ()) (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN
0 0 DIMENSION POINTER) UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0))
(321 'ghost2unn' (INTEGER 4 0 0 INTEGER ()) (2 0 DEFERRED () () () ()) (
UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 DIMENSION
POINTER) UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0)) (322
'refcount' (DERIVED 55 0 0 DERIVED ()) () (UNKNOWN-FL UNKNOWN-INTENT
UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 POINTER) UNKNOWN-ACCESS (NULL (UNKNOWN
0 0 0 UNKNOWN ()) 0)) (323 'name' (CHARACTER 1 0 0 CHARACTER ((CONSTANT
(INTEGER 4 0 0 INTEGER ()) 0 '101'))) () (UNKNOWN-FL UNKNOWN-INTENT
UNKNOWN-PROC UNKNOWN UNKNOWN 0 0) UNKNOWN-ACCESS (CONSTANT (CHARACTER 1
0 0 CHARACTER (())) 0 101
'                                                                                                     ')))
PUBLIC (() () () ()) () 0 0 16529124)
314 'integer_set' 'integer_set_module' 'integer_set' 1 ((DERIVED
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0) (UNKNOWN 0 0 0 UNKNOWN
()) 0 0 () () 0 ((324 'address' (DERIVED 325 0 0 DERIVED ()) () (
UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0)
UNKNOWN-ACCESS ())) PUBLIC (() () () ()) () 0 0 67217964)
312 'picker_type' 'picker_data_types' 'picker_type' 1 ((DERIVED
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 POINTER_COMP) (UNKNOWN 0
0 0 UNKNOWN ()) 0 0 () () 0 ((326 'name' (CHARACTER 1 0 0 CHARACTER ((
CONSTANT (INTEGER 4 0 0 INTEGER ()) 0 '101'))) () (UNKNOWN-FL
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0) UNKNOWN-ACCESS ()) (
327 'refcount' (DERIVED 55 0 0 DERIVED ()) () (UNKNOWN-FL UNKNOWN-INTENT
UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 POINTER) UNKNOWN-ACCESS (NULL (UNKNOWN
0 0 0 UNKNOWN ()) 0)) (328 'picker_id' (INTEGER 4 0 0 INTEGER ()) () (
UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0)
UNKNOWN-ACCESS (CONSTANT (INTEGER 4 0 0 INTEGER ()) 0 '0')) (329
'last_mesh_movement' (INTEGER 4 0 0 INTEGER ()) () (UNKNOWN-FL
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0) UNKNOWN-ACCESS (
CONSTANT (INTEGER 4 0 0 INTEGER ()) 0 '0'))) PUBLIC (() () () ()) () 0 0
8821665)
325 'c_ptr' '__iso_c_binding' '' 1 ((DERIVED UNKNOWN-INTENT UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0 IS_BIND_C IS_C_INTEROP IS_ISO_C) (DERIVED 325 1 1
UNKNOWN ()) 0 0 () () 0 ((330 '__c_ptr_c_address' (INTEGER 8 1 0 INTEGER
()) () (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0)
UNKNOWN-ACCESS ())) UNKNOWN-ACCESS () () 2 39 0)
)

('limit_fpn' 0 19 'limit_slope_dg' 0 23 'limit_vb' 0 30 'limiter_cockburn'
0 34 'limiter_fpn' 0 35 'limiter_hermite_weno' 0 36 'limiter_minimal' 0
37 'limiter_vb' 0 38)
