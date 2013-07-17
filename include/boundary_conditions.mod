GFORTRAN module version '6' created from Boundary_Conditions.F90 on Mon Jun 10 19:47:29 2013
MD5:740e79d3532e6508bf46d6c137fa88cc -- If you edit this, you'll get what you deserve.

(() () () () () () () () () () () () () () () () () () () () () () ()
() () () ())

()

(('add_boundary_condition' 'boundary_conditions' 2 3) (
'add_boundary_condition_surface_elements' 'boundary_conditions' 4 5) (
'apply_dirichlet_conditions' 'boundary_conditions' 6 7 8 9 10 11 12 13)
('extract_scalar_surface_field' 'boundary_conditions' 14 15) (
'extract_surface_field' 'boundary_conditions' 16 17 18 19) (
'get_boundary_condition_count' 'boundary_conditions' 20 21) (
'get_entire_boundary_condition' 'boundary_conditions' 22 23) (
'get_boundary_condition_nodes' 'boundary_conditions' 24 25) (
'get_boundary_condition' 'boundary_conditions' 26 27 28 29) (
'has_boundary_condition_name' 'boundary_conditions' 30 31) (
'has_surface_field' 'boundary_conditions' 32 33) ('insert_surface_field'
'boundary_conditions' 34 35 36 37 38 39) ('has_boundary_condition'
'boundary_conditions' 40 41) ('remove_boundary_condition'
'boundary_conditions' 42 43) ('set_reference_node' 'boundary_conditions'
44 45) ('set_dirichlet_consistent' 'boundary_conditions' 46 47 48))

(('mpi_fortran_argv_null' 49 0 0 'mpi_fortran_argv_null') (
'mpi_fortran_errcodes_ignore' 50 0 0 'mpi_fortran_errcodes_ignore') (
'mpi_fortran_in_place' 51 0 0 'mpi_fortran_in_place') (
'mpi_fortran_status_ignore' 52 0 0 'mpi_fortran_status_ignore') (
'mpi_fortran_statuses_ignore' 53 0 0 'mpi_fortran_statuses_ignore') (
'petscfortran1' 54 0 0 'petscfortran1') ('petscfortran10' 55 0 0
'petscfortran10') ('petscfortran2' 56 0 0 'petscfortran2') (
'petscfortran3' 57 0 0 'petscfortran3') ('mpi_fortran_argvs_null' 58 0 0
'mpi_fortran_argvs_null') ('mpi_fortran_bottom' 59 0 0
'mpi_fortran_bottom') ('petscfortran4' 60 0 0 'petscfortran4') (
'petscfortran5' 61 0 0 'petscfortran5') ('petscfortran6' 62 0 0
'petscfortran6') ('petscfortran7' 63 0 0 'petscfortran7') (
'petscfortran8' 64 0 0 'petscfortran8') ('petscfortran9' 65 0 0
'petscfortran9'))

()

()

(66 'derive_collapsed_bcs' 'boundary_conditions' 'derive_collapsed_bcs'
1 ((PROCEDURE UNKNOWN-INTENT MODULE-PROC DECL UNKNOWN 0 0 SUBROUTINE
ALWAYS_EXPLICIT) (UNKNOWN 0 0 0 UNKNOWN ()) 67 0 (68 69 70) () 0 () () ()
0 0)
71 'get_dg_surface_mesh' 'boundary_conditions' 'get_dg_surface_mesh' 1 (
(PROCEDURE UNKNOWN-INTENT MODULE-PROC DECL UNKNOWN 0 0 POINTER FUNCTION
ALWAYS_EXPLICIT) (DERIVED 72 0 0 DERIVED ()) 73 0 (74) () 71 () () () 0
0)
75 'get_periodic_boundary_condition' 'boundary_conditions'
'get_periodic_boundary_condition' 1 ((PROCEDURE UNKNOWN-INTENT
MODULE-PROC DECL UNKNOWN 0 0 SUBROUTINE ALWAYS_EXPLICIT) (UNKNOWN 0 0 0
UNKNOWN ()) 76 0 (77 78) () 0 () () () 0 0)
79 'has_scalar_surface_field' 'boundary_conditions'
'has_scalar_surface_field' 1 ((PROCEDURE UNKNOWN-INTENT MODULE-PROC DECL
UNKNOWN 0 0 FUNCTION) (LOGICAL 4 0 0 LOGICAL ()) 80 0 (81 82 83) () 79 ()
() () 0 0)
48 'set_dirichlet_consistent' 'boundary_conditions'
'set_dirichlet_consistent' 1 ((PROCEDURE UNKNOWN-INTENT MODULE-PROC DECL
UNKNOWN 0 0 SUBROUTINE GENERIC ALWAYS_EXPLICIT) (UNKNOWN 0 0 0 UNKNOWN ())
84 0 (85) () 0 () () () 0 0)
72 'mesh_type' 'fields_data_types' 'mesh_type' 1 ((DERIVED
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 POINTER_COMP) (UNKNOWN 0
0 0 UNKNOWN ()) 0 0 () () 0 ((86 'ndglno' (INTEGER 4 0 0 INTEGER ()) (1
0 DEFERRED () ()) (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN
UNKNOWN 0 0 DIMENSION POINTER) UNKNOWN-ACCESS ()) (87 'wrapped' (
LOGICAL 4 0 0 LOGICAL ()) () (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0) UNKNOWN-ACCESS (CONSTANT (LOGICAL 4 0 0 LOGICAL ())
0 1)) (88 'shape' (DERIVED 89 0 0 DERIVED ()) () (UNKNOWN-FL
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0) UNKNOWN-ACCESS (
STRUCTURE (DERIVED 89 0 0 DERIVED ()) 0 ((() ()) (() ()) (() ()) (() ())
((NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0) ()) ((NULL (UNKNOWN 0 0 0 UNKNOWN ())
0) ()) ((NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0) ()) ((NULL (UNKNOWN 0 0 0
UNKNOWN ()) 0) ()) ((NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0) ()) ((NULL (
UNKNOWN 0 0 0 UNKNOWN ()) 0) ()) ((NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0) ())
((STRUCTURE (DERIVED 90 0 0 DERIVED ()) 0 ((() ()) (() ()) (() ()) (() ())
((NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0) ()) ((NULL (UNKNOWN 0 0 0 UNKNOWN ())
0) ()) (() ()) ((NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0) ()) (() ())) ()) ())
((NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0) ()) ((NULL (UNKNOWN 0 0 0 UNKNOWN ())
0) ()) ((NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0) ()) ((NULL (UNKNOWN 0 0 0
UNKNOWN ()) 0) ()) (() ())) ())) (91 'elements' (INTEGER 4 0 0 INTEGER ())
() (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0)
UNKNOWN-ACCESS ()) (92 'nodes' (INTEGER 4 0 0 INTEGER ()) () (
UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0)
UNKNOWN-ACCESS ()) (93 'name' (CHARACTER 1 0 0 CHARACTER ((CONSTANT (
INTEGER 4 0 0 INTEGER ()) 0 '101'))) () (UNKNOWN-FL UNKNOWN-INTENT
UNKNOWN-PROC UNKNOWN UNKNOWN 0 0) UNKNOWN-ACCESS ()) (94 'option_path' (
CHARACTER 1 0 0 CHARACTER ((CONSTANT (INTEGER 4 0 0 INTEGER ()) 0 '8192')))
() (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0)
UNKNOWN-ACCESS (CONSTANT (CHARACTER 1 0 0 CHARACTER (())) 0 8192
'/uninitialised_path/                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            '))
(95 'continuity' (INTEGER 4 0 0 INTEGER ()) () (UNKNOWN-FL
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0) UNKNOWN-ACCESS (
CONSTANT (INTEGER 4 0 0 INTEGER ()) 0 '0')) (96 'refcount' (DERIVED 97 0
0 DERIVED ()) () (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN
0 0 POINTER) UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0)) (98
'faces' (DERIVED 99 0 0 DERIVED ()) () (UNKNOWN-FL UNKNOWN-INTENT
UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 POINTER) UNKNOWN-ACCESS (NULL (UNKNOWN
0 0 0 UNKNOWN ()) 0)) (100 'subdomain_mesh' (DERIVED 101 0 0 DERIVED ())
() (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 POINTER)
UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0)) (102 'adj_lists' (
DERIVED 103 0 0 DERIVED ()) () (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0 POINTER) UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN
()) 0)) (104 'columns' (INTEGER 4 0 0 INTEGER ()) (1 0 DEFERRED () ()) (
UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 DIMENSION
POINTER) UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0)) (105
'element_columns' (INTEGER 4 0 0 INTEGER ()) (1 0 DEFERRED () ()) (
UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 DIMENSION
POINTER) UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0)) (106
'region_ids' (INTEGER 4 0 0 INTEGER ()) (1 0 DEFERRED () ()) (
UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 DIMENSION
POINTER) UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0)) (107 'halos'
(DERIVED 108 0 0 DERIVED ()) (1 0 DEFERRED () ()) (UNKNOWN-FL
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 DIMENSION POINTER)
UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0)) (109 'element_halos'
(DERIVED 108 0 0 DERIVED ()) (1 0 DEFERRED () ()) (UNKNOWN-FL
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 DIMENSION POINTER)
UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0)) (110 'colourings' (
DERIVED 111 0 0 DERIVED ()) (1 0 DEFERRED () ()) (UNKNOWN-FL
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 DIMENSION POINTER)
UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0)) (112 'periodic' (
LOGICAL 4 0 0 LOGICAL ()) () (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0) UNKNOWN-ACCESS (CONSTANT (LOGICAL 4 0 0 LOGICAL ())
0 0))) PUBLIC (() () () ()) () 0 0 2572791)
57 'petsc_null' 'petscsysdef' 'petsc_null' 1 ((VARIABLE UNKNOWN-INTENT
UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 IN_COMMON) (INTEGER 4 0 0 INTEGER ()) 0
0 () () 0 () () () 0 0)
49 'mpi_argv_null' 'mpi_interfaces' 'mpi_argv_null' 1 ((VARIABLE
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 DIMENSION IN_COMMON) (
CHARACTER 1 0 0 CHARACTER ((CONSTANT (INTEGER 4 0 0 INTEGER ()) 0 '1')))
0 0 () (1 0 EXPLICIT (CONSTANT (INTEGER 4 0 0 INTEGER ()) 0 '1') (
CONSTANT (INTEGER 4 0 0 INTEGER ()) 0 '1')) 0 () () () 0 0)
60 'petsc_null_scalar' 'petscsysdef' 'petsc_null_scalar' 1 ((VARIABLE
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 IN_COMMON) (REAL 8 0 0
REAL ()) 0 0 () () 0 () () () 0 0)
63 'petsc_null_truth' 'petscsysdef' 'petsc_null_truth' 1 ((VARIABLE
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 IN_COMMON) (LOGICAL 4 0
0 LOGICAL ()) 0 0 () () 0 () () () 0 0)
64 'petsc_null_object' 'petscsysdef' 'petsc_null_object' 1 ((VARIABLE
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 IN_COMMON) (INTEGER 8 0
0 INTEGER ()) 0 0 () () 0 () () () 0 0)
58 'mpi_argvs_null' 'mpi_interfaces' 'mpi_argvs_null' 1 ((VARIABLE
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 IN_COMMON) (REAL 8 0 0
REAL ()) 0 0 () () 0 () () () 0 0)
65 'petsc_comm_world' 'petscsysdef' 'petsc_comm_world' 1 ((VARIABLE
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 IN_COMMON) (INTEGER 4 0
0 INTEGER ()) 0 0 () () 0 () () () 0 0)
62 'petsc_null_real' 'petscsysdef' 'petsc_null_real' 1 ((VARIABLE
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 IN_COMMON) (REAL 8 0 0
REAL ()) 0 0 () () 0 () () () 0 0)
59 'mpi_bottom' 'mpi_interfaces' 'mpi_bottom' 1 ((VARIABLE
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 IN_COMMON) (INTEGER 4 0
0 INTEGER ()) 0 0 () () 0 () () () 0 0)
61 'petsc_null_double' 'petscsysdef' 'petsc_null_double' 1 ((VARIABLE
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 IN_COMMON) (REAL 8 0 0
REAL ()) 0 0 () () 0 () () () 0 0)
52 'mpi_status_ignore' 'mpi_interfaces' 'mpi_status_ignore' 1 ((
VARIABLE UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 DIMENSION
IN_COMMON) (INTEGER 4 0 0 INTEGER ()) 0 0 () (1 0 EXPLICIT (CONSTANT (
INTEGER 4 0 0 INTEGER ()) 0 '1') (CONSTANT (INTEGER 4 0 0 INTEGER ()) 0
'5')) 0 () () () 0 0)
51 'mpi_in_place' 'mpi_interfaces' 'mpi_in_place' 1 ((VARIABLE
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 IN_COMMON) (INTEGER 4 0
0 INTEGER ()) 0 0 () () 0 () () () 0 0)
54 'petsc_null_character' 'petscsysdef' 'petsc_null_character' 1 ((
VARIABLE UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN IMPLICIT-SAVE 0 0 IN_COMMON)
(CHARACTER 1 0 0 CHARACTER ((CONSTANT (INTEGER 4 0 0 INTEGER ()) 0 '80')))
0 0 () () 0 () () () 0 0)
53 'mpi_statuses_ignore' 'mpi_interfaces' 'mpi_statuses_ignore' 1 ((
VARIABLE UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 IN_COMMON) (
REAL 8 0 0 REAL ()) 0 0 () () 0 () () () 0 0)
56 'petsc_null_integer' 'petscsysdef' 'petsc_null_integer' 1 ((VARIABLE
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 IN_COMMON) (INTEGER 4 0
0 INTEGER ()) 0 0 () () 0 () () () 0 0)
89 'element_type' 'elements' 'element_type' 1 ((DERIVED UNKNOWN-INTENT
UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 POINTER_COMP) (UNKNOWN 0 0 0 UNKNOWN ())
0 0 () () 0 ((113 'dim' (INTEGER 4 0 0 INTEGER ()) () (UNKNOWN-FL
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0) UNKNOWN-ACCESS ()) (
114 'loc' (INTEGER 4 0 0 INTEGER ()) () (UNKNOWN-FL UNKNOWN-INTENT
UNKNOWN-PROC UNKNOWN UNKNOWN 0 0) UNKNOWN-ACCESS ()) (115 'ngi' (
INTEGER 4 0 0 INTEGER ()) () (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0) UNKNOWN-ACCESS ()) (116 'degree' (INTEGER 4 0 0
INTEGER ()) () (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0
0) UNKNOWN-ACCESS ()) (117 'n' (REAL 8 0 0 REAL ()) (2 0 DEFERRED () ()
() ()) (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DIMENSION POINTER) UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0)) (
118 'dn' (REAL 8 0 0 REAL ()) (3 0 DEFERRED () () () () () ()) (
UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 DIMENSION
POINTER) UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0)) (119 'n_s'
(REAL 8 0 0 REAL ()) (3 0 DEFERRED () () () () () ()) (UNKNOWN-FL
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 DIMENSION POINTER)
UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0)) (120 'dn_s' (REAL 8
0 0 REAL ()) (4 0 DEFERRED () () () () () () () ()) (UNKNOWN-FL
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 DIMENSION POINTER)
UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0)) (121 'spoly' (
DERIVED 122 0 0 DERIVED ()) (2 0 DEFERRED () () () ()) (UNKNOWN-FL
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 DIMENSION POINTER)
UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0)) (123 'dspoly' (
DERIVED 122 0 0 DERIVED ()) (2 0 DEFERRED () () () ()) (UNKNOWN-FL
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 DIMENSION POINTER)
UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0)) (124 'numbering' (
DERIVED 125 0 0 DERIVED ()) () (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0 POINTER) UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN
()) 0)) (126 'quadrature' (DERIVED 90 0 0 DERIVED ()) () (UNKNOWN-FL
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0) UNKNOWN-ACCESS (
STRUCTURE (DERIVED 90 0 0 DERIVED ()) 0 ((() ()) (() ()) (() ()) (() ())
((NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0) ()) ((NULL (UNKNOWN 0 0 0 UNKNOWN ())
0) ()) (() ()) ((NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0) ()) (() ())) ())) (
127 'surface_quadrature' (DERIVED 90 0 0 DERIVED ()) () (UNKNOWN-FL
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 POINTER) UNKNOWN-ACCESS
(NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0)) (128 'superconvergence' (DERIVED
129 0 0 DERIVED ()) () (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN
UNKNOWN 0 0 POINTER) UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0))
(130 'constraints' (DERIVED 131 0 0 DERIVED ()) () (UNKNOWN-FL
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 POINTER) UNKNOWN-ACCESS
(NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0)) (132 'refcount' (DERIVED 97 0 0
DERIVED ()) () (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0
0 POINTER) UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0)) (133 'name'
(CHARACTER 1 0 0 CHARACTER ((CONSTANT (INTEGER 4 0 0 INTEGER ()) 0 '0')))
() (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0)
UNKNOWN-ACCESS ())) PUBLIC (() () () ()) () 0 0 79461029)
3 'add_scalar_boundary_condition' 'boundary_conditions'
'add_scalar_boundary_condition' 1 ((PROCEDURE UNKNOWN-INTENT MODULE-PROC
DECL UNKNOWN 0 0 SUBROUTINE ALWAYS_EXPLICIT) (UNKNOWN 0 0 0 UNKNOWN ())
134 0 (135 136 137 138 139) () 0 () () () 0 0)
2 'add_vector_boundary_condition' 'boundary_conditions'
'add_vector_boundary_condition' 1 ((PROCEDURE UNKNOWN-INTENT MODULE-PROC
DECL UNKNOWN 0 0 SUBROUTINE ALWAYS_EXPLICIT) (UNKNOWN 0 0 0 UNKNOWN ())
140 0 (141 142 143 144 145 146) () 0 () () () 0 0)
43 'remove_scalar_boundary_condition' 'boundary_conditions'
'remove_scalar_boundary_condition' 1 ((PROCEDURE UNKNOWN-INTENT
MODULE-PROC DECL UNKNOWN 0 0 SUBROUTINE) (UNKNOWN 0 0 0 UNKNOWN ()) 147
0 (148 149) () 0 () () () 0 0)
42 'remove_vector_boundary_condition' 'boundary_conditions'
'remove_vector_boundary_condition' 1 ((PROCEDURE UNKNOWN-INTENT
MODULE-PROC DECL UNKNOWN 0 0 SUBROUTINE) (UNKNOWN 0 0 0 UNKNOWN ()) 150
0 (151 152) () 0 () () () 0 0)
5 'add_scalar_boundary_condition_surface_elements' 'boundary_conditions'
'add_scalar_boundary_condition_surface_elements' 1 ((PROCEDURE
UNKNOWN-INTENT MODULE-PROC DECL UNKNOWN 0 0 SUBROUTINE ALWAYS_EXPLICIT)
(UNKNOWN 0 0 0 UNKNOWN ()) 153 0 (154 155 156 157 158) () 0 () () () 0 0)
4 'add_vector_boundary_condition_surface_elements' 'boundary_conditions'
'add_vector_boundary_condition_surface_elements' 1 ((PROCEDURE
UNKNOWN-INTENT MODULE-PROC DECL UNKNOWN 0 0 SUBROUTINE ALWAYS_EXPLICIT)
(UNKNOWN 0 0 0 UNKNOWN ()) 159 0 (160 161 162 163 164 165) () 0 () () ()
0 0)
29 'get_scalar_boundary_condition_by_number' 'boundary_conditions'
'get_scalar_boundary_condition_by_number' 1 ((PROCEDURE UNKNOWN-INTENT
MODULE-PROC DECL UNKNOWN 0 0 SUBROUTINE ALWAYS_EXPLICIT) (UNKNOWN 0 0 0
UNKNOWN ()) 166 0 (167 168 169 170 171 172 173 174) () 0 () () () 0 0)
28 'get_vector_boundary_condition_by_number' 'boundary_conditions'
'get_vector_boundary_condition_by_number' 1 ((PROCEDURE UNKNOWN-INTENT
MODULE-PROC DECL UNKNOWN 0 0 SUBROUTINE ALWAYS_EXPLICIT) (UNKNOWN 0 0 0
UNKNOWN ()) 175 0 (176 177 178 179 180 181 182 183 184) () 0 () () () 0
0)
27 'get_scalar_boundary_condition_by_name' 'boundary_conditions'
'get_scalar_boundary_condition_by_name' 1 ((PROCEDURE UNKNOWN-INTENT
MODULE-PROC DECL UNKNOWN 0 0 SUBROUTINE ALWAYS_EXPLICIT) (UNKNOWN 0 0 0
UNKNOWN ()) 185 0 (186 187 188 189 190 191 192) () 0 () () () 0 0)
26 'get_vector_boundary_condition_by_name' 'boundary_conditions'
'get_vector_boundary_condition_by_name' 1 ((PROCEDURE UNKNOWN-INTENT
MODULE-PROC DECL UNKNOWN 0 0 SUBROUTINE ALWAYS_EXPLICIT) (UNKNOWN 0 0 0
UNKNOWN ()) 193 0 (194 195 196 197 198 199 200 201) () 0 () () () 0 0)
39 'insert_scalar_surface_field' 'boundary_conditions'
'insert_scalar_surface_field' 1 ((PROCEDURE UNKNOWN-INTENT MODULE-PROC
DECL UNKNOWN 0 0 SUBROUTINE) (UNKNOWN 0 0 0 UNKNOWN ()) 202 0 (203 204
205) () 0 () () () 0 0)
38 'insert_vector_surface_field' 'boundary_conditions'
'insert_vector_surface_field' 1 ((PROCEDURE UNKNOWN-INTENT MODULE-PROC
DECL UNKNOWN 0 0 SUBROUTINE) (UNKNOWN 0 0 0 UNKNOWN ()) 206 0 (207 208
209) () 0 () () () 0 0)
37 'insert_scalar_surface_field_by_name' 'boundary_conditions'
'insert_scalar_surface_field_by_name' 1 ((PROCEDURE UNKNOWN-INTENT
MODULE-PROC DECL UNKNOWN 0 0 SUBROUTINE) (UNKNOWN 0 0 0 UNKNOWN ()) 210
0 (211 212 213) () 0 () () () 0 0)
36 'insert_vector_surface_field_by_name' 'boundary_conditions'
'insert_vector_surface_field_by_name' 1 ((PROCEDURE UNKNOWN-INTENT
MODULE-PROC DECL UNKNOWN 0 0 SUBROUTINE) (UNKNOWN 0 0 0 UNKNOWN ()) 214
0 (215 216 217) () 0 () () () 0 0)
35 'insert_vector_scalar_surface_field' 'boundary_conditions'
'insert_vector_scalar_surface_field' 1 ((PROCEDURE UNKNOWN-INTENT
MODULE-PROC DECL UNKNOWN 0 0 SUBROUTINE) (UNKNOWN 0 0 0 UNKNOWN ()) 218
0 (219 220 221) () 0 () () () 0 0)
34 'insert_vector_scalar_surface_field_by_name' 'boundary_conditions'
'insert_vector_scalar_surface_field_by_name' 1 ((PROCEDURE
UNKNOWN-INTENT MODULE-PROC DECL UNKNOWN 0 0 SUBROUTINE) (UNKNOWN 0 0 0
UNKNOWN ()) 222 0 (223 224 225) () 0 () () () 0 0)
19 'extract_scalar_surface_field_by_number' 'boundary_conditions'
'extract_scalar_surface_field_by_number' 1 ((PROCEDURE UNKNOWN-INTENT
MODULE-PROC DECL UNKNOWN 0 0 POINTER FUNCTION ALWAYS_EXPLICIT) (DERIVED
226 0 0 DERIVED ()) 227 0 (228 229 230) () 231 () () () 0 0)
18 'extract_vector_surface_field_by_number' 'boundary_conditions'
'extract_vector_surface_field_by_number' 1 ((PROCEDURE UNKNOWN-INTENT
MODULE-PROC DECL UNKNOWN 0 0 POINTER FUNCTION ALWAYS_EXPLICIT) (DERIVED
232 0 0 DERIVED ()) 233 0 (234 235 236) () 237 () () () 0 0)
17 'extract_scalar_surface_field_by_name' 'boundary_conditions'
'extract_scalar_surface_field_by_name' 1 ((PROCEDURE UNKNOWN-INTENT
MODULE-PROC DECL UNKNOWN 0 0 POINTER FUNCTION ALWAYS_EXPLICIT) (DERIVED
226 0 0 DERIVED ()) 238 0 (239 240 241) () 242 () () () 0 0)
16 'extract_vector_surface_field_by_name' 'boundary_conditions'
'extract_vector_surface_field_by_name' 1 ((PROCEDURE UNKNOWN-INTENT
MODULE-PROC DECL UNKNOWN 0 0 POINTER FUNCTION ALWAYS_EXPLICIT) (DERIVED
232 0 0 DERIVED ()) 243 0 (244 245 246) () 247 () () () 0 0)
15 'extract_vector_scalar_surface_field' 'boundary_conditions'
'extract_vector_scalar_surface_field' 1 ((PROCEDURE UNKNOWN-INTENT
MODULE-PROC DECL UNKNOWN 0 0 POINTER FUNCTION ALWAYS_EXPLICIT) (DERIVED
226 0 0 DERIVED ()) 248 0 (249 250 251) () 252 () () () 0 0)
14 'extract_vector_scalar_surface_field_by_name' 'boundary_conditions'
'extract_vector_scalar_surface_field_by_name' 1 ((PROCEDURE
UNKNOWN-INTENT MODULE-PROC DECL UNKNOWN 0 0 POINTER FUNCTION
ALWAYS_EXPLICIT) (DERIVED 226 0 0 DERIVED ()) 253 0 (254 255 256) () 257
() () () 0 0)
33 'has_scalar_surface_field_specific' 'boundary_conditions'
'has_scalar_surface_field_specific' 1 ((PROCEDURE UNKNOWN-INTENT
MODULE-PROC DECL UNKNOWN 0 0 FUNCTION) (LOGICAL 4 0 0 LOGICAL ()) 258 0
(259 260 261) () 33 () () () 0 0)
32 'has_vector_surface_field_specific' 'boundary_conditions'
'has_vector_surface_field_specific' 1 ((PROCEDURE UNKNOWN-INTENT
MODULE-PROC DECL UNKNOWN 0 0 FUNCTION) (LOGICAL 4 0 0 LOGICAL ()) 262 0
(263 264 265) () 32 () () () 0 0)
21 'get_scalar_boundary_condition_count' 'boundary_conditions'
'get_scalar_boundary_condition_count' 1 ((PROCEDURE UNKNOWN-INTENT
MODULE-PROC DECL UNKNOWN 0 0 FUNCTION IMPLICIT_PURE) (INTEGER 4 0 0
INTEGER ()) 266 0 (267) () 21 () () () 0 0)
20 'get_vector_boundary_condition_count' 'boundary_conditions'
'get_vector_boundary_condition_count' 1 ((PROCEDURE UNKNOWN-INTENT
MODULE-PROC DECL UNKNOWN 0 0 FUNCTION IMPLICIT_PURE) (INTEGER 4 0 0
INTEGER ()) 268 0 (269) () 20 () () () 0 0)
23 'get_entire_scalar_boundary_condition' 'boundary_conditions'
'get_entire_scalar_boundary_condition' 1 ((PROCEDURE UNKNOWN-INTENT
MODULE-PROC DECL UNKNOWN 0 0 SUBROUTINE ALWAYS_EXPLICIT) (UNKNOWN 0 0 0
UNKNOWN ()) 270 0 (271 272 273 274 275 276) () 0 () () () 0 0)
22 'get_entire_vector_boundary_condition' 'boundary_conditions'
'get_entire_vector_boundary_condition' 1 ((PROCEDURE UNKNOWN-INTENT
MODULE-PROC DECL UNKNOWN 0 0 SUBROUTINE ALWAYS_EXPLICIT) (UNKNOWN 0 0 0
UNKNOWN ()) 277 0 (278 279 280 281 282) () 0 () () () 0 0)
25 'get_scalar_boundary_condition_nodes' 'boundary_conditions'
'get_scalar_boundary_condition_nodes' 1 ((PROCEDURE UNKNOWN-INTENT
MODULE-PROC DECL UNKNOWN 0 0 SUBROUTINE ALWAYS_EXPLICIT) (UNKNOWN 0 0 0
UNKNOWN ()) 283 0 (284 285 286) () 0 () () () 0 0)
24 'get_vector_boundary_condition_nodes' 'boundary_conditions'
'get_vector_boundary_condition_nodes' 1 ((PROCEDURE UNKNOWN-INTENT
MODULE-PROC DECL UNKNOWN 0 0 SUBROUTINE ALWAYS_EXPLICIT) (UNKNOWN 0 0 0
UNKNOWN ()) 287 0 (288 289 290) () 0 () () () 0 0)
45 'set_reference_node_scalar' 'boundary_conditions'
'set_reference_node_scalar' 1 ((PROCEDURE UNKNOWN-INTENT MODULE-PROC
DECL UNKNOWN 0 0 SUBROUTINE ALWAYS_EXPLICIT) (UNKNOWN 0 0 0 UNKNOWN ())
291 0 (292 293 294 295 296) () 0 () () () 0 0)
44 'set_reference_node_vector_petsc' 'boundary_conditions'
'set_reference_node_vector_petsc' 1 ((PROCEDURE UNKNOWN-INTENT
MODULE-PROC DECL UNKNOWN 0 0 SUBROUTINE ALWAYS_EXPLICIT) (UNKNOWN 0 0 0
UNKNOWN ()) 297 0 (298 299 300 301 302) () 0 () () () 0 0)
41 'has_boundary_condition_scalar' 'boundary_conditions'
'has_boundary_condition_scalar' 1 ((PROCEDURE UNKNOWN-INTENT MODULE-PROC
DECL UNKNOWN 0 0 FUNCTION IMPLICIT_PURE) (LOGICAL 4 0 0 LOGICAL ()) 303
0 (304 305) () 41 () () () 0 0)
40 'has_boundary_condition_vector' 'boundary_conditions'
'has_boundary_condition_vector' 1 ((PROCEDURE UNKNOWN-INTENT MODULE-PROC
DECL UNKNOWN 0 0 FUNCTION IMPLICIT_PURE) (LOGICAL 4 0 0 LOGICAL ()) 306
0 (307 308) () 40 () () () 0 0)
31 'has_boundary_condition_name_scalar' 'boundary_conditions'
'has_boundary_condition_name_scalar' 1 ((PROCEDURE UNKNOWN-INTENT
MODULE-PROC DECL UNKNOWN 0 0 FUNCTION IMPLICIT_PURE) (LOGICAL 4 0 0
LOGICAL ()) 309 0 (310 311) () 31 () () () 0 0)
30 'has_boundary_condition_name_vector' 'boundary_conditions'
'has_boundary_condition_name_vector' 1 ((PROCEDURE UNKNOWN-INTENT
MODULE-PROC DECL UNKNOWN 0 0 FUNCTION IMPLICIT_PURE) (LOGICAL 4 0 0
LOGICAL ()) 312 0 (313 314) () 30 () () () 0 0)
47 'set_dirichlet_consistent_scalar' 'boundary_conditions'
'set_dirichlet_consistent_scalar' 1 ((PROCEDURE UNKNOWN-INTENT
MODULE-PROC DECL UNKNOWN 0 0 SUBROUTINE) (UNKNOWN 0 0 0 UNKNOWN ()) 315
0 (316) () 0 () () () 0 0)
46 'set_dirichlet_consistent_vector' 'boundary_conditions'
'set_dirichlet_consistent_vector' 1 ((PROCEDURE UNKNOWN-INTENT
MODULE-PROC DECL UNKNOWN 0 0 SUBROUTINE) (UNKNOWN 0 0 0 UNKNOWN ()) 317
0 (318) () 0 () () () 0 0)
13 'apply_dirichlet_conditions_scalar' 'boundary_conditions'
'apply_dirichlet_conditions_scalar' 1 ((PROCEDURE UNKNOWN-INTENT
MODULE-PROC DECL UNKNOWN 0 0 SUBROUTINE ALWAYS_EXPLICIT) (UNKNOWN 0 0 0
UNKNOWN ()) 319 0 (320 321 322 323) () 0 () () () 0 0)
12 'apply_dirichlet_conditions_scalar_lumped' 'boundary_conditions'
'apply_dirichlet_conditions_scalar_lumped' 1 ((PROCEDURE UNKNOWN-INTENT
MODULE-PROC DECL UNKNOWN 0 0 SUBROUTINE ALWAYS_EXPLICIT) (UNKNOWN 0 0 0
UNKNOWN ()) 324 0 (325 326 327 328) () 0 () () () 0 0)
11 'apply_dirichlet_conditions_scalars_petsc_csr' 'boundary_conditions'
'apply_dirichlet_conditions_scalars_petsc_csr' 1 ((PROCEDURE
UNKNOWN-INTENT MODULE-PROC DECL UNKNOWN 0 0 SUBROUTINE ALWAYS_EXPLICIT)
(UNKNOWN 0 0 0 UNKNOWN ()) 329 0 (330 331 332 333) () 0 () () () 0 0)
10 'apply_dirichlet_conditions_vector' 'boundary_conditions'
'apply_dirichlet_conditions_vector' 1 ((PROCEDURE UNKNOWN-INTENT
MODULE-PROC DECL UNKNOWN 0 0 SUBROUTINE ALWAYS_EXPLICIT) (UNKNOWN 0 0 0
UNKNOWN ()) 334 0 (335 336 337 338) () 0 () () () 0 0)
9 'apply_dirichlet_conditions_vector_petsc_csr' 'boundary_conditions'
'apply_dirichlet_conditions_vector_petsc_csr' 1 ((PROCEDURE
UNKNOWN-INTENT MODULE-PROC DECL UNKNOWN 0 0 SUBROUTINE ALWAYS_EXPLICIT)
(UNKNOWN 0 0 0 UNKNOWN ()) 339 0 (340 341 342 343) () 0 () () () 0 0)
8 'apply_dirichlet_conditions_vector_component' 'boundary_conditions'
'apply_dirichlet_conditions_vector_component' 1 ((PROCEDURE
UNKNOWN-INTENT MODULE-PROC DECL UNKNOWN 0 0 SUBROUTINE ALWAYS_EXPLICIT)
(UNKNOWN 0 0 0 UNKNOWN ()) 344 0 (345 346 347 348 349) () 0 () () () 0 0)
7 'apply_dirichlet_conditions_vector_lumped' 'boundary_conditions'
'apply_dirichlet_conditions_vector_lumped' 1 ((PROCEDURE UNKNOWN-INTENT
MODULE-PROC DECL UNKNOWN 0 0 SUBROUTINE ALWAYS_EXPLICIT) (UNKNOWN 0 0 0
UNKNOWN ()) 350 0 (351 352 353 354) () 0 () () () 0 0)
6 'apply_dirichlet_conditions_vector_component_lumped'
'boundary_conditions' 'apply_dirichlet_conditions_vector_component_lumped'
1 ((PROCEDURE UNKNOWN-INTENT MODULE-PROC DECL UNKNOWN 0 0 SUBROUTINE
ALWAYS_EXPLICIT) (UNKNOWN 0 0 0 UNKNOWN ()) 355 0 (356 357 358 359 360)
() 0 () () () 0 0)
135 'field' '' 'field' 134 ((VARIABLE INOUT UNKNOWN-PROC UNKNOWN UNKNOWN
0 0 DUMMY) (DERIVED 226 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
138 'boundary_ids' '' 'boundary_ids' 134 ((VARIABLE IN UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0 DIMENSION DUMMY) (INTEGER 4 0 0 INTEGER ()) 0 0 () (
1 0 ASSUMED_SHAPE (CONSTANT (INTEGER 4 0 0 INTEGER ()) 0 '1') ()) 0 () ()
() 0 0)
139 'option_path' '' 'option_path' 134 ((VARIABLE IN UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0 OPTIONAL DUMMY) (CHARACTER 1 0 0 CHARACTER (())) 0 0
() () 0 () () () 0 0)
136 'name' '' 'name' 134 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DUMMY) (CHARACTER 1 0 0 CHARACTER (())) 0 0 () () 0 () () () 0 0)
137 'type' '' 'type' 134 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DUMMY) (CHARACTER 1 0 0 CHARACTER (())) 0 0 () () 0 () () () 0 0)
154 'field' '' 'field' 153 ((VARIABLE INOUT UNKNOWN-PROC UNKNOWN UNKNOWN
0 0 DUMMY) (DERIVED 226 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
155 'name' '' 'name' 153 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DUMMY) (CHARACTER 1 0 0 CHARACTER (())) 0 0 () () 0 () () () 0 0)
156 'type' '' 'type' 153 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DUMMY) (CHARACTER 1 0 0 CHARACTER (())) 0 0 () () 0 () () () 0 0)
157 'surface_element_list' '' 'surface_element_list' 153 ((VARIABLE IN
UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 DIMENSION DUMMY) (INTEGER 4 0 0 INTEGER
()) 0 0 () (1 0 ASSUMED_SHAPE (CONSTANT (INTEGER 4 0 0 INTEGER ()) 0 '1')
()) 0 () () () 0 0)
158 'option_path' '' 'option_path' 153 ((VARIABLE IN UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0 OPTIONAL DUMMY) (CHARACTER 1 0 0 CHARACTER (())) 0 0
() () 0 () () () 0 0)
143 'type' '' 'type' 140 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DUMMY) (CHARACTER 1 0 0 CHARACTER (())) 0 0 () () 0 () () () 0 0)
142 'name' '' 'name' 140 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DUMMY) (CHARACTER 1 0 0 CHARACTER (())) 0 0 () () 0 () () () 0 0)
141 'field' '' 'field' 140 ((VARIABLE INOUT UNKNOWN-PROC UNKNOWN UNKNOWN
0 0 DUMMY) (DERIVED 232 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
144 'boundary_ids' '' 'boundary_ids' 140 ((VARIABLE IN UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0 DIMENSION DUMMY) (INTEGER 4 0 0 INTEGER ()) 0 0 () (
1 0 ASSUMED_SHAPE (CONSTANT (INTEGER 4 0 0 INTEGER ()) 0 '1') ()) 0 () ()
() 0 0)
145 'applies' '' 'applies' 140 ((VARIABLE IN UNKNOWN-PROC UNKNOWN
UNKNOWN 0 0 DIMENSION OPTIONAL DUMMY) (LOGICAL 4 0 0 LOGICAL ()) 0 0 ()
(1 0 ASSUMED_SHAPE (CONSTANT (INTEGER 4 0 0 INTEGER ()) 0 '1') ()) 0 ()
() () 0 0)
146 'option_path' '' 'option_path' 140 ((VARIABLE IN UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0 OPTIONAL DUMMY) (CHARACTER 1 0 0 CHARACTER (())) 0 0
() () 0 () () () 0 0)
160 'field' '' 'field' 159 ((VARIABLE INOUT UNKNOWN-PROC UNKNOWN UNKNOWN
0 0 DUMMY) (DERIVED 232 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
161 'name' '' 'name' 159 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DUMMY) (CHARACTER 1 0 0 CHARACTER (())) 0 0 () () 0 () () () 0 0)
162 'type' '' 'type' 159 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DUMMY) (CHARACTER 1 0 0 CHARACTER (())) 0 0 () () 0 () () () 0 0)
163 'surface_element_list' '' 'surface_element_list' 159 ((VARIABLE IN
UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 DIMENSION DUMMY) (INTEGER 4 0 0 INTEGER
()) 0 0 () (1 0 ASSUMED_SHAPE (CONSTANT (INTEGER 4 0 0 INTEGER ()) 0 '1')
()) 0 () () () 0 0)
164 'applies' '' 'applies' 159 ((VARIABLE IN UNKNOWN-PROC UNKNOWN
UNKNOWN 0 0 DIMENSION OPTIONAL DUMMY) (LOGICAL 4 0 0 LOGICAL ()) 0 0 ()
(1 0 ASSUMED_SHAPE (CONSTANT (INTEGER 4 0 0 INTEGER ()) 0 '1') ()) 0 ()
() () 0 0)
165 'option_path' '' 'option_path' 159 ((VARIABLE IN UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0 OPTIONAL DUMMY) (CHARACTER 1 0 0 CHARACTER (())) 0 0
() () 0 () () () 0 0)
149 'name' '' 'name' 147 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DUMMY) (CHARACTER 1 0 0 CHARACTER (())) 0 0 () () 0 () () () 0 0)
148 'field' '' 'field' 147 ((VARIABLE INOUT UNKNOWN-PROC UNKNOWN UNKNOWN
0 0 DUMMY) (DERIVED 226 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
151 'field' '' 'field' 150 ((VARIABLE INOUT UNKNOWN-PROC UNKNOWN UNKNOWN
0 0 DUMMY) (DERIVED 232 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
152 'name' '' 'name' 150 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DUMMY) (CHARACTER 1 0 0 CHARACTER (())) 0 0 () () 0 () () () 0 0)
203 'field' '' 'field' 202 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0
0 DUMMY) (DERIVED 226 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
204 'n' '' 'n' 202 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 DUMMY)
(INTEGER 4 0 0 INTEGER ()) 0 0 () () 0 () () () 0 0)
205 'surface_field' '' 'surface_field' 202 ((VARIABLE IN UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0 DUMMY) (DERIVED 226 0 0 DERIVED ()) 0 0 () () 0 () ()
() 0 0)
207 'field' '' 'field' 206 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0
0 DUMMY) (DERIVED 232 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
208 'n' '' 'n' 206 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 DUMMY)
(INTEGER 4 0 0 INTEGER ()) 0 0 () () 0 () () () 0 0)
209 'surface_field' '' 'surface_field' 206 ((VARIABLE IN UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0 DUMMY) (DERIVED 232 0 0 DERIVED ()) 0 0 () () 0 () ()
() 0 0)
219 'field' '' 'field' 218 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0
0 DUMMY) (DERIVED 232 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
220 'n' '' 'n' 218 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 DUMMY)
(INTEGER 4 0 0 INTEGER ()) 0 0 () () 0 () () () 0 0)
221 'surface_field' '' 'surface_field' 218 ((VARIABLE IN UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0 DUMMY) (DERIVED 226 0 0 DERIVED ()) 0 0 () () 0 () ()
() 0 0)
211 'field' '' 'field' 210 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0
0 DUMMY) (DERIVED 226 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
212 'name' '' 'name' 210 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DUMMY) (CHARACTER 1 0 0 CHARACTER (())) 0 0 () () 0 () () () 0 0)
213 'surface_field' '' 'surface_field' 210 ((VARIABLE IN UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0 DUMMY) (DERIVED 226 0 0 DERIVED ()) 0 0 () () 0 () ()
() 0 0)
215 'field' '' 'field' 214 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0
0 DUMMY) (DERIVED 232 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
216 'name' '' 'name' 214 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DUMMY) (CHARACTER 1 0 0 CHARACTER (())) 0 0 () () 0 () () () 0 0)
217 'surface_field' '' 'surface_field' 214 ((VARIABLE IN UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0 DUMMY) (DERIVED 232 0 0 DERIVED ()) 0 0 () () 0 () ()
() 0 0)
223 'field' '' 'field' 222 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0
0 DUMMY) (DERIVED 232 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
224 'name' '' 'name' 222 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DUMMY) (CHARACTER 1 0 0 CHARACTER (())) 0 0 () () 0 () () () 0 0)
225 'surface_field' '' 'surface_field' 222 ((VARIABLE IN UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0 DUMMY) (DERIVED 226 0 0 DERIVED ()) 0 0 () () 0 () ()
() 0 0)
228 'field' '' 'field' 227 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0
0 DUMMY) (DERIVED 226 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
229 'n' '' 'n' 227 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 DUMMY)
(INTEGER 4 0 0 INTEGER ()) 0 0 () () 0 () () () 0 0)
197 'surface_element_list' '' 'surface_element_list' 193 ((VARIABLE
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 DIMENSION OPTIONAL
POINTER DUMMY) (INTEGER 4 0 0 INTEGER ()) 0 0 () (1 0 DEFERRED () ()) 0
() () () 0 0)
198 'surface_node_list' '' 'surface_node_list' 193 ((VARIABLE
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 DIMENSION OPTIONAL
POINTER DUMMY) (INTEGER 4 0 0 INTEGER ()) 0 0 () (1 0 DEFERRED () ()) 0
() () () 0 0)
199 'applies' '' 'applies' 193 ((VARIABLE OUT UNKNOWN-PROC UNKNOWN
UNKNOWN 0 0 DIMENSION OPTIONAL DUMMY) (LOGICAL 4 0 0 LOGICAL ()) 0 0 ()
(1 0 ASSUMED_SHAPE (CONSTANT (INTEGER 4 0 0 INTEGER ()) 0 '1') ()) 0 ()
() () 0 0)
200 'surface_mesh' '' 'surface_mesh' 193 ((VARIABLE UNKNOWN-INTENT
UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 OPTIONAL POINTER DUMMY) (DERIVED 72 0 0
DERIVED ()) 0 0 () () 0 () () () 0 0)
201 'option_path' '' 'option_path' 193 ((VARIABLE OUT UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0 OPTIONAL DUMMY) (CHARACTER 1 0 0 CHARACTER (())) 0 0
() () 0 () () () 0 0)
77 'mesh' '' 'mesh' 76 ((VARIABLE INOUT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DUMMY) (DERIVED 72 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
78 'periodic_bc_list' '' 'periodic_bc_list' 76 ((VARIABLE OUT
UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 DIMENSION DUMMY) (LOGICAL 4 0 0 LOGICAL
()) 0 0 () (1 0 ASSUMED_SHAPE (CONSTANT (INTEGER 4 0 0 INTEGER ()) 0 '1')
()) 0 () () () 0 0)
271 'field' '' 'field' 270 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0
0 TARGET DUMMY) (DERIVED 226 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
272 'types' '' 'types' 270 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0
0 DIMENSION DUMMY) (CHARACTER 1 0 0 CHARACTER (())) 0 0 () (1 0
ASSUMED_SHAPE (CONSTANT (INTEGER 4 0 0 INTEGER ()) 0 '1') ()) 0 () () ()
0 0)
273 'boundary_value' '' 'boundary_value' 270 ((VARIABLE OUT UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0 DUMMY) (DERIVED 226 0 0 DERIVED ()) 0 0 () () 0 () ()
() 0 0)
274 'bc_type_list' '' 'bc_type_list' 270 ((VARIABLE OUT UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0 DIMENSION DUMMY) (INTEGER 4 0 0 INTEGER ()) 0 0 () (
1 0 ASSUMED_SHAPE (CONSTANT (INTEGER 4 0 0 INTEGER ()) 0 '1') ()) 0 () ()
() 0 0)
275 'bc_number_list' '' 'bc_number_list' 270 ((VARIABLE OUT UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0 DIMENSION OPTIONAL DUMMY) (INTEGER 4 0 0 INTEGER ())
0 0 () (1 0 ASSUMED_SHAPE (CONSTANT (INTEGER 4 0 0 INTEGER ()) 0 '1') ())
0 () () () 0 0)
276 'boundary_second_value' '' 'boundary_second_value' 270 ((VARIABLE
OUT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 OPTIONAL DUMMY) (DERIVED 226 0 0
DERIVED ()) 0 0 () () 0 () () () 0 0)
285 'types' '' 'types' 283 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0
0 DIMENSION DUMMY) (CHARACTER 1 0 0 CHARACTER (())) 0 0 () (1 0
ASSUMED_SHAPE (CONSTANT (INTEGER 4 0 0 INTEGER ()) 0 '1') ()) 0 () () ()
0 0)
278 'field' '' 'field' 277 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0
0 TARGET DUMMY) (DERIVED 232 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
279 'types' '' 'types' 277 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0
0 DIMENSION DUMMY) (CHARACTER 1 0 0 CHARACTER (())) 0 0 () (1 0
ASSUMED_SHAPE (CONSTANT (INTEGER 4 0 0 INTEGER ()) 0 '1') ()) 0 () () ()
0 0)
280 'boundary_value' '' 'boundary_value' 277 ((VARIABLE OUT UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0 DUMMY) (DERIVED 232 0 0 DERIVED ()) 0 0 () () 0 () ()
() 0 0)
281 'bc_type_list' '' 'bc_type_list' 277 ((VARIABLE OUT UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0 DIMENSION DUMMY) (INTEGER 4 0 0 INTEGER ()) 0 0 () (
2 0 ASSUMED_SHAPE (CONSTANT (INTEGER 4 0 0 INTEGER ()) 0 '1') () (
CONSTANT (INTEGER 4 0 0 INTEGER ()) 0 '1') ()) 0 () () () 0 0)
282 'bc_number_list' '' 'bc_number_list' 277 ((VARIABLE OUT UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0 DIMENSION OPTIONAL DUMMY) (INTEGER 4 0 0 INTEGER ())
0 0 () (2 0 ASSUMED_SHAPE (CONSTANT (INTEGER 4 0 0 INTEGER ()) 0 '1') ()
(CONSTANT (INTEGER 4 0 0 INTEGER ()) 0 '1') ()) 0 () () () 0 0)
284 'field' '' 'field' 283 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0
0 DUMMY) (DERIVED 226 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
74 'mesh' '' 'mesh' 73 ((VARIABLE INOUT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DUMMY) (DERIVED 72 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
286 'bc_type_node_list' '' 'bc_type_node_list' 283 ((VARIABLE OUT
UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 DIMENSION DUMMY) (INTEGER 4 0 0 INTEGER
()) 0 0 () (1 0 ASSUMED_SHAPE (CONSTANT (INTEGER 4 0 0 INTEGER ()) 0 '1')
()) 0 () () () 0 0)
289 'types' '' 'types' 287 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0
0 DIMENSION DUMMY) (CHARACTER 1 0 0 CHARACTER (())) 0 0 () (1 0
ASSUMED_SHAPE (CONSTANT (INTEGER 4 0 0 INTEGER ()) 0 '1') ()) 0 () () ()
0 0)
290 'bc_type_node_list' '' 'bc_type_node_list' 287 ((VARIABLE OUT
UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 DIMENSION DUMMY) (INTEGER 4 0 0 INTEGER
()) 0 0 () (2 0 ASSUMED_SHAPE (CONSTANT (INTEGER 4 0 0 INTEGER ()) 0 '1')
() (CONSTANT (INTEGER 4 0 0 INTEGER ()) 0 '1') ()) 0 () () () 0 0)
288 'field' '' 'field' 287 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0
0 DUMMY) (DERIVED 232 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
293 'node' '' 'node' 291 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DUMMY) (INTEGER 4 0 0 INTEGER ()) 0 0 () () 0 () () () 0 0)
103 'adjacency_cache' 'fields_data_types' 'adjacency_cache' 1 ((DERIVED
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 POINTER_COMP) (UNKNOWN 0
0 0 UNKNOWN ()) 0 0 () () 0 ((361 'nnlist' (DERIVED 362 0 0 DERIVED ())
() (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 POINTER)
UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0)) (363 'nelist' (
DERIVED 362 0 0 DERIVED ()) () (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0 POINTER) UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN
()) 0)) (364 'eelist' (DERIVED 362 0 0 DERIVED ()) () (UNKNOWN-FL
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 POINTER) UNKNOWN-ACCESS
(NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0))) PUBLIC (() () () ()) () 0 0
37419158)
99 'mesh_faces' 'fields_data_types' 'mesh_faces' 1 ((DERIVED
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 POINTER_COMP) (UNKNOWN 0
0 0 UNKNOWN ()) 0 0 () () 0 ((365 'shape' (DERIVED 89 0 0 DERIVED ()) ()
(UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 POINTER)
UNKNOWN-ACCESS ()) (366 'face_list' (DERIVED 367 0 0 DERIVED ()) () (
UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0)
UNKNOWN-ACCESS (STRUCTURE (DERIVED 367 0 0 DERIVED ()) 0 (((STRUCTURE (
DERIVED 362 0 0 DERIVED ()) 0 (((NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0) ())
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
())) ())) (368 'face_lno' (INTEGER 4 0 0 INTEGER ()) (1 0 DEFERRED () ())
(UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 DIMENSION
POINTER) UNKNOWN-ACCESS ()) (369 'surface_mesh' (DERIVED 72 0 0 DERIVED
()) () (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0)
UNKNOWN-ACCESS (STRUCTURE (DERIVED 72 0 0 DERIVED ()) 0 ((() ()) ((
CONSTANT (LOGICAL 4 0 0 LOGICAL ()) 0 1) ()) ((STRUCTURE (DERIVED 89 0 0
DERIVED ()) 0 ((() ()) (() ()) (() ()) (() ()) ((NULL (UNKNOWN 0 0 0
UNKNOWN ()) 0) ()) ((NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0) ()) ((NULL (
UNKNOWN 0 0 0 UNKNOWN ()) 0) ()) ((NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0) ())
((NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0) ()) ((NULL (UNKNOWN 0 0 0 UNKNOWN ())
0) ()) ((NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0) ()) ((STRUCTURE (DERIVED 90
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
0) ())) ())) (370 'surface_node_list' (INTEGER 4 0 0 INTEGER ()) (1 0
DEFERRED () ()) (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN
0 0 DIMENSION POINTER) UNKNOWN-ACCESS ()) (371 'face_element_list' (
INTEGER 4 0 0 INTEGER ()) (1 0 DEFERRED () ()) (UNKNOWN-FL
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 DIMENSION POINTER)
UNKNOWN-ACCESS ()) (372 'boundary_ids' (INTEGER 4 0 0 INTEGER ()) (1 0
DEFERRED () ()) (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN
0 0 DIMENSION POINTER) UNKNOWN-ACCESS ()) (373 'coplanar_ids' (INTEGER 4
0 0 INTEGER ()) (1 0 DEFERRED () ()) (UNKNOWN-FL UNKNOWN-INTENT
UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 DIMENSION POINTER) UNKNOWN-ACCESS (
NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0)) (374 'dg_surface_mesh' (DERIVED 72 0
0 DERIVED ()) () (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN
0 0 POINTER) UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0)) (375
'has_internal_boundaries' (LOGICAL 4 0 0 LOGICAL ()) () (UNKNOWN-FL
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0) UNKNOWN-ACCESS (
CONSTANT (LOGICAL 4 0 0 LOGICAL ()) 0 0))) PUBLIC (() () () ()) () 0 0
11936185)
101 'mesh_subdomain_mesh' 'fields_data_types' 'mesh_subdomain_mesh' 1 (
(DERIVED UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 POINTER_COMP) (
UNKNOWN 0 0 0 UNKNOWN ()) 0 0 () () 0 ((376 'element_list' (INTEGER 4 0
0 INTEGER ()) (1 0 DEFERRED () ()) (UNKNOWN-FL UNKNOWN-INTENT
UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 DIMENSION POINTER) UNKNOWN-ACCESS ()) (
377 'node_list' (INTEGER 4 0 0 INTEGER ()) (1 0 DEFERRED () ()) (
UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 DIMENSION
POINTER) UNKNOWN-ACCESS ())) PUBLIC (() () () ()) () 0 0 85302181)
226 'scalar_field' 'fields_data_types' 'scalar_field' 1 ((DERIVED
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 POINTER_COMP) (UNKNOWN 0
0 0 UNKNOWN ()) 0 0 () () 0 ((378 'val' (REAL 8 0 0 REAL ()) (1 0
DEFERRED () ()) (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN
0 0 DIMENSION POINTER) UNKNOWN-ACCESS ()) (379 'val_stride' (INTEGER 4 0
0 INTEGER ()) () (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN
0 0) UNKNOWN-ACCESS (CONSTANT (INTEGER 4 0 0 INTEGER ()) 0 '1')) (380
'wrapped' (LOGICAL 4 0 0 LOGICAL ()) () (UNKNOWN-FL UNKNOWN-INTENT
UNKNOWN-PROC UNKNOWN UNKNOWN 0 0) UNKNOWN-ACCESS (CONSTANT (LOGICAL 4 0
0 LOGICAL ()) 0 1)) (381 'field_type' (INTEGER 4 0 0 INTEGER ()) () (
UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0)
UNKNOWN-ACCESS (CONSTANT (INTEGER 4 0 0 INTEGER ()) 0 '0')) (382 'bc' (
DERIVED 383 0 0 DERIVED ()) () (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0 POINTER) UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN
()) 0)) (384 'name' (CHARACTER 1 0 0 CHARACTER ((CONSTANT (INTEGER 4 0 0
INTEGER ()) 0 '101'))) () (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0) UNKNOWN-ACCESS ()) (385 'option_path' (CHARACTER 1
0 0 CHARACTER ((CONSTANT (INTEGER 4 0 0 INTEGER ()) 0 '8192'))) () (
UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0)
UNKNOWN-ACCESS (CONSTANT (CHARACTER 1 0 0 CHARACTER (())) 0 8192
'/uninitialised_path/                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            '))
(386 'mesh' (DERIVED 72 0 0 DERIVED ()) () (UNKNOWN-FL UNKNOWN-INTENT
UNKNOWN-PROC UNKNOWN UNKNOWN 0 0) UNKNOWN-ACCESS (STRUCTURE (DERIVED 72
0 0 DERIVED ()) 0 ((() ()) ((CONSTANT (LOGICAL 4 0 0 LOGICAL ()) 0 1) ())
((STRUCTURE (DERIVED 89 0 0 DERIVED ()) 0 ((() ()) (() ()) (() ()) (() ())
((NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0) ()) ((NULL (UNKNOWN 0 0 0 UNKNOWN ())
0) ()) ((NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0) ()) ((NULL (UNKNOWN 0 0 0
UNKNOWN ()) 0) ()) ((NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0) ()) ((NULL (
UNKNOWN 0 0 0 UNKNOWN ()) 0) ()) ((NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0) ())
((STRUCTURE (DERIVED 90 0 0 DERIVED ()) 0 ((() ()) (() ()) (() ()) (() ())
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
0) ())) ())) (387 'refcount' (DERIVED 97 0 0 DERIVED ()) () (UNKNOWN-FL
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 POINTER) UNKNOWN-ACCESS
(NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0)) (388 'aliased' (LOGICAL 4 0 0
LOGICAL ()) () (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0
0) UNKNOWN-ACCESS (CONSTANT (LOGICAL 4 0 0 LOGICAL ()) 0 0)) (389
'py_locweight' (REAL 8 0 0 REAL ()) (2 0 DEFERRED () () () ()) (
UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 DIMENSION
POINTER) UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0)) (390
'py_func' (CHARACTER 1 0 0 CHARACTER ((CONSTANT (INTEGER 4 0 0 INTEGER ())
0 '8192'))) () (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0
0) UNKNOWN-ACCESS ()) (391 'py_positions' (DERIVED 232 0 0 DERIVED ()) ()
(UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 POINTER)
UNKNOWN-ACCESS ()) (392 'py_positions_same_mesh' (LOGICAL 4 0 0 LOGICAL
()) () (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0)
UNKNOWN-ACCESS ()) (393 'py_dim' (INTEGER 4 0 0 INTEGER ()) () (
UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0)
UNKNOWN-ACCESS ()) (394 'py_positions_shape' (DERIVED 89 0 0 DERIVED ())
() (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 POINTER)
UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0))) PUBLIC (() () () ())
() 0 0 64912956)
232 'vector_field' 'fields_data_types' 'vector_field' 1 ((DERIVED
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 POINTER_COMP) (UNKNOWN 0
0 0 UNKNOWN ()) 0 0 () () 0 ((395 'val' (REAL 8 0 0 REAL ()) (2 0
DEFERRED () () () ()) (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN
UNKNOWN 0 0 DIMENSION POINTER) UNKNOWN-ACCESS ()) (396 'wrapped' (
LOGICAL 4 0 0 LOGICAL ()) () (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0) UNKNOWN-ACCESS (CONSTANT (LOGICAL 4 0 0 LOGICAL ())
0 1)) (397 'field_type' (INTEGER 4 0 0 INTEGER ()) () (UNKNOWN-FL
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0) UNKNOWN-ACCESS (
CONSTANT (INTEGER 4 0 0 INTEGER ()) 0 '0')) (398 'bc' (DERIVED 399 0 0
DERIVED ()) () (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0
0 POINTER) UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0)) (400 'name'
(CHARACTER 1 0 0 CHARACTER ((CONSTANT (INTEGER 4 0 0 INTEGER ()) 0 '101')))
() (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0)
UNKNOWN-ACCESS ()) (401 'dim' (INTEGER 4 0 0 INTEGER ()) () (UNKNOWN-FL
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0) UNKNOWN-ACCESS ()) (
402 'option_path' (CHARACTER 1 0 0 CHARACTER ((CONSTANT (INTEGER 4 0 0
INTEGER ()) 0 '8192'))) () (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0) UNKNOWN-ACCESS (CONSTANT (CHARACTER 1 0 0 CHARACTER
(())) 0 8192
'/uninitialised_path/                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            '))
(403 'mesh' (DERIVED 72 0 0 DERIVED ()) () (UNKNOWN-FL UNKNOWN-INTENT
UNKNOWN-PROC UNKNOWN UNKNOWN 0 0) UNKNOWN-ACCESS (STRUCTURE (DERIVED 72
0 0 DERIVED ()) 0 ((() ()) ((CONSTANT (LOGICAL 4 0 0 LOGICAL ()) 0 1) ())
((STRUCTURE (DERIVED 89 0 0 DERIVED ()) 0 ((() ()) (() ()) (() ()) (() ())
((NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0) ()) ((NULL (UNKNOWN 0 0 0 UNKNOWN ())
0) ()) ((NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0) ()) ((NULL (UNKNOWN 0 0 0
UNKNOWN ()) 0) ()) ((NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0) ()) ((NULL (
UNKNOWN 0 0 0 UNKNOWN ()) 0) ()) ((NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0) ())
((STRUCTURE (DERIVED 90 0 0 DERIVED ()) 0 ((() ()) (() ()) (() ()) (() ())
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
0) ())) ())) (404 'refcount' (DERIVED 97 0 0 DERIVED ()) () (UNKNOWN-FL
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 POINTER) UNKNOWN-ACCESS
(NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0)) (405 'aliased' (LOGICAL 4 0 0
LOGICAL ()) () (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0
0) UNKNOWN-ACCESS (CONSTANT (LOGICAL 4 0 0 LOGICAL ()) 0 0)) (406 'picker'
(DERIVED 407 0 0 DERIVED ()) () (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0 POINTER) UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN
()) 0))) PUBLIC (() () () ()) () 0 0 43598963)
90 'quadrature_type' 'quadrature' 'quadrature_type' 1 ((DERIVED
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 POINTER_COMP) (UNKNOWN 0
0 0 UNKNOWN ()) 0 0 () () 0 ((408 'dim' (INTEGER 4 0 0 INTEGER ()) () (
UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0)
UNKNOWN-ACCESS ()) (409 'degree' (INTEGER 4 0 0 INTEGER ()) () (
UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0)
UNKNOWN-ACCESS ()) (410 'vertices' (INTEGER 4 0 0 INTEGER ()) () (
UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0)
UNKNOWN-ACCESS ()) (411 'ngi' (INTEGER 4 0 0 INTEGER ()) () (UNKNOWN-FL
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0) UNKNOWN-ACCESS ()) (
412 'weight' (REAL 8 0 0 REAL ()) (1 0 DEFERRED () ()) (UNKNOWN-FL
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 DIMENSION POINTER)
UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0)) (413 'l' (REAL 8 0 0
REAL ()) (2 0 DEFERRED () () () ()) (UNKNOWN-FL UNKNOWN-INTENT
UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 DIMENSION POINTER) UNKNOWN-ACCESS (
NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0)) (414 'name' (CHARACTER 1 0 0
CHARACTER ((CONSTANT (INTEGER 4 0 0 INTEGER ()) 0 '0'))) () (UNKNOWN-FL
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0) UNKNOWN-ACCESS ()) (
415 'refcount' (DERIVED 97 0 0 DERIVED ()) () (UNKNOWN-FL UNKNOWN-INTENT
UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 POINTER) UNKNOWN-ACCESS (NULL (UNKNOWN
0 0 0 UNKNOWN ()) 0)) (416 'family' (INTEGER 4 0 0 INTEGER ()) () (
UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0)
UNKNOWN-ACCESS ())) PUBLIC (() () () ()) () 0 0 59837722)
367 'csr_matrix' 'sparse_tools' 'csr_matrix' 1 ((DERIVED UNKNOWN-INTENT
UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 POINTER_COMP) (UNKNOWN 0 0 0 UNKNOWN ())
0 0 () () 0 ((417 'sparsity' (DERIVED 362 0 0 DERIVED ()) () (
UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0)
UNKNOWN-ACCESS (STRUCTURE (DERIVED 362 0 0 DERIVED ()) 0 (((NULL (
UNKNOWN 0 0 0 UNKNOWN ()) 0) ()) ((NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0) ())
((NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0) ()) (() ()) ((NULL (UNKNOWN 0 0 0
UNKNOWN ()) 0) ()) ((NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0) ()) ((NULL (
UNKNOWN 0 0 0 UNKNOWN ()) 0) ()) ((CONSTANT (CHARACTER 1 0 0 CHARACTER (
())) 0 101
'                                                                                                     ')
()) ((CONSTANT (LOGICAL 4 0 0 LOGICAL ()) 0 0) ()) ((CONSTANT (LOGICAL 4
0 0 LOGICAL ()) 0 0) ())) ())) (418 'val' (REAL 8 0 0 REAL ()) (1 0
DEFERRED () ()) (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN
0 0 DIMENSION POINTER) UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0))
(419 'ival' (INTEGER 4 0 0 INTEGER ()) (1 0 DEFERRED () ()) (UNKNOWN-FL
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 DIMENSION POINTER)
UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0)) (420 'clone' (
LOGICAL 4 0 0 LOGICAL ()) () (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0) UNKNOWN-ACCESS (CONSTANT (LOGICAL 4 0 0 LOGICAL ())
0 0)) (421 'external_val' (LOGICAL 4 0 0 LOGICAL ()) () (UNKNOWN-FL
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0) UNKNOWN-ACCESS (
CONSTANT (LOGICAL 4 0 0 LOGICAL ()) 0 0)) (422 'inactive' (DERIVED 423 0
0 DERIVED ()) () (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN
0 0 POINTER) UNKNOWN-ACCESS ()) (424 'ksp' (INTEGER 8 0 0 INTEGER ()) ()
(UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 POINTER)
UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0)) (425 'refcount' (
DERIVED 97 0 0 DERIVED ()) () (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0 POINTER) UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN
()) 0)) (426 'name' (CHARACTER 1 0 0 CHARACTER ((CONSTANT (INTEGER 4 0 0
INTEGER ()) 0 '101'))) () (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0) UNKNOWN-ACCESS (CONSTANT (CHARACTER 1 0 0 CHARACTER
(())) 0 101
'                                                                                                     ')))
PUBLIC (() () () ()) () 0 0 66371681)
50 'mpi_errcodes_ignore' 'mpi_interfaces' 'mpi_errcodes_ignore' 1 ((
VARIABLE UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 DIMENSION
IN_COMMON) (INTEGER 4 0 0 INTEGER ()) 0 0 () (1 0 EXPLICIT (CONSTANT (
INTEGER 4 0 0 INTEGER ()) 0 '1') (CONSTANT (INTEGER 4 0 0 INTEGER ()) 0
'1')) 0 () () () 0 0)
362 'csr_sparsity' 'sparse_tools' 'csr_sparsity' 1 ((DERIVED
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 POINTER_COMP) (UNKNOWN 0
0 0 UNKNOWN ()) 0 0 () () 0 ((427 'findrm' (INTEGER 4 0 0 INTEGER ()) (
1 0 DEFERRED () ()) (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN
UNKNOWN 0 0 DIMENSION POINTER) UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0
UNKNOWN ()) 0)) (428 'centrm' (INTEGER 4 0 0 INTEGER ()) (1 0 DEFERRED ()
()) (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DIMENSION POINTER) UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0)) (
429 'colm' (INTEGER 4 0 0 INTEGER ()) (1 0 DEFERRED () ()) (UNKNOWN-FL
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 DIMENSION POINTER)
UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0)) (430 'columns' (
INTEGER 4 0 0 INTEGER ()) () (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0) UNKNOWN-ACCESS ()) (431 'row_halo' (DERIVED 108 0 0
DERIVED ()) () (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0
0 POINTER) UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0)) (432
'column_halo' (DERIVED 108 0 0 DERIVED ()) () (UNKNOWN-FL UNKNOWN-INTENT
UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 POINTER) UNKNOWN-ACCESS (NULL (UNKNOWN
0 0 0 UNKNOWN ()) 0)) (433 'refcount' (DERIVED 97 0 0 DERIVED ()) () (
UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 POINTER)
UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0)) (434 'name' (
CHARACTER 1 0 0 CHARACTER ((CONSTANT (INTEGER 4 0 0 INTEGER ()) 0 '101')))
() (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0)
UNKNOWN-ACCESS (CONSTANT (CHARACTER 1 0 0 CHARACTER (())) 0 101
'                                                                                                     '))
(435 'wrapped' (LOGICAL 4 0 0 LOGICAL ()) () (UNKNOWN-FL UNKNOWN-INTENT
UNKNOWN-PROC UNKNOWN UNKNOWN 0 0) UNKNOWN-ACCESS (CONSTANT (LOGICAL 4 0
0 LOGICAL ()) 0 0)) (436 'sorted_rows' (LOGICAL 4 0 0 LOGICAL ()) () (
UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0)
UNKNOWN-ACCESS (CONSTANT (LOGICAL 4 0 0 LOGICAL ()) 0 0))) PUBLIC (() ()
() ()) () 0 0 78882345)
97 'refcount_type' 'reference_counting' 'refcount_type' 1 ((DERIVED
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 POINTER_COMP) (UNKNOWN 0
0 0 UNKNOWN ()) 0 0 () () 0 ((437 'prev' (DERIVED 97 0 0 DERIVED ()) ()
(UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 POINTER)
UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0)) (438 'next' (
DERIVED 97 0 0 DERIVED ()) () (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0 POINTER) UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN
()) 0)) (439 'count' (INTEGER 4 0 0 INTEGER ()) () (UNKNOWN-FL
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0) UNKNOWN-ACCESS (
CONSTANT (INTEGER 4 0 0 INTEGER ()) 0 '0')) (440 'id' (INTEGER 4 0 0
INTEGER ()) () (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0
0) UNKNOWN-ACCESS ()) (441 'name' (CHARACTER 1 0 0 CHARACTER ((CONSTANT
(INTEGER 4 0 0 INTEGER ()) 0 '101'))) () (UNKNOWN-FL UNKNOWN-INTENT
UNKNOWN-PROC UNKNOWN UNKNOWN 0 0) UNKNOWN-ACCESS ()) (442 'type' (
CHARACTER 1 0 0 CHARACTER ((CONSTANT (INTEGER 4 0 0 INTEGER ()) 0 '101')))
() (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0)
UNKNOWN-ACCESS ()) (443 'tagged' (LOGICAL 4 0 0 LOGICAL ()) () (
UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0)
UNKNOWN-ACCESS (CONSTANT (LOGICAL 4 0 0 LOGICAL ()) 0 0))) PUBLIC (() ()
() ()) () 0 0 25948645)
230 'name' '' 'name' 227 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DUMMY) (CHARACTER 1 0 0 CHARACTER (())) 0 0 () () 0 () () () 0 0)
231 'surface_field' '' 'surface_field' 227 ((VARIABLE UNKNOWN-INTENT
UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 POINTER RESULT ALWAYS_EXPLICIT) (
DERIVED 226 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
234 'field' '' 'field' 233 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0
0 DUMMY) (DERIVED 232 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
235 'n' '' 'n' 233 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 DUMMY)
(INTEGER 4 0 0 INTEGER ()) 0 0 () () 0 () () () 0 0)
236 'name' '' 'name' 233 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DUMMY) (CHARACTER 1 0 0 CHARACTER (())) 0 0 () () 0 () () () 0 0)
237 'surface_field' '' 'surface_field' 233 ((VARIABLE UNKNOWN-INTENT
UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 POINTER RESULT ALWAYS_EXPLICIT) (
DERIVED 232 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
249 'field' '' 'field' 248 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0
0 DUMMY) (DERIVED 232 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
250 'n' '' 'n' 248 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 DUMMY)
(INTEGER 4 0 0 INTEGER ()) 0 0 () () 0 () () () 0 0)
251 'name' '' 'name' 248 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DUMMY) (CHARACTER 1 0 0 CHARACTER (())) 0 0 () () 0 () () () 0 0)
252 'surface_field' '' 'surface_field' 248 ((VARIABLE UNKNOWN-INTENT
UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 POINTER RESULT ALWAYS_EXPLICIT) (
DERIVED 226 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
239 'field' '' 'field' 238 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0
0 DUMMY) (DERIVED 226 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
240 'bc_name' '' 'bc_name' 238 ((VARIABLE IN UNKNOWN-PROC UNKNOWN
UNKNOWN 0 0 DUMMY) (CHARACTER 1 0 0 CHARACTER (())) 0 0 () () 0 () () ()
0 0)
241 'name' '' 'name' 238 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DUMMY) (CHARACTER 1 0 0 CHARACTER (())) 0 0 () () 0 () () () 0 0)
242 'surface_field' '' 'surface_field' 238 ((VARIABLE UNKNOWN-INTENT
UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 POINTER RESULT ALWAYS_EXPLICIT) (
DERIVED 226 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
244 'field' '' 'field' 243 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0
0 DUMMY) (DERIVED 232 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
245 'bc_name' '' 'bc_name' 243 ((VARIABLE IN UNKNOWN-PROC UNKNOWN
UNKNOWN 0 0 DUMMY) (CHARACTER 1 0 0 CHARACTER (())) 0 0 () () 0 () () ()
0 0)
246 'name' '' 'name' 243 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DUMMY) (CHARACTER 1 0 0 CHARACTER (())) 0 0 () () 0 () () () 0 0)
247 'surface_field' '' 'surface_field' 243 ((VARIABLE UNKNOWN-INTENT
UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 POINTER RESULT ALWAYS_EXPLICIT) (
DERIVED 232 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
254 'field' '' 'field' 253 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0
0 DUMMY) (DERIVED 232 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
255 'bc_name' '' 'bc_name' 253 ((VARIABLE IN UNKNOWN-PROC UNKNOWN
UNKNOWN 0 0 DUMMY) (CHARACTER 1 0 0 CHARACTER (())) 0 0 () () 0 () () ()
0 0)
256 'name' '' 'name' 253 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DUMMY) (CHARACTER 1 0 0 CHARACTER (())) 0 0 () () 0 () () () 0 0)
257 'surface_field' '' 'surface_field' 253 ((VARIABLE UNKNOWN-INTENT
UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 POINTER RESULT ALWAYS_EXPLICIT) (
DERIVED 226 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
259 'field' '' 'field' 258 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0
0 DUMMY) (DERIVED 226 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
260 'n' '' 'n' 258 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 DUMMY)
(INTEGER 4 0 0 INTEGER ()) 0 0 () () 0 () () () 0 0)
261 'name' '' 'name' 258 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DUMMY) (CHARACTER 1 0 0 CHARACTER (())) 0 0 () () 0 () () () 0 0)
263 'field' '' 'field' 262 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0
0 DUMMY) (DERIVED 232 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
264 'n' '' 'n' 262 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 DUMMY)
(INTEGER 4 0 0 INTEGER ()) 0 0 () () 0 () () () 0 0)
265 'name' '' 'name' 262 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DUMMY) (CHARACTER 1 0 0 CHARACTER (())) 0 0 () () 0 () () () 0 0)
81 'field' '' 'field' 80 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DUMMY) (DERIVED 232 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
82 'n' '' 'n' 80 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 DUMMY) (
INTEGER 4 0 0 INTEGER ()) 0 0 () () 0 () () () 0 0)
83 'name' '' 'name' 80 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DUMMY) (CHARACTER 1 0 0 CHARACTER (())) 0 0 () () 0 () () () 0 0)
267 'field' '' 'field' 266 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0
0 DUMMY) (DERIVED 226 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
269 'field' '' 'field' 268 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0
0 DUMMY) (DERIVED 232 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
167 'field' '' 'field' 166 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0
0 DUMMY) (DERIVED 226 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
168 'n' '' 'n' 166 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 DUMMY)
(INTEGER 4 0 0 INTEGER ()) 0 0 () () 0 () () () 0 0)
169 'name' '' 'name' 166 ((VARIABLE OUT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
OPTIONAL DUMMY) (CHARACTER 1 0 0 CHARACTER (())) 0 0 () () 0 () () () 0
0)
170 'type' '' 'type' 166 ((VARIABLE OUT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
OPTIONAL DUMMY) (CHARACTER 1 0 0 CHARACTER (())) 0 0 () () 0 () () () 0
0)
171 'surface_element_list' '' 'surface_element_list' 166 ((VARIABLE
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 DIMENSION OPTIONAL
POINTER DUMMY) (INTEGER 4 0 0 INTEGER ()) 0 0 () (1 0 DEFERRED () ()) 0
() () () 0 0)
172 'surface_node_list' '' 'surface_node_list' 166 ((VARIABLE
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 DIMENSION OPTIONAL
POINTER DUMMY) (INTEGER 4 0 0 INTEGER ()) 0 0 () (1 0 DEFERRED () ()) 0
() () () 0 0)
173 'surface_mesh' '' 'surface_mesh' 166 ((VARIABLE UNKNOWN-INTENT
UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 OPTIONAL POINTER DUMMY) (DERIVED 72 0 0
DERIVED ()) 0 0 () () 0 () () () 0 0)
174 'option_path' '' 'option_path' 166 ((VARIABLE OUT UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0 OPTIONAL DUMMY) (CHARACTER 1 0 0 CHARACTER (())) 0 0
() () 0 () () () 0 0)
176 'field' '' 'field' 175 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0
0 DUMMY) (DERIVED 232 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
177 'n' '' 'n' 175 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 DUMMY)
(INTEGER 4 0 0 INTEGER ()) 0 0 () () 0 () () () 0 0)
178 'name' '' 'name' 175 ((VARIABLE OUT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
OPTIONAL DUMMY) (CHARACTER 1 0 0 CHARACTER (())) 0 0 () () 0 () () () 0
0)
179 'type' '' 'type' 175 ((VARIABLE OUT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
OPTIONAL DUMMY) (CHARACTER 1 0 0 CHARACTER (())) 0 0 () () 0 () () () 0
0)
180 'surface_element_list' '' 'surface_element_list' 175 ((VARIABLE
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 DIMENSION OPTIONAL
POINTER DUMMY) (INTEGER 4 0 0 INTEGER ()) 0 0 () (1 0 DEFERRED () ()) 0
() () () 0 0)
181 'surface_node_list' '' 'surface_node_list' 175 ((VARIABLE
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 DIMENSION OPTIONAL
POINTER DUMMY) (INTEGER 4 0 0 INTEGER ()) 0 0 () (1 0 DEFERRED () ()) 0
() () () 0 0)
182 'applies' '' 'applies' 175 ((VARIABLE OUT UNKNOWN-PROC UNKNOWN
UNKNOWN 0 0 DIMENSION OPTIONAL DUMMY) (LOGICAL 4 0 0 LOGICAL ()) 0 0 ()
(1 0 ASSUMED_SHAPE (CONSTANT (INTEGER 4 0 0 INTEGER ()) 0 '1') ()) 0 ()
() () 0 0)
183 'surface_mesh' '' 'surface_mesh' 175 ((VARIABLE UNKNOWN-INTENT
UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 OPTIONAL POINTER DUMMY) (DERIVED 72 0 0
DERIVED ()) 0 0 () () 0 () () () 0 0)
184 'option_path' '' 'option_path' 175 ((VARIABLE OUT UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0 OPTIONAL DUMMY) (CHARACTER 1 0 0 CHARACTER (())) 0 0
() () 0 () () () 0 0)
186 'field' '' 'field' 185 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0
0 DUMMY) (DERIVED 226 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
187 'name' '' 'name' 185 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DUMMY) (CHARACTER 1 0 0 CHARACTER (())) 0 0 () () 0 () () () 0 0)
188 'type' '' 'type' 185 ((VARIABLE OUT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
OPTIONAL DUMMY) (CHARACTER 1 0 0 CHARACTER (())) 0 0 () () 0 () () () 0
0)
189 'surface_node_list' '' 'surface_node_list' 185 ((VARIABLE
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 DIMENSION OPTIONAL
POINTER DUMMY) (INTEGER 4 0 0 INTEGER ()) 0 0 () (1 0 DEFERRED () ()) 0
() () () 0 0)
190 'surface_element_list' '' 'surface_element_list' 185 ((VARIABLE
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 DIMENSION OPTIONAL
POINTER DUMMY) (INTEGER 4 0 0 INTEGER ()) 0 0 () (1 0 DEFERRED () ()) 0
() () () 0 0)
191 'surface_mesh' '' 'surface_mesh' 185 ((VARIABLE UNKNOWN-INTENT
UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 OPTIONAL POINTER DUMMY) (DERIVED 72 0 0
DERIVED ()) 0 0 () () 0 () () () 0 0)
192 'option_path' '' 'option_path' 185 ((VARIABLE OUT UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0 OPTIONAL DUMMY) (CHARACTER 1 0 0 CHARACTER (())) 0 0
() () 0 () () () 0 0)
194 'field' '' 'field' 193 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0
0 DUMMY) (DERIVED 232 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
195 'name' '' 'name' 193 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DUMMY) (CHARACTER 1 0 0 CHARACTER (())) 0 0 () () 0 () () () 0 0)
196 'type' '' 'type' 193 ((VARIABLE OUT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
OPTIONAL DUMMY) (CHARACTER 1 0 0 CHARACTER (())) 0 0 () () 0 () () () 0
0)
294 'rhs' '' 'rhs' 291 ((VARIABLE INOUT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
OPTIONAL DUMMY) (DERIVED 226 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
292 'matrix' '' 'matrix' 291 ((VARIABLE INOUT UNKNOWN-PROC UNKNOWN
UNKNOWN 0 0 DUMMY) (DERIVED 367 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
295 'reference_value' '' 'reference_value' 291 ((VARIABLE IN
UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 OPTIONAL DUMMY) (REAL 8 0 0 REAL ()) 0
0 () () 0 () () () 0 0)
296 'reference_node_owned' '' 'reference_node_owned' 291 ((VARIABLE IN
UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 OPTIONAL DUMMY) (LOGICAL 4 0 0 LOGICAL
()) 0 0 () () 0 () () () 0 0)
298 'matrix' '' 'matrix' 297 ((VARIABLE INOUT UNKNOWN-PROC UNKNOWN
UNKNOWN 0 0 DUMMY) (DERIVED 444 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
299 'node' '' 'node' 297 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DUMMY) (INTEGER 4 0 0 INTEGER ()) 0 0 () () 0 () () () 0 0)
300 'rhs' '' 'rhs' 297 ((VARIABLE INOUT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
OPTIONAL DUMMY) (DERIVED 232 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
301 'reference_value' '' 'reference_value' 297 ((VARIABLE IN
UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 DIMENSION OPTIONAL DUMMY) (REAL 8 0 0
REAL ()) 0 0 () (1 0 EXPLICIT (CONSTANT (INTEGER 4 0 0 INTEGER ()) 0 '1')
(FUNCTION (INTEGER 4 0 0 INTEGER ()) 0 445 (('' (VARIABLE (DERIVED 444 0
0 DERIVED ()) 0 298 ())) ('' (CONSTANT (INTEGER 4 0 0 INTEGER ()) 0 '1')))
'petsc_csr_blocks_withdim' 1 446)) 0 () () () 0 0)
302 'reference_node_owned' '' 'reference_node_owned' 297 ((VARIABLE IN
UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 OPTIONAL DUMMY) (LOGICAL 4 0 0 LOGICAL
()) 0 0 () () 0 () () () 0 0)
304 'field' '' 'field' 303 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0
0 DUMMY) (DERIVED 226 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
305 'type' '' 'type' 303 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DUMMY) (CHARACTER 1 0 0 CHARACTER (())) 0 0 () () 0 () () () 0 0)
307 'field' '' 'field' 306 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0
0 DUMMY) (DERIVED 232 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
308 'type' '' 'type' 306 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DUMMY) (CHARACTER 1 0 0 CHARACTER (())) 0 0 () () 0 () () () 0 0)
310 'field' '' 'field' 309 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0
0 DUMMY) (DERIVED 226 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
311 'name' '' 'name' 309 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DUMMY) (CHARACTER 1 0 0 CHARACTER (())) 0 0 () () 0 () () () 0 0)
313 'field' '' 'field' 312 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0
0 DUMMY) (DERIVED 232 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
314 'name' '' 'name' 312 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DUMMY) (CHARACTER 1 0 0 CHARACTER (())) 0 0 () () 0 () () () 0 0)
85 'states' '' 'states' 84 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0
0 DIMENSION DUMMY) (DERIVED 447 0 0 DERIVED ()) 0 0 () (1 0
ASSUMED_SHAPE (CONSTANT (INTEGER 4 0 0 INTEGER ()) 0 '1') ()) 0 () () ()
0 0)
316 'field' '' 'field' 315 ((VARIABLE INOUT UNKNOWN-PROC UNKNOWN UNKNOWN
0 0 DUMMY) (DERIVED 226 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
318 'field' '' 'field' 317 ((VARIABLE INOUT UNKNOWN-PROC UNKNOWN UNKNOWN
0 0 DUMMY) (DERIVED 232 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
320 'matrix' '' 'matrix' 319 ((VARIABLE INOUT UNKNOWN-PROC UNKNOWN
UNKNOWN 0 0 DUMMY) (DERIVED 367 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
321 'rhs' '' 'rhs' 319 ((VARIABLE INOUT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DUMMY) (DERIVED 226 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
322 'field' '' 'field' 319 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0
0 DUMMY) (DERIVED 226 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
323 'dt' '' 'dt' 319 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
OPTIONAL DUMMY) (REAL 8 0 0 REAL ()) 0 0 () () 0 () () () 0 0)
325 'lhs' '' 'lhs' 324 ((VARIABLE INOUT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DUMMY) (DERIVED 226 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
326 'rhs' '' 'rhs' 324 ((VARIABLE INOUT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DUMMY) (DERIVED 226 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
327 'field' '' 'field' 324 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0
0 DUMMY) (DERIVED 226 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
328 'dt' '' 'dt' 324 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
OPTIONAL DUMMY) (REAL 8 0 0 REAL ()) 0 0 () () 0 () () () 0 0)
330 'matrix' '' 'matrix' 329 ((VARIABLE INOUT UNKNOWN-PROC UNKNOWN
UNKNOWN 0 0 DUMMY) (DERIVED 444 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
331 'rhs' '' 'rhs' 329 ((VARIABLE INOUT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DIMENSION OPTIONAL DUMMY) (DERIVED 448 0 0 DERIVED ()) 0 0 () (1 0
ASSUMED_SHAPE (CONSTANT (INTEGER 4 0 0 INTEGER ()) 0 '1') ()) 0 () () ()
0 0)
332 'sfields' '' 'sfields' 329 ((VARIABLE IN UNKNOWN-PROC UNKNOWN
UNKNOWN 0 0 DIMENSION DUMMY) (DERIVED 448 0 0 DERIVED ()) 0 0 () (1 0
ASSUMED_SHAPE (CONSTANT (INTEGER 4 0 0 INTEGER ()) 0 '1') ()) 0 () () ()
0 0)
333 'dt' '' 'dt' 329 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
OPTIONAL DUMMY) (REAL 8 0 0 REAL ()) 0 0 () () 0 () () () 0 0)
335 'matrix' '' 'matrix' 334 ((VARIABLE INOUT UNKNOWN-PROC UNKNOWN
UNKNOWN 0 0 DUMMY) (DERIVED 449 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
336 'rhs' '' 'rhs' 334 ((VARIABLE INOUT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
OPTIONAL DUMMY) (DERIVED 232 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
337 'field' '' 'field' 334 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0
0 DUMMY) (DERIVED 232 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
338 'dt' '' 'dt' 334 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
OPTIONAL DUMMY) (REAL 8 0 0 REAL ()) 0 0 () () 0 () () () 0 0)
340 'matrix' '' 'matrix' 339 ((VARIABLE INOUT UNKNOWN-PROC UNKNOWN
UNKNOWN 0 0 DUMMY) (DERIVED 444 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
341 'rhs' '' 'rhs' 339 ((VARIABLE INOUT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
OPTIONAL DUMMY) (DERIVED 232 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
342 'field' '' 'field' 339 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0
0 DUMMY) (DERIVED 232 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
343 'dt' '' 'dt' 339 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
OPTIONAL DUMMY) (REAL 8 0 0 REAL ()) 0 0 () () 0 () () () 0 0)
345 'matrix' '' 'matrix' 344 ((VARIABLE INOUT UNKNOWN-PROC UNKNOWN
UNKNOWN 0 0 DUMMY) (DERIVED 367 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
346 'rhs' '' 'rhs' 344 ((VARIABLE INOUT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
OPTIONAL DUMMY) (DERIVED 232 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
347 'field' '' 'field' 344 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0
0 DUMMY) (DERIVED 232 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
348 'dt' '' 'dt' 344 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
OPTIONAL DUMMY) (REAL 8 0 0 REAL ()) 0 0 () () 0 () () () 0 0)
349 'dim' '' 'dim' 344 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DUMMY) (INTEGER 4 0 0 INTEGER ()) 0 0 () () 0 () () () 0 0)
353 'field' '' 'field' 350 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0
0 DUMMY) (DERIVED 232 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
352 'rhs' '' 'rhs' 350 ((VARIABLE INOUT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
OPTIONAL DUMMY) (DERIVED 232 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
351 'lhs' '' 'lhs' 350 ((VARIABLE INOUT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DUMMY) (DERIVED 232 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
354 'dt' '' 'dt' 350 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
OPTIONAL DUMMY) (REAL 8 0 0 REAL ()) 0 0 () () 0 () () () 0 0)
356 'lhs' '' 'lhs' 355 ((VARIABLE INOUT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DUMMY) (DERIVED 226 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
357 'rhs' '' 'rhs' 355 ((VARIABLE INOUT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
OPTIONAL DUMMY) (DERIVED 232 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
358 'field' '' 'field' 355 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0
0 DUMMY) (DERIVED 232 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
359 'dt' '' 'dt' 355 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
OPTIONAL DUMMY) (REAL 8 0 0 REAL ()) 0 0 () () 0 () () () 0 0)
360 'dim' '' 'dim' 355 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DUMMY) (INTEGER 4 0 0 INTEGER ()) 0 0 () () 0 () () () 0 0)
68 'input_states' '' 'input_states' 67 ((VARIABLE IN UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0 DIMENSION DUMMY) (DERIVED 447 0 0 DERIVED ()) 0 0 ()
(1 0 ASSUMED_SHAPE (CONSTANT (INTEGER 4 0 0 INTEGER ()) 0 '1') ()) 0 ()
() () 0 0)
69 'collapsed_states' '' 'collapsed_states' 67 ((VARIABLE INOUT
UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 DIMENSION DUMMY) (DERIVED 447 0 0
DERIVED ()) 0 0 () (1 0 EXPLICIT (CONSTANT (INTEGER 4 0 0 INTEGER ()) 0
'1') (FUNCTION (INTEGER 4 0 0 INTEGER ()) 0 450 (('' (VARIABLE (DERIVED
447 0 0 DERIVED ()) 1 68 ((ARRAY (FULL 0))))) ('' ()) ('' ())) '' 0 'size'))
0 () () () 0 0)
70 'bctype' '' 'bctype' 67 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0
0 OPTIONAL DUMMY) (CHARACTER 1 0 0 CHARACTER (())) 0 0 () () 0 () () ()
0 0)
444 'petsc_csr_matrix' 'sparse_tools_petsc' 'petsc_csr_matrix' 1 ((
DERIVED UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 POINTER_COMP) (
UNKNOWN 0 0 0 UNKNOWN ()) 0 0 () () 0 ((451 'm' (INTEGER 8 0 0 INTEGER ())
() (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0)
UNKNOWN-ACCESS ()) (452 'row_numbering' (DERIVED 453 0 0 DERIVED ()) ()
(UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0)
UNKNOWN-ACCESS (STRUCTURE (DERIVED 453 0 0 DERIVED ()) 0 (((NULL (
UNKNOWN 0 0 0 UNKNOWN ()) 0) ()) (() ()) (() ()) (() ()) (() ()) ((NULL
(UNKNOWN 0 0 0 UNKNOWN ()) 0) ()) ((NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0) ())
((NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0) ()) ((CONSTANT (CHARACTER 1 0 0
CHARACTER (())) 0 101
'                                                                                                     ')
())) ())) (454 'column_numbering' (DERIVED 453 0 0 DERIVED ()) () (
UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0)
UNKNOWN-ACCESS (STRUCTURE (DERIVED 453 0 0 DERIVED ()) 0 (((NULL (
UNKNOWN 0 0 0 UNKNOWN ()) 0) ()) (() ()) (() ()) (() ()) (() ()) ((NULL
(UNKNOWN 0 0 0 UNKNOWN ()) 0) ()) ((NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0) ())
((NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0) ()) ((CONSTANT (CHARACTER 1 0 0
CHARACTER (())) 0 101
'                                                                                                     ')
())) ())) (455 'row_halo' (DERIVED 108 0 0 DERIVED ()) () (UNKNOWN-FL
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 POINTER) UNKNOWN-ACCESS
(NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0)) (456 'column_halo' (DERIVED 108 0 0
DERIVED ()) () (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0
0 POINTER) UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0)) (457
'refcount' (DERIVED 97 0 0 DERIVED ()) () (UNKNOWN-FL UNKNOWN-INTENT
UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 POINTER) UNKNOWN-ACCESS (NULL (UNKNOWN
0 0 0 UNKNOWN ()) 0)) (458 'name' (CHARACTER 1 0 0 CHARACTER ((CONSTANT
(INTEGER 4 0 0 INTEGER ()) 0 '101'))) () (UNKNOWN-FL UNKNOWN-INTENT
UNKNOWN-PROC UNKNOWN UNKNOWN 0 0) UNKNOWN-ACCESS (CONSTANT (CHARACTER 1
0 0 CHARACTER (())) 0 101
'                                                                                                     '))
(459 'is_assembled' (LOGICAL 4 0 0 LOGICAL ()) () (UNKNOWN-FL
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0) UNKNOWN-ACCESS (
CONSTANT (LOGICAL 4 0 0 LOGICAL ()) 0 0))) PUBLIC (() () () ()) () 0 0
75554753)
383 'scalar_boundary_conditions_ptr' 'fields_data_types'
'scalar_boundary_conditions_ptr' 1 ((DERIVED UNKNOWN-INTENT UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0 POINTER_COMP) (UNKNOWN 0 0 0 UNKNOWN ()) 0 0 () () 0
((460 'boundary_condition' (DERIVED 461 0 0 DERIVED ()) (1 0 DEFERRED ()
()) (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DIMENSION POINTER) UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0)))
PUBLIC (() () () ()) () 0 0 47502942)
448 'scalar_field_pointer' 'fields_data_types' 'scalar_field_pointer' 1
((DERIVED UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 POINTER_COMP)
(UNKNOWN 0 0 0 UNKNOWN ()) 0 0 () () 0 ((462 'ptr' (DERIVED 226 0 0
DERIVED ()) () (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0
0 POINTER) UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0))) PUBLIC (
() () () ()) () 0 0 54828570)
399 'vector_boundary_conditions_ptr' 'fields_data_types'
'vector_boundary_conditions_ptr' 1 ((DERIVED UNKNOWN-INTENT UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0 POINTER_COMP) (UNKNOWN 0 0 0 UNKNOWN ()) 0 0 () () 0
((463 'boundary_condition' (DERIVED 464 0 0 DERIVED ()) (1 0 DEFERRED ()
()) (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DIMENSION POINTER) UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0)))
PUBLIC (() () () ()) () 0 0 84476181)
55 'petsc_comm_self' 'petscsysdef' 'petsc_comm_self' 1 ((VARIABLE
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 IN_COMMON) (INTEGER 4 0
0 INTEGER ()) 0 0 () () 0 () () () 0 0)
423 'logical_array_ptr' 'sparse_tools' 'logical_array_ptr' 1 ((DERIVED
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 POINTER_COMP) (UNKNOWN 0
0 0 UNKNOWN ()) 0 0 () () 0 ((465 'ptr' (LOGICAL 4 0 0 LOGICAL ()) (1 0
DEFERRED () ()) (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN
0 0 DIMENSION POINTER) UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0)))
PUBLIC (() () () ()) () 0 0 30974511)
108 'halo_type' 'halo_data_types' 'halo_type' 1 ((DERIVED UNKNOWN-INTENT
UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 POINTER_COMP) (UNKNOWN 0 0 0 UNKNOWN ())
0 0 () () 0 ((466 'name' (CHARACTER 1 0 0 CHARACTER ((CONSTANT (INTEGER
4 0 0 INTEGER ()) 0 '101'))) () (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0) UNKNOWN-ACCESS ()) (467 'refcount' (DERIVED 97 0 0
DERIVED ()) () (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0
0 POINTER) UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0)) (468
'data_type' (INTEGER 4 0 0 INTEGER ()) () (UNKNOWN-FL UNKNOWN-INTENT
UNKNOWN-PROC UNKNOWN UNKNOWN 0 0) UNKNOWN-ACCESS (CONSTANT (INTEGER 4 0
0 INTEGER ()) 0 '0')) (469 'ordering_scheme' (INTEGER 4 0 0 INTEGER ())
() (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0)
UNKNOWN-ACCESS (CONSTANT (INTEGER 4 0 0 INTEGER ()) 0 '0')) (470
'communicator' (INTEGER 4 0 0 INTEGER ()) () (UNKNOWN-FL UNKNOWN-INTENT
UNKNOWN-PROC UNKNOWN UNKNOWN 0 0) UNKNOWN-ACCESS ()) (471 'nprocs' (
INTEGER 4 0 0 INTEGER ()) () (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0) UNKNOWN-ACCESS (CONSTANT (INTEGER 4 0 0 INTEGER ())
0 '0')) (472 'sends' (DERIVED 473 0 0 DERIVED ()) (1 0 DEFERRED () ()) (
UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 DIMENSION
POINTER) UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0)) (474
'receives' (DERIVED 473 0 0 DERIVED ()) (1 0 DEFERRED () ()) (
UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 DIMENSION
POINTER) UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0)) (475
'nowned_nodes' (INTEGER 4 0 0 INTEGER ()) () (UNKNOWN-FL UNKNOWN-INTENT
UNKNOWN-PROC UNKNOWN UNKNOWN 0 0) UNKNOWN-ACCESS (CONSTANT (INTEGER 4 0
0 INTEGER ()) 0 '-1')) (476 'owners' (INTEGER 4 0 0 INTEGER ()) (1 0
DEFERRED () ()) (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN
0 0 DIMENSION POINTER) UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0))
(477 'unn_count' (INTEGER 4 0 0 INTEGER ()) () (UNKNOWN-FL
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0) UNKNOWN-ACCESS (
CONSTANT (INTEGER 4 0 0 INTEGER ()) 0 '-1')) (478 'owned_nodes_unn_base'
(INTEGER 4 0 0 INTEGER ()) (1 0 DEFERRED () ()) (UNKNOWN-FL
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 DIMENSION POINTER)
UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0)) (479
'my_owned_nodes_unn_base' (INTEGER 4 0 0 INTEGER ()) () (UNKNOWN-FL
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0) UNKNOWN-ACCESS (
CONSTANT (INTEGER 4 0 0 INTEGER ()) 0 '-1')) (480 'receives_gnn_to_unn'
(INTEGER 4 0 0 INTEGER ()) (1 0 DEFERRED () ()) (UNKNOWN-FL
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 DIMENSION POINTER)
UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0)) (481 'gnn_to_unn' (
INTEGER 4 0 0 INTEGER ()) (1 0 DEFERRED () ()) (UNKNOWN-FL
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 DIMENSION POINTER)
UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0))) PUBLIC (() () () ())
() 0 0 63994053)
111 'integer_set_vector' 'integer_set_module' 'integer_set_vector' 1 ((
DERIVED UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 POINTER_COMP) (
UNKNOWN 0 0 0 UNKNOWN ()) 0 0 () () 0 ((482 'sets' (DERIVED 483 0 0
DERIVED ()) (1 0 DEFERRED () ()) (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0 DIMENSION POINTER) UNKNOWN-ACCESS ())) PUBLIC (() ()
() ()) () 0 0 40688950)
122 'polynomial' 'polynomials' 'polynomial' 1 ((DERIVED UNKNOWN-INTENT
UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 POINTER_COMP) (UNKNOWN 0 0 0 UNKNOWN ())
0 0 () () 0 ((484 'coefs' (REAL 8 0 0 REAL ()) (1 0 DEFERRED () ()) (
UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 DIMENSION
POINTER) UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0)) (485 'degree'
(INTEGER 4 0 0 INTEGER ()) () (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0) UNKNOWN-ACCESS (CONSTANT (INTEGER 4 0 0 INTEGER ())
0 '-1'))) PUBLIC (() () () ()) () 0 0 87989236)
125 'ele_numbering_type' 'element_numbering' 'ele_numbering_type' 1 ((
DERIVED UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 POINTER_COMP) (
UNKNOWN 0 0 0 UNKNOWN ()) 0 0 () () 0 ((486 'faces' (INTEGER 4 0 0
INTEGER ()) () (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0
0) UNKNOWN-ACCESS ()) (487 'vertices' (INTEGER 4 0 0 INTEGER ()) () (
UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0)
UNKNOWN-ACCESS ()) (488 'edges' (INTEGER 4 0 0 INTEGER ()) () (
UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0)
UNKNOWN-ACCESS ()) (489 'boundaries' (INTEGER 4 0 0 INTEGER ()) () (
UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0)
UNKNOWN-ACCESS ()) (490 'degree' (INTEGER 4 0 0 INTEGER ()) () (
UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0)
UNKNOWN-ACCESS ()) (491 'dimension' (INTEGER 4 0 0 INTEGER ()) () (
UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0)
UNKNOWN-ACCESS ()) (492 'nodes' (INTEGER 4 0 0 INTEGER ()) () (
UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0)
UNKNOWN-ACCESS ()) (493 'type' (INTEGER 4 0 0 INTEGER ()) () (
UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0)
UNKNOWN-ACCESS (CONSTANT (INTEGER 4 0 0 INTEGER ()) 0 '1')) (494 'family'
(INTEGER 4 0 0 INTEGER ()) () (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0) UNKNOWN-ACCESS ()) (495 'count2number' (INTEGER 4 0
0 INTEGER ()) (3 0 DEFERRED () () () () () ()) (UNKNOWN-FL
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 DIMENSION POINTER)
UNKNOWN-ACCESS ()) (496 'number2count' (INTEGER 4 0 0 INTEGER ()) (2 0
DEFERRED () () () ()) (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN
UNKNOWN 0 0 DIMENSION POINTER) UNKNOWN-ACCESS ()) (497 'boundary_coord'
(INTEGER 4 0 0 INTEGER ()) (1 0 DEFERRED () ()) (UNKNOWN-FL
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 DIMENSION POINTER)
UNKNOWN-ACCESS ()) (498 'boundary_val' (INTEGER 4 0 0 INTEGER ()) (1 0
DEFERRED () ()) (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN
0 0 DIMENSION POINTER) UNKNOWN-ACCESS ())) PUBLIC (() () () ()) () 0 0
96431082)
129 'superconvergence_type' 'elements' 'superconvergence_type' 1 ((
DERIVED UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 POINTER_COMP) (
UNKNOWN 0 0 0 UNKNOWN ()) 0 0 () () 0 ((499 'nsp' (INTEGER 4 0 0 INTEGER
()) () (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0)
UNKNOWN-ACCESS ()) (500 'l' (REAL 8 0 0 REAL ()) (2 0 DEFERRED () () ()
()) (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DIMENSION POINTER) UNKNOWN-ACCESS ()) (501 'n' (REAL 8 0 0 REAL ()) (2 0
DEFERRED () () () ()) (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN
UNKNOWN 0 0 DIMENSION POINTER) UNKNOWN-ACCESS ()) (502 'dn' (REAL 8 0 0
REAL ()) (3 0 DEFERRED () () () () () ()) (UNKNOWN-FL UNKNOWN-INTENT
UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 DIMENSION POINTER) UNKNOWN-ACCESS ()))
PUBLIC (() () () ()) () 0 0 18282395)
131 'constraints_type' 'elements' 'constraints_type' 1 ((DERIVED
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 POINTER_COMP) (UNKNOWN 0
0 0 UNKNOWN ()) 0 0 () () 0 ((503 'type' (INTEGER 4 0 0 INTEGER ()) () (
UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0)
UNKNOWN-ACCESS ()) (504 'dim' (INTEGER 4 0 0 INTEGER ()) () (UNKNOWN-FL
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0) UNKNOWN-ACCESS ()) (
505 'degree' (INTEGER 4 0 0 INTEGER ()) () (UNKNOWN-FL UNKNOWN-INTENT
UNKNOWN-PROC UNKNOWN UNKNOWN 0 0) UNKNOWN-ACCESS ()) (506 'loc' (
INTEGER 4 0 0 INTEGER ()) () (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0) UNKNOWN-ACCESS ()) (507 'n_constraints' (INTEGER 4
0 0 INTEGER ()) () (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN
UNKNOWN 0 0) UNKNOWN-ACCESS ()) (508 'orthogonal' (REAL 8 0 0 REAL ()) (
3 0 DEFERRED () () () () () ()) (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0 DIMENSION POINTER) UNKNOWN-ACCESS (NULL (UNKNOWN 0 0
0 UNKNOWN ()) 0))) PUBLIC (() () () ()) () 0 0 57548971)
407 'picker_ptr' 'picker_data_types' 'picker_ptr' 1 ((DERIVED
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 POINTER_COMP) (UNKNOWN 0
0 0 UNKNOWN ()) 0 0 () () 0 ((509 'ptr' (DERIVED 510 0 0 DERIVED ()) ()
(UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 POINTER)
UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0))) PUBLIC (() () () ())
() 0 0 71120007)
510 'picker_type' 'picker_data_types' 'picker_type' 1 ((DERIVED
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 POINTER_COMP) (UNKNOWN 0
0 0 UNKNOWN ()) 0 0 () () 0 ((511 'name' (CHARACTER 1 0 0 CHARACTER ((
CONSTANT (INTEGER 4 0 0 INTEGER ()) 0 '101'))) () (UNKNOWN-FL
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0) UNKNOWN-ACCESS ()) (
512 'refcount' (DERIVED 97 0 0 DERIVED ()) () (UNKNOWN-FL UNKNOWN-INTENT
UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 POINTER) UNKNOWN-ACCESS (NULL (UNKNOWN
0 0 0 UNKNOWN ()) 0)) (513 'picker_id' (INTEGER 4 0 0 INTEGER ()) () (
UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0)
UNKNOWN-ACCESS (CONSTANT (INTEGER 4 0 0 INTEGER ()) 0 '0')) (514
'last_mesh_movement' (INTEGER 4 0 0 INTEGER ()) () (UNKNOWN-FL
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0) UNKNOWN-ACCESS (
CONSTANT (INTEGER 4 0 0 INTEGER ()) 0 '0'))) PUBLIC (() () () ()) () 0 0
8821665)
446 'petsc_csr_blocks_withdim' 'sparse_tools_petsc'
'petsc_csr_blocks_withdim' 1 ((PROCEDURE UNKNOWN-INTENT MODULE-PROC DECL
UNKNOWN 0 0 FUNCTION PURE ALWAYS_EXPLICIT) (INTEGER 4 0 0 INTEGER ())
515 0 (516 517) () 518 () () () 0 0)
516 'matrix' '' 'matrix' 515 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN
0 0 DUMMY) (DERIVED 444 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
453 'petsc_numbering_type' 'petsc_tools' 'petsc_numbering_type' 1 ((
DERIVED UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 POINTER_COMP) (
UNKNOWN 0 0 0 UNKNOWN ()) 0 0 () () 0 ((519 'halo' (DERIVED 108 0 0
DERIVED ()) () (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0
0 POINTER) UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0)) (520
'nprivatenodes' (INTEGER 4 0 0 INTEGER ()) () (UNKNOWN-FL UNKNOWN-INTENT
UNKNOWN-PROC UNKNOWN UNKNOWN 0 0) UNKNOWN-ACCESS ()) (521
'universal_length' (INTEGER 4 0 0 INTEGER ()) () (UNKNOWN-FL
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0) UNKNOWN-ACCESS ()) (
522 'offset' (INTEGER 4 0 0 INTEGER ()) () (UNKNOWN-FL UNKNOWN-INTENT
UNKNOWN-PROC UNKNOWN UNKNOWN 0 0) UNKNOWN-ACCESS ()) (523 'gnn2unn' (
INTEGER 4 0 0 INTEGER ()) (2 0 DEFERRED () () () ()) (UNKNOWN-FL
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 DIMENSION POINTER)
UNKNOWN-ACCESS ()) (524 'ghost_nodes' (INTEGER 4 0 0 INTEGER ()) (1 0
DEFERRED () ()) (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN
0 0 DIMENSION POINTER) UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0))
(525 'ghost2unn' (INTEGER 4 0 0 INTEGER ()) (2 0 DEFERRED () () () ()) (
UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 DIMENSION
POINTER) UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0)) (526
'refcount' (DERIVED 97 0 0 DERIVED ()) () (UNKNOWN-FL UNKNOWN-INTENT
UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 POINTER) UNKNOWN-ACCESS (NULL (UNKNOWN
0 0 0 UNKNOWN ()) 0)) (527 'name' (CHARACTER 1 0 0 CHARACTER ((CONSTANT
(INTEGER 4 0 0 INTEGER ()) 0 '101'))) () (UNKNOWN-FL UNKNOWN-INTENT
UNKNOWN-PROC UNKNOWN UNKNOWN 0 0) UNKNOWN-ACCESS (CONSTANT (CHARACTER 1
0 0 CHARACTER (())) 0 101
'                                                                                                     ')))
PUBLIC (() () () ()) () 0 0 16529124)
447 'state_type' 'state_module' 'state_type' 1 ((DERIVED UNKNOWN-INTENT
UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 POINTER_COMP) (UNKNOWN 0 0 0 UNKNOWN ())
0 0 () () 0 ((528 'name' (CHARACTER 1 0 0 CHARACTER ((CONSTANT (INTEGER
4 0 0 INTEGER ()) 0 '101'))) () (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0) UNKNOWN-ACCESS (CONSTANT (CHARACTER 1 0 0 CHARACTER
(())) 0 101
'                                                                                                     '))
(529 'option_path' (CHARACTER 1 0 0 CHARACTER ((CONSTANT (INTEGER 4 0 0
INTEGER ()) 0 '8192'))) () (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0) UNKNOWN-ACCESS (CONSTANT (CHARACTER 1 0 0 CHARACTER
(())) 0 8192
'/uninitialised_path/                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            '))
(530 'vector_names' (CHARACTER 1 0 0 CHARACTER ((CONSTANT (INTEGER 4 0 0
INTEGER ()) 0 '101'))) (1 0 DEFERRED () ()) (UNKNOWN-FL UNKNOWN-INTENT
UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 DIMENSION POINTER) UNKNOWN-ACCESS (
NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0)) (531 'scalar_names' (CHARACTER 1 0 0
CHARACTER ((CONSTANT (INTEGER 4 0 0 INTEGER ()) 0 '101'))) (1 0 DEFERRED
() ()) (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DIMENSION POINTER) UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0)) (
532 'mesh_names' (CHARACTER 1 0 0 CHARACTER ((CONSTANT (INTEGER 4 0 0
INTEGER ()) 0 '101'))) (1 0 DEFERRED () ()) (UNKNOWN-FL UNKNOWN-INTENT
UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 DIMENSION POINTER) UNKNOWN-ACCESS (
NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0)) (533 'halo_names' (CHARACTER 1 0 0
CHARACTER ((CONSTANT (INTEGER 4 0 0 INTEGER ()) 0 '101'))) (1 0 DEFERRED
() ()) (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DIMENSION POINTER) UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0)) (
534 'tensor_names' (CHARACTER 1 0 0 CHARACTER ((CONSTANT (INTEGER 4 0 0
INTEGER ()) 0 '101'))) (1 0 DEFERRED () ()) (UNKNOWN-FL UNKNOWN-INTENT
UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 DIMENSION POINTER) UNKNOWN-ACCESS (
NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0)) (535 'csr_sparsity_names' (
CHARACTER 1 0 0 CHARACTER ((CONSTANT (INTEGER 4 0 0 INTEGER ()) 0 '101')))
(1 0 DEFERRED () ()) (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN
UNKNOWN 0 0 DIMENSION POINTER) UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0
UNKNOWN ()) 0)) (536 'csr_matrix_names' (CHARACTER 1 0 0 CHARACTER ((
CONSTANT (INTEGER 4 0 0 INTEGER ()) 0 '101'))) (1 0 DEFERRED () ()) (
UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 DIMENSION
POINTER) UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0)) (537
'block_csr_matrix_names' (CHARACTER 1 0 0 CHARACTER ((CONSTANT (INTEGER
4 0 0 INTEGER ()) 0 '101'))) (1 0 DEFERRED () ()) (UNKNOWN-FL
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 DIMENSION POINTER)
UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0)) (538
'petsc_csr_matrix_names' (CHARACTER 1 0 0 CHARACTER ((CONSTANT (INTEGER
4 0 0 INTEGER ()) 0 '101'))) (1 0 DEFERRED () ()) (UNKNOWN-FL
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 DIMENSION POINTER)
UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0)) (539 'vector_fields'
(DERIVED 540 0 0 DERIVED ()) (1 0 DEFERRED () ()) (UNKNOWN-FL
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 DIMENSION POINTER)
UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0)) (541 'tensor_fields'
(DERIVED 542 0 0 DERIVED ()) (1 0 DEFERRED () ()) (UNKNOWN-FL
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 DIMENSION POINTER)
UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0)) (543 'scalar_fields'
(DERIVED 448 0 0 DERIVED ()) (1 0 DEFERRED () ()) (UNKNOWN-FL
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 DIMENSION POINTER)
UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0)) (544 'meshes' (
DERIVED 545 0 0 DERIVED ()) (1 0 DEFERRED () ()) (UNKNOWN-FL
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 DIMENSION POINTER)
UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0)) (546 'halos' (
DERIVED 547 0 0 DERIVED ()) (1 0 DEFERRED () ()) (UNKNOWN-FL
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 DIMENSION POINTER)
UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0)) (548 'csr_sparsities'
(DERIVED 549 0 0 DERIVED ()) (1 0 DEFERRED () ()) (UNKNOWN-FL
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 DIMENSION POINTER)
UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0)) (550 'csr_matrices'
(DERIVED 551 0 0 DERIVED ()) (1 0 DEFERRED () ()) (UNKNOWN-FL
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 DIMENSION POINTER)
UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0)) (552
'block_csr_matrices' (DERIVED 553 0 0 DERIVED ()) (1 0 DEFERRED () ()) (
UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 DIMENSION
POINTER) UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0)) (554
'petsc_csr_matrices' (DERIVED 555 0 0 DERIVED ()) (1 0 DEFERRED () ()) (
UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 DIMENSION
POINTER) UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0))) PUBLIC (()
() () ()) () 0 0 7575341)
449 'block_csr_matrix' 'sparse_tools' 'block_csr_matrix' 1 ((DERIVED
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 POINTER_COMP) (UNKNOWN 0
0 0 UNKNOWN ()) 0 0 () () 0 ((556 'sparsity' (DERIVED 362 0 0 DERIVED ())
() (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0)
UNKNOWN-ACCESS (STRUCTURE (DERIVED 362 0 0 DERIVED ()) 0 (((NULL (
UNKNOWN 0 0 0 UNKNOWN ()) 0) ()) ((NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0) ())
((NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0) ()) (() ()) ((NULL (UNKNOWN 0 0 0
UNKNOWN ()) 0) ()) ((NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0) ()) ((NULL (
UNKNOWN 0 0 0 UNKNOWN ()) 0) ()) ((CONSTANT (CHARACTER 1 0 0 CHARACTER (
())) 0 101
'                                                                                                     ')
()) ((CONSTANT (LOGICAL 4 0 0 LOGICAL ()) 0 0) ()) ((CONSTANT (LOGICAL 4
0 0 LOGICAL ()) 0 0) ())) ())) (557 'val' (DERIVED 558 0 0 DERIVED ()) (
2 0 DEFERRED () () () ()) (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0 DIMENSION POINTER) UNKNOWN-ACCESS (NULL (UNKNOWN 0 0
0 UNKNOWN ()) 0)) (559 'contiguous_val' (REAL 8 0 0 REAL ()) (1 0
DEFERRED () ()) (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN
0 0 DIMENSION POINTER) UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0))
(560 'ival' (DERIVED 473 0 0 DERIVED ()) (2 0 DEFERRED () () () ()) (
UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 DIMENSION
POINTER) UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0)) (561 'blocks'
(INTEGER 4 0 0 INTEGER ()) (1 0 EXPLICIT (CONSTANT (INTEGER 4 0 0
INTEGER ()) 0 '1') (CONSTANT (INTEGER 4 0 0 INTEGER ()) 0 '2')) (
UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 DIMENSION)
UNKNOWN-ACCESS (ARRAY (INTEGER 4 0 0 INTEGER ()) 1 (((CONSTANT (INTEGER
4 0 0 INTEGER ()) 0 '0') ()) ((CONSTANT (INTEGER 4 0 0 INTEGER ()) 0 '0')
())) ('2'))) (562 'clone' (LOGICAL 4 0 0 LOGICAL ()) () (UNKNOWN-FL
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0) UNKNOWN-ACCESS (
CONSTANT (LOGICAL 4 0 0 LOGICAL ()) 0 0)) (563 'external_val' (LOGICAL 4
0 0 LOGICAL ()) () (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN
UNKNOWN 0 0) UNKNOWN-ACCESS (CONSTANT (LOGICAL 4 0 0 LOGICAL ()) 0 0)) (
564 'columns' (INTEGER 4 0 0 INTEGER ()) () (UNKNOWN-FL UNKNOWN-INTENT
UNKNOWN-PROC UNKNOWN UNKNOWN 0 0) UNKNOWN-ACCESS ()) (565 'refcount' (
DERIVED 97 0 0 DERIVED ()) () (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0 POINTER) UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN
()) 0)) (566 'name' (CHARACTER 1 0 0 CHARACTER ((CONSTANT (INTEGER 4 0 0
INTEGER ()) 0 '101'))) () (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0) UNKNOWN-ACCESS (CONSTANT (CHARACTER 1 0 0 CHARACTER
(())) 0 101
'                                                                                                     '))
(567 'diagonal' (LOGICAL 4 0 0 LOGICAL ()) () (UNKNOWN-FL UNKNOWN-INTENT
UNKNOWN-PROC UNKNOWN UNKNOWN 0 0) UNKNOWN-ACCESS (CONSTANT (LOGICAL 4 0
0 LOGICAL ()) 0 0)) (568 'equal_diagonal_blocks' (LOGICAL 4 0 0 LOGICAL
()) () (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0)
UNKNOWN-ACCESS (CONSTANT (LOGICAL 4 0 0 LOGICAL ()) 0 0)) (569 'ksp' (
INTEGER 8 0 0 INTEGER ()) () (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0 POINTER) UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN
()) 0))) PUBLIC (() () () ()) () 0 0 61921235)
553 'block_csr_matrix_pointer' 'sparse_tools' 'block_csr_matrix_pointer'
1 ((DERIVED UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 POINTER_COMP)
(UNKNOWN 0 0 0 UNKNOWN ()) 0 0 () () 0 ((570 'ptr' (DERIVED 449 0 0
DERIVED ()) () (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0
0 POINTER) UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0))) PUBLIC (
() () () ()) () 0 0 83606449)
551 'csr_matrix_pointer' 'sparse_tools' 'csr_matrix_pointer' 1 ((
DERIVED UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 POINTER_COMP) (
UNKNOWN 0 0 0 UNKNOWN ()) 0 0 () () 0 ((571 'ptr' (DERIVED 367 0 0
DERIVED ()) () (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0
0 POINTER) UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0))) PUBLIC (
() () () ()) () 0 0 82770239)
549 'csr_sparsity_pointer' 'sparse_tools' 'csr_sparsity_pointer' 1 ((
DERIVED UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 POINTER_COMP) (
UNKNOWN 0 0 0 UNKNOWN ()) 0 0 () () 0 ((572 'ptr' (DERIVED 362 0 0
DERIVED ()) () (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0
0 POINTER) UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0))) PUBLIC (
() () () ()) () 0 0 9687815)
558 'real_vector' 'futils' 'real_vector' 1 ((DERIVED UNKNOWN-INTENT
UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 POINTER_COMP) (UNKNOWN 0 0 0 UNKNOWN ())
0 0 () () 0 ((573 'ptr' (REAL 8 0 0 REAL ()) (1 0 DEFERRED () ()) (
UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 DIMENSION
POINTER) UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0))) PUBLIC (()
() () ()) () 0 0 72870256)
445 'blocks' '' 'blocks' 297 ((PROCEDURE UNKNOWN-INTENT UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0 FUNCTION) (UNKNOWN 0 0 0 UNKNOWN ()) 0 0 () () 445 ()
() () 0 0)
450 'size' '(intrinsic)' 'size' 67 ((PROCEDURE UNKNOWN-INTENT
UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 FUNCTION) (UNKNOWN 0 0 0 UNKNOWN ()) 0
0 () () 450 () () () 0 0)
555 'petsc_csr_matrix_pointer' 'sparse_tools_petsc'
'petsc_csr_matrix_pointer' 1 ((DERIVED UNKNOWN-INTENT UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0 POINTER_COMP) (UNKNOWN 0 0 0 UNKNOWN ()) 0 0 () () 0
((574 'ptr' (DERIVED 444 0 0 DERIVED ()) () (UNKNOWN-FL UNKNOWN-INTENT
UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 POINTER) UNKNOWN-ACCESS (NULL (UNKNOWN
0 0 0 UNKNOWN ()) 0))) PUBLIC (() () () ()) () 0 0 85445535)
545 'mesh_pointer' 'fields_data_types' 'mesh_pointer' 1 ((DERIVED
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 POINTER_COMP) (UNKNOWN 0
0 0 UNKNOWN ()) 0 0 () () 0 ((575 'ptr' (DERIVED 72 0 0 DERIVED ()) () (
UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 POINTER)
UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0))) PUBLIC (() () () ())
() 0 0 65432384)
461 'scalar_boundary_condition' 'fields_data_types'
'scalar_boundary_condition' 1 ((DERIVED UNKNOWN-INTENT UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0 POINTER_COMP) (UNKNOWN 0 0 0 UNKNOWN ()) 0 0 () () 0
((576 'name' (CHARACTER 1 0 0 CHARACTER ((CONSTANT (INTEGER 4 0 0
INTEGER ()) 0 '101'))) () (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0) UNKNOWN-ACCESS ()) (577 'type' (CHARACTER 1 0 0
CHARACTER ((CONSTANT (INTEGER 4 0 0 INTEGER ()) 0 '101'))) () (
UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0)
UNKNOWN-ACCESS (CONSTANT (CHARACTER 1 0 0 CHARACTER (())) 0 101
'                                                                                                     '))
(578 'surface_element_list' (INTEGER 4 0 0 INTEGER ()) (1 0 DEFERRED ()
()) (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DIMENSION POINTER) UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0)) (
579 'surface_node_list' (INTEGER 4 0 0 INTEGER ()) (1 0 DEFERRED () ())
(UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 DIMENSION
POINTER) UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0)) (580
'surface_mesh' (DERIVED 72 0 0 DERIVED ()) () (UNKNOWN-FL UNKNOWN-INTENT
UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 POINTER) UNKNOWN-ACCESS ()) (581
'surface_fields' (DERIVED 226 0 0 DERIVED ()) (1 0 DEFERRED () ()) (
UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 DIMENSION
POINTER) UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0)) (582
'option_path' (CHARACTER 1 0 0 CHARACTER ((CONSTANT (INTEGER 4 0 0
INTEGER ()) 0 '8192'))) () (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0) UNKNOWN-ACCESS (CONSTANT (CHARACTER 1 0 0 CHARACTER
(())) 0 8192
'/uninitialised_path/                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            ')))
PUBLIC (() () () ()) () 0 0 50740068)
542 'tensor_field_pointer' 'fields_data_types' 'tensor_field_pointer' 1
((DERIVED UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 POINTER_COMP)
(UNKNOWN 0 0 0 UNKNOWN ()) 0 0 () () 0 ((583 'ptr' (DERIVED 584 0 0
DERIVED ()) () (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0
0 POINTER) UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0))) PUBLIC (
() () () ()) () 0 0 16722823)
464 'vector_boundary_condition' 'fields_data_types'
'vector_boundary_condition' 1 ((DERIVED UNKNOWN-INTENT UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0 POINTER_COMP) (UNKNOWN 0 0 0 UNKNOWN ()) 0 0 () () 0
((585 'name' (CHARACTER 1 0 0 CHARACTER ((CONSTANT (INTEGER 4 0 0
INTEGER ()) 0 '101'))) () (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0) UNKNOWN-ACCESS ()) (586 'type' (CHARACTER 1 0 0
CHARACTER ((CONSTANT (INTEGER 4 0 0 INTEGER ()) 0 '101'))) () (
UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0)
UNKNOWN-ACCESS (CONSTANT (CHARACTER 1 0 0 CHARACTER (())) 0 101
'                                                                                                     '))
(587 'applies' (LOGICAL 4 0 0 LOGICAL ()) (1 0 EXPLICIT (CONSTANT (
INTEGER 4 0 0 INTEGER ()) 0 '1') (CONSTANT (INTEGER 4 0 0 INTEGER ()) 0
'3')) (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DIMENSION) UNKNOWN-ACCESS ()) (588 'surface_element_list' (INTEGER 4 0 0
INTEGER ()) (1 0 DEFERRED () ()) (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0 DIMENSION POINTER) UNKNOWN-ACCESS (NULL (UNKNOWN 0 0
0 UNKNOWN ()) 0)) (589 'surface_node_list' (INTEGER 4 0 0 INTEGER ()) (
1 0 DEFERRED () ()) (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN
UNKNOWN 0 0 DIMENSION POINTER) UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0
UNKNOWN ()) 0)) (590 'surface_mesh' (DERIVED 72 0 0 DERIVED ()) () (
UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 POINTER)
UNKNOWN-ACCESS ()) (591 'surface_fields' (DERIVED 232 0 0 DERIVED ()) (
1 0 DEFERRED () ()) (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN
UNKNOWN 0 0 DIMENSION POINTER) UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0
UNKNOWN ()) 0)) (592 'scalar_surface_fields' (DERIVED 226 0 0 DERIVED ())
(1 0 DEFERRED () ()) (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN
UNKNOWN 0 0 DIMENSION POINTER) UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0
UNKNOWN ()) 0)) (593 'option_path' (CHARACTER 1 0 0 CHARACTER ((
CONSTANT (INTEGER 4 0 0 INTEGER ()) 0 '8192'))) () (UNKNOWN-FL
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0) UNKNOWN-ACCESS (
CONSTANT (CHARACTER 1 0 0 CHARACTER (())) 0 8192
'/uninitialised_path/                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            ')))
PUBLIC (() () () ()) () 0 0 42205165)
540 'vector_field_pointer' 'fields_data_types' 'vector_field_pointer' 1
((DERIVED UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 POINTER_COMP)
(UNKNOWN 0 0 0 UNKNOWN ()) 0 0 () () 0 ((594 'ptr' (DERIVED 232 0 0
DERIVED ()) () (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0
0 POINTER) UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0))) PUBLIC (
() () () ()) () 0 0 74448721)
518 'blocks' '' 'blocks' 515 ((VARIABLE UNKNOWN-INTENT UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0 RESULT ALWAYS_EXPLICIT) (INTEGER 4 0 0 INTEGER ()) 0
0 () () 0 () () () 0 0)
517 'dim' '' 'dim' 515 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
OPTIONAL DUMMY) (INTEGER 4 0 0 INTEGER ()) 0 0 () () 0 () () () 0 0)
547 'halo_pointer' 'halo_data_types' 'halo_pointer' 1 ((DERIVED
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 POINTER_COMP) (UNKNOWN 0
0 0 UNKNOWN ()) 0 0 () () 0 ((595 'ptr' (DERIVED 108 0 0 DERIVED ()) ()
(UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 POINTER)
UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0))) PUBLIC (() () () ())
() 0 0 72023282)
584 'tensor_field' 'fields_data_types' 'tensor_field' 1 ((DERIVED
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 POINTER_COMP) (UNKNOWN 0
0 0 UNKNOWN ()) 0 0 () () 0 ((596 'val' (REAL 8 0 0 REAL ()) (3 0
DEFERRED () () () () () ()) (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0 DIMENSION POINTER) UNKNOWN-ACCESS ()) (597 'wrapped'
(LOGICAL 4 0 0 LOGICAL ()) () (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0) UNKNOWN-ACCESS (CONSTANT (LOGICAL 4 0 0 LOGICAL ())
0 1)) (598 'field_type' (INTEGER 4 0 0 INTEGER ()) () (UNKNOWN-FL
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0) UNKNOWN-ACCESS (
CONSTANT (INTEGER 4 0 0 INTEGER ()) 0 '0')) (599 'name' (CHARACTER 1 0 0
CHARACTER ((CONSTANT (INTEGER 4 0 0 INTEGER ()) 0 '101'))) () (
UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0)
UNKNOWN-ACCESS ()) (600 'dim' (INTEGER 4 0 0 INTEGER ()) (1 0 EXPLICIT (
CONSTANT (INTEGER 4 0 0 INTEGER ()) 0 '1') (CONSTANT (INTEGER 4 0 0
INTEGER ()) 0 '2')) (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN
UNKNOWN 0 0 DIMENSION) UNKNOWN-ACCESS ()) (601 'option_path' (CHARACTER
1 0 0 CHARACTER ((CONSTANT (INTEGER 4 0 0 INTEGER ()) 0 '8192'))) () (
UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0)
UNKNOWN-ACCESS (CONSTANT (CHARACTER 1 0 0 CHARACTER (())) 0 8192
'/uninitialised_path/                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            '))
(602 'mesh' (DERIVED 72 0 0 DERIVED ()) () (UNKNOWN-FL UNKNOWN-INTENT
UNKNOWN-PROC UNKNOWN UNKNOWN 0 0) UNKNOWN-ACCESS (STRUCTURE (DERIVED 72
0 0 DERIVED ()) 0 ((() ()) ((CONSTANT (LOGICAL 4 0 0 LOGICAL ()) 0 1) ())
((STRUCTURE (DERIVED 89 0 0 DERIVED ()) 0 ((() ()) (() ()) (() ()) (() ())
((NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0) ()) ((NULL (UNKNOWN 0 0 0 UNKNOWN ())
0) ()) ((NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0) ()) ((NULL (UNKNOWN 0 0 0
UNKNOWN ()) 0) ()) ((NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0) ()) ((NULL (
UNKNOWN 0 0 0 UNKNOWN ()) 0) ()) ((NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0) ())
((STRUCTURE (DERIVED 90 0 0 DERIVED ()) 0 ((() ()) (() ()) (() ()) (() ())
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
0) ())) ())) (603 'refcount' (DERIVED 97 0 0 DERIVED ()) () (UNKNOWN-FL
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 POINTER) UNKNOWN-ACCESS
(NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0)) (604 'aliased' (LOGICAL 4 0 0
LOGICAL ()) () (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0
0) UNKNOWN-ACCESS (CONSTANT (LOGICAL 4 0 0 LOGICAL ()) 0 0))) PUBLIC (()
() () ()) () 0 0 25110185)
473 'integer_vector' 'futils' 'integer_vector' 1 ((DERIVED
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 POINTER_COMP) (UNKNOWN 0
0 0 UNKNOWN ()) 0 0 () () 0 ((605 'ptr' (INTEGER 4 0 0 INTEGER ()) (1 0
DEFERRED () ()) (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN
0 0 DIMENSION POINTER) UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0)))
PUBLIC (() () () ()) () 0 0 9661976)
483 'integer_set' 'integer_set_module' 'integer_set' 1 ((DERIVED
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0) (UNKNOWN 0 0 0 UNKNOWN
()) 0 0 () () 0 ((606 'address' (DERIVED 607 0 0 DERIVED ()) () (
UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0)
UNKNOWN-ACCESS ())) PUBLIC (() () () ()) () 0 0 67217964)
607 'c_ptr' '__iso_c_binding' '' 1 ((DERIVED UNKNOWN-INTENT UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0 IS_BIND_C IS_C_INTEROP IS_ISO_C) (DERIVED 607 1 1
UNKNOWN ()) 0 0 () () 0 ((608 '__c_ptr_c_address' (INTEGER 8 1 0 INTEGER
()) () (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0)
UNKNOWN-ACCESS ())) UNKNOWN-ACCESS () () 2 39 0)
)

('derive_collapsed_bcs' 0 66 'get_dg_surface_mesh' 0 71
'get_periodic_boundary_condition' 0 75 'has_scalar_surface_field' 0 79
'set_dirichlet_consistent' 0 48)
