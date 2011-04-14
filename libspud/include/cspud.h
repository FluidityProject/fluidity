/*  Copyright (C) 2006 Imperial College London and others.

    Please see the AUTHORS file in the main source directory for a full list
    of copyright holders.

    Prof. C Pain
    Applied Modelling and Computation Group
    Department of Earth Science and Engineering
    Imperial College London

    C.Pain@Imperial.ac.uk

    This library is free software; you can redistribute it and/or
    modify it under the terms of the GNU Lesser General Public
    License as published by the Free Software Foundation,
    version 2.1 of the License.

    This library is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    Lesser General Public License for more details.

    You should have received a copy of the GNU Lesser General Public
    License along with this library; if not, write to the Free Software
    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307
    USA
*/

#ifndef CSPUD_H
#define CSPUD_H

#include "spud_config.h"
#include "spud_enums.h"

#define cspud_clear_options F77_FUNC(cspud_clear_options, CSPUD_CLEAR_OPTIONS)
#define cspud_load_options F77_FUNC(cspud_load_options, CSPUD_LOAD_OPTIONS)
#define cspud_write_options F77_FUNC(cspud_write_options, CSPUD_WRITE_OPTIONS)
#define cspud_get_child_name F77_FUNC(cspud_get_child_name, CSPUD_GET_CHILD_NAME)
#define cspud_number_of_children F77_FUNC(cspud_number_of_children, CSPUD_GET_NUMBER_OF_CHILDREN)
#define cspud_option_count F77_FUNC(cspud_option_count, CSPUD_OPTION_COUNT)
#define cspud_have_option F77_FUNC(cspud_have_option, CSPUD_HAVE_OPTION)
#define cspud_get_option_type F77_FUNC(cspud_get_option_type, CSPUD_GET_OPTION_TYPE)
#define cspud_get_option_rank F77_FUNC(cspud_get_option_rank, CSPUD_GET_OPTION_RANK)
#define cspud_get_option_shape F77_FUNC(cspud_get_option_shape, CSPUD_GET_OPTION_SHAPE)
#define cspud_get_option F77_FUNC(cspud_get_option, CSPUD_GET_OPTION)
#define cspud_add_option F77_FUNC(cspud_add_option, CSPUD_ADD_OPTION)
#define cspud_set_option F77_FUNC(cspud_set_option, CSPUD_SET_OPTION)
#define cspud_set_option_attribute F77_FUNC(cspud_set_option_attribute, CSPUD_SET_OPTION_ATTRIBUTE)
#define cspud_move_option F77_FUNC(cspud_move_option, CSPUD_MOVE_OPTION)
#define cspud_copy_option F77_FUNC(cspud_copy_option, CSPUD_COPY_OPTION)
#define cspud_delete_option F77_FUNC(cspud_delete_option, CSPUD_DELETE_OPTION)
#define cspud_print_options F77_FUNC(cspud_print_options, CSPUD_PRINT_OPTIONS)

#ifdef __cplusplus
extern "C" {
#endif

  void cspud_clear_options();

  void cspud_load_options(const char* filename, const int* filename_len);
  void cspud_write_options(const char* filename, const int* filename_len);

  int cspud_get_child_name(const char* key, const int* key_len, const int* index, char* child_name, const int* child_name_len);

  int cspud_number_of_children(const char* key, const int* key_len);

  int cspud_option_count(const char* key, const int* key_len);

  int cspud_have_option(const char* key, const int* key_len);

  int cspud_get_option_type(const char* key, const int* key_len, int* type);
  int cspud_get_option_rank(const char* key, const int* key_len, int* rank);
  int cspud_get_option_shape(const char* key, const int* key_len, int* shape);

  int cspud_get_option(const char* key, const int* key_len, void* val);

  int cspud_add_option(const char* key, const int* key_len);

  int cspud_set_option(const char* key, const int* key_len, const void* val, const int* type, const int* rank, const int* shape);

  int cspud_set_option_attribute(const char* key, const int* key_len, const char* val, const int* val_len);

  int cspud_move_option(const char* key1, const int* key1_len, const char* key2, const int* key2_len);
  int cspud_copy_option(const char* key1, const int* key1_len, const char* key2, const int* key2_len);
   
  int cspud_delete_option(const char* key, const int* key_len);

  void cspud_print_options();

#ifdef __cplusplus
}
#endif

#endif
