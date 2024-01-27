/* ---------------------------------------------------------------------------


  This file is part of the ``utils'' package of NuSMV version 2.
  Copyright (C) 1998-2001 by CMU and FBK-irst.

  NuSMV version 2 is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public
  License as published by the Free Software Foundation; either
  version 2 of the License, or (at your option) any later version.

  NuSMV version 2 is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  Lesser General Public License for more details.

  You should have received a copy of the GNU Lesser General Public
  License along with this library; if not, write to the Free Software
  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307  USA.

  For more information on NuSMV see <http://nusmv.fbk.eu>
  or email to <nusmv-users@fbk.eu>.
  Please report bugs to <nusmv-users@fbk.eu>.

  To contact the NuSMV development board, email to <nusmv@fbk.eu>. 

-----------------------------------------------------------------------------*/

/*!
  \author Originated from glu library of VIS
  \brief Generic array manipulator

  \todo: Missing description

*/


#include "nusmv/core/utils/StreamMgr.h"
#include <stdio.h>
#include "cudd/util.h"
#include "nusmv/core/utils/array.h"
#include "nusmv/core/utils/error.h"

/*!
  \brief \todo Missing synopsis

  \todo Missing description
*/
#define INIT_SIZE       3

unsigned int array_global_index;
int array_global_insert;


/* Allocate an array of 'number' elements,each of which require 'size' bytes */
array_t* array_do_alloc(int size, int number)
{
  array_t *array;

  array = ALLOC(array_t, 1);
  if (array == NIL(array_t)) {
    return NIL(array_t);
  }
  array->num = 0;
  array->n_size = MAX(number, INIT_SIZE);
  array->obj_size = size;
  array->index = -size;
  array->e_index = 0;
  array->e_insert = 0;
  array->space = ALLOC(char, array->n_size * array->obj_size);

  if (array->space == NIL(char)) {
    FREE(array);
    return NIL(array_t);
  }
  (void) memset(array->space, 0, array->n_size * array->obj_size);
  return array;
}


void array_free(array_t *array)
{
  if (array == NIL(array_t)) return;
  if (array->index >= 0) array_abort(array,4);
  FREE(array->space);
  FREE(array);
}


array_t* array_dup(array_t *old)
{
  array_t *new;

  new = ALLOC(array_t, 1);
  if (new == NIL(array_t)) {
    return NIL(array_t);
  }
  new->num = old->num;
  new->n_size = old->num;
  new->obj_size = old->obj_size;
  new->index = -new->obj_size;
  new->e_index = 0;
  new->e_insert = 0;
  new->space = ALLOC(char, new->n_size * new->obj_size);

  if (new->space == NIL(char)) {
    FREE(new);
    return NIL(array_t);
  }
  (void) memcpy(new->space, old->space, old->num * old->obj_size);
  return new;
}


/* append the elements of array2 to the end of array1 */
int array_append(array_t* array1, array_t* array2)
{
  char *pos;

  if (array1->index >= 0) array_abort(array1,4);
  if (array1->obj_size != array2->obj_size) {
    array_abort(array1,2);
    /* NOTREACHED */
  }

  /* make sure array1 has enough room */
  if (array1->n_size < array1->num + array2->num) {
    if (array_resize(array1, array1->num + array2->num) == ARRAY_OUT_OF_MEM) {
      return ARRAY_OUT_OF_MEM;
    }
  }
  pos = array1->space + array1->num * array1->obj_size;
  (void) memcpy(pos, array2->space, array2->num * array2->obj_size);
  array1->num += array2->num;

  return 1;
}


/* join array1 and array2, returning a new array */
array_t* array_join(array_t* array1, array_t* array2)
{
  array_t *array;
  char *pos;

  if (array1->obj_size != array2->obj_size) {
    array_abort(array1,3);
    fail("array: join not defined for arrays of different sizes\n");
    /* NOTREACHED */
  }
  array = ALLOC(array_t, 1);
  if (array == NIL(array_t)) {
    return NIL(array_t);
  }
  array->num = array1->num + array2->num;
  array->n_size = array->num;
  array->obj_size = array1->obj_size;
  array->index = -array->obj_size;
  array->e_index = 0;
  array->e_insert = 0;
  array->space = ALLOC(char, array->n_size * array->obj_size);
  if (array->space == NIL(char)) {
    FREE(array);
    return NIL(array_t);
  }
  (void) memcpy(array->space, array1->space, array1->num * array1->obj_size);
  pos = array->space + array1->num * array1->obj_size;
  (void) memcpy(pos, array2->space, array2->num * array2->obj_size);
  return array;
}


char* array_do_data(array_t *array)
{
  char *data;

  data = ALLOC(char, array->num * array->obj_size);
  if (data == NIL(char)) {
    return NIL(char);
  }
  (void) memcpy(data, array->space, array->num * array->obj_size);
  return data;
}


int array_resize(array_t* array, int new_size)
{
  int old_size;
  char *pos, *newspace;

  /* Note that this is not an exported function, and does not check if
     the array is locked since that is already done by the caller. */
  old_size = array->n_size;
  array->n_size = MAX(array->n_size * 2, new_size);
  newspace = REALLOC(char, array->space, array->n_size * array->obj_size);
  if (newspace == NIL(char)) {
    array->n_size = old_size;
    return ARRAY_OUT_OF_MEM;
  } else {
    array->space = newspace;
  }
  pos = array->space + old_size * array->obj_size;
  (void) memset(pos, 0, (array->n_size - old_size)*array->obj_size);
  return 1;
}

void array_sort(array_t* array, int (*compare)(const void*, const void*))
{
  qsort((void *)array->space, array->num, array->obj_size, compare);
}

void array_uniq(array_t* array,
                int (*compare)(char**, char**),
                void (*free_func)(char*))
{
  int i, last;
  char *dest, *obj1, *obj2;

  dest = array->space;
  obj1 = array->space;
  obj2 = array->space + array->obj_size;
  last = array->num;

  for(i = 1; i < last; i++) {
    if ((*compare)((char **) obj1, (char **) obj2) != 0) {
      if (dest != obj1) {
        (void) memcpy(dest, obj1, array->obj_size);
      }
      dest += array->obj_size;
    } else {
      if (free_func != 0) (*free_func)(obj1);
      array->num--;
    }
    obj1 += array->obj_size;
    obj2 += array->obj_size;
  }
  if (dest != obj1) {
    (void) memcpy(dest, obj1, array->obj_size);
  }
}

int array_abort(const array_t* a, int i)
{
  switch (i) {

  case 0:           /* index error on insert */
    error_unreachable_code_msg("array_t: insert of %d\n", a->index);
    break;

  case 1:           /* index error on fetch */
    error_unreachable_code_msg("array_t: fetch index %d not in [0,%d]\n",
                               a->e_index, a->num-1);
    break;

  case 2:           /* append with different element sizes */
    error_unreachable_code_msg("array_t: append undefined for arrays of different sizes\n");
    break;

  case 3:           /* join with different element sizes */
    error_unreachable_code_msg("array_t: join not defined for arrays of different sizes\n");
    break;

  case 4:           /* size error or locked error */
    if (a->index >= 0) {
      /* Since array_insert is a macro, it is not allowed to nest a
         call to any routine which might move the array space through
         a realloc or free inside an array_insert call. */
      error_unreachable_code_msg("array_t: nested insert, append, or free operations\n");
    }
    else {
      error_unreachable_code_msg("array_t: object size mismatch\n");
    }
    break;
  default:
    error_unreachable_code_msg("array_t: unknown error");
    break;
  }
}
