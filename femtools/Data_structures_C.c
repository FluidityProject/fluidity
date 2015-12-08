#include "Judy.h"
#include "confdefs.h"
#include "stdio.h"
#include "assert.h"
#include "string.h"

/* To understand these, read
   http://judy.sourceforge.net/doc/Judy1_3x.htm and
   http://judy.sourceforge.net/doc/JudyL_3x.htm */

void integer_set_create_c(Pvoid_t* i)
{
  *i = (Pvoid_t) NULL;
}

void integer_set_delete_c(Pvoid_t* i)
{
  Word_t mem_freed;
  Pvoid_t ptr = (Pvoid_t) *i;
  J1FA(mem_freed, ptr);
  *i = ptr;
}

void integer_set_insert_c(Pvoid_t* i, int* v, int* c)
{
  int changed;
  Word_t index = *v;
  Pvoid_t ptr = (Pvoid_t) *i;
  J1S(changed, ptr, index);
  *c = changed;
  *i = ptr;
}

void integer_set_length_c(Pvoid_t* i, int* l)
{
  Word_t len;
  Pvoid_t ptr = (Pvoid_t) *i;
  J1C(len, ptr, 0, -1);
  *l = len;
  *i = ptr;
}

void integer_set_fetch_c(Pvoid_t* i, int* idx, int* val)
{
  Word_t index = *idx, value;
  int worked;
  Pvoid_t ptr = (Pvoid_t) *i;
  J1BC(worked, ptr, index, value); 
  assert(worked == 1);
  *val = value;
  *i = ptr;
}

void integer_set_remove_c(Pvoid_t* i, int* idx, int* status)
{
  Word_t index = *idx;
  int stat;
  Pvoid_t ptr = (Pvoid_t) *i;
  J1U(stat, ptr, index); 
  *status = stat;
  *i = ptr;
}

void integer_set_has_value_c(Pvoid_t* i, int* val, int* present)
{
  Word_t value = *val;
  int wpresent;
  Pvoid_t ptr = (Pvoid_t) *i;
  J1T(wpresent, ptr, value); 
  *present = wpresent;
}

void integer_hash_table_create_c(Pvoid_t* i)
{
  assert(sizeof(int*) == sizeof(Pvoid_t));
  *i = (Pvoid_t) NULL;
}

void integer_hash_table_delete_c(Pvoid_t* i)
{
  Word_t mem_freed;
  Pvoid_t ptr = (Pvoid_t) *i;
  JLFA(mem_freed, ptr);
  *i = ptr;
}

void integer_hash_table_insert_c(Pvoid_t* i, int* k, int* v)
{
  Word_t key = *k;
  PWord_t pvalue;
  Pvoid_t ptr = (Pvoid_t) *i;
  JLI(pvalue, ptr, key);
  *pvalue = *v;
  *i = ptr;
}

void integer_hash_table_length_c(Pvoid_t* i, int* l)
{
  Word_t len;
  Pvoid_t ptr = (Pvoid_t) *i;
  JLC(len, ptr, 0, -1);
  *l = len;
  *i = ptr;
}

void integer_hash_table_fetch_c(Pvoid_t* i, int* k, int* v)
{
  Word_t key = *k, value;
  PWord_t pvalue = &value;
  Pvoid_t ptr = (Pvoid_t) *i;
  JLG(pvalue, ptr, key); 
  if (pvalue == NULL)
  {
    fprintf(stderr, "Error: hash table has no key %d\n", *k);
    assert(pvalue != NULL);
  }
  *v = *pvalue;
  *i = ptr;
}

void integer_hash_table_remove_c(Pvoid_t* i, int* k, int* status)
{
  Word_t key = *k;
  int stat;
  Pvoid_t ptr = (Pvoid_t) *i;
  JLD(stat, ptr, key); 
  *status = stat;
  *i = ptr;
}

void integer_hash_table_has_key_c(Pvoid_t* i, int* k, int* present)
{
  Word_t key = *k, value;
  PWord_t pvalue = &value;
  Pvoid_t ptr = (Pvoid_t) *i;
  JLG(pvalue, ptr, key); 
  *present = (pvalue != NULL);
  *i = ptr;
}

void integer_hash_table_fetch_pair_c(Pvoid_t* i, int* idx, int* key, int* val)
{
  Word_t nth = *idx; /* what Judy calls nth is what I am calling index */
  Word_t index; /* what Judy calls index is what I am calling key */
  PWord_t pvalue;
  Pvoid_t ptr = (Pvoid_t) *i;
  JLBC(pvalue, ptr, nth, index); 
  assert(pvalue != NULL);
  *key = index;
  *val = *pvalue;
  *i = ptr;
}

void string_hash_table_create_c(Pvoid_t* i)
{
  assert(sizeof(int*) == sizeof(Pvoid_t));
  *i = (Pvoid_t) NULL;
}

void string_hash_table_delete_c(Pvoid_t* i)
{
  Word_t mem_freed;
  Pvoid_t ptr = (Pvoid_t) *i;
  JSLFA(mem_freed, ptr);
  *i = ptr;
}

void string_hash_table_insert_c(Pvoid_t* i, uint8_t* key, int* v)
{
  PWord_t pvalue;
  Pvoid_t ptr = (Pvoid_t) *i;
  JSLI(pvalue, ptr, key);
  *pvalue = *v;
  *i = ptr;
}

void string_hash_table_insert_pointer_c(Pvoid_t* i, uint8_t* key, void* p)
{
  PWord_t pvalue;
  Pvoid_t ptr = (Pvoid_t) *i;
  JSLI(pvalue, ptr, key);
  *pvalue = *(PWord_t)p;
  *i = ptr;
}

void string_hash_table_fetch_c(Pvoid_t* i, uint8_t* key, int* v)
{
  PWord_t pvalue;
  Pvoid_t ptr = (Pvoid_t) *i;
  JSLG(pvalue, ptr, key); 
  if (pvalue == NULL)
  {
    fprintf(stderr, "Error: hash table has no key %s\n", (char*) key);
    assert(pvalue != NULL);
  }
  *v = *pvalue;
  *i = ptr;
}

void string_hash_table_fetch_pointer_c(Pvoid_t* i, uint8_t* key, void* p)
{
  PWord_t pvalue;
  Pvoid_t ptr = (Pvoid_t) *i;
  JSLG(pvalue, ptr, key); 
  if (pvalue == NULL)
  {
    fprintf(stderr, "Error: hash table has no key %s\n", (char*) key);
    assert(pvalue != NULL);
  }
  *(PWord_t) p = *pvalue;
  *i = ptr;
}

void string_hash_table_remove_c(Pvoid_t* i, uint8_t* key, int* status)
{
  int stat;
  Pvoid_t ptr = (Pvoid_t) *i;
  JSLD(stat, ptr, key); 
  *status = stat;
  *i = ptr;
}

void string_hash_table_has_key_c(Pvoid_t* i, uint8_t* key, int* present)
{
  PWord_t pvalue ;
  Pvoid_t ptr = (Pvoid_t) *i;

  JSLG(pvalue, ptr, key); 
  *present = (pvalue != NULL);
  *i = ptr;
}

void string_hash_table_get_first_c(Pvoid_t* i, uint8_t* key, int* key_len, int* v, int* status)
{
  key[0] = '\0';

  PWord_t pvalue ;
  Pvoid_t ptr = (Pvoid_t) *i;
  JSLF(pvalue, ptr, key); 
  *status = (pvalue != NULL);
  *v = *pvalue;
  *i = ptr;
  if (*status){ *key_len = strlen((char*) key);} else { *key_len=0;};
}

void string_hash_table_get_first_pointer_c(Pvoid_t* i, uint8_t* key, int* key_len, void* p, int* status)
{
  key[0] = '\0';

  PWord_t pvalue ;
  Pvoid_t ptr = (Pvoid_t) *i;
  JSLF(pvalue, ptr, key); 
  *status = (pvalue != NULL);
  *i = ptr;
  if (*status){ 
    *key_len = strlen((char*) key);
    *(PWord_t) p = *pvalue;
  } else { 
    *key_len=0;
    *(PWord_t) p = NULL;
  };
}

void string_hash_table_get_next_c(Pvoid_t* i, uint8_t* key, int* key_len, int* v, int* status)
{
  PWord_t pvalue ;
  Pvoid_t ptr = (Pvoid_t) *i;
  JSLN(pvalue, ptr, key); 
  *status = (pvalue != NULL);
  *v = *pvalue;
  *i = ptr;
  if (*status){ *key_len = strlen((char*) key);} else { *key_len=0;};
}

void string_hash_table_get_next_pointer_c(Pvoid_t* i, uint8_t* key, int* key_len, void* p, int* status)
{
  PWord_t pvalue ;
  Pvoid_t ptr = (Pvoid_t) *i;
  JSLN(pvalue, ptr, key); 
  *status = (pvalue != NULL);
  *i = ptr;
  if (*status){ 
    *key_len = strlen((char*) key);
    *(PWord_t) p = *pvalue;
  } else { 
    *key_len=0;
    *(PWord_t) p = NULL;
  };
}

