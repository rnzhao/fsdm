/*
   Copyright (c) 2020 Roy Zhao <roy.zhao@uci.edu>
   Mozilla Public License Version 2.0
*/

#ifndef BC_HASH_H
#define BC_HASH_H

#include <stdbool.h>
#include <stddef.h>
#include <stdint.h>

typedef struct hash_kv hash_kv;
typedef struct bc_hash_table bc_hash_table;

bc_hash_table *init_hash_table(size_t num_items);

void hash_table_insert(bc_hash_table *self,
                       const char *key,
                       int8_t value,
                       size_t bc_index,
                       bool overwrite);

int hash_table_lookup(const bc_hash_table *self,
                      const char *key,
                      size_t bc_index);

void prune_hash_table(bc_hash_table **ht_double_ptr);

void destroy_hash_table(bc_hash_table **ht_double_ptr);

#endif
