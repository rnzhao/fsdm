/*
   Copyright (c) 2020 Roy Zhao <roy.zhao@uci.edu>
   Mozilla Public License Version 2.0
*/

#include "bc_hash.h"

#include <assert.h>
#include <errno.h>
#include <stdbool.h>
#include <stddef.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#if defined __clang__ || defined __GNUC__ || defined __INTEL_COMPILER
    #define LIKELY(x) __builtin_expect(!!(x), 1)
#else
    #define LIKELY(x) (x)
#endif


/* The value is one of:
    0 -> empty slot
   -1 -> ambiguous/duplicated entry
    1 to num_barcodes -> unique entry */
struct hash_kv {
    char key[6];
    int8_t value[2];
};

struct bc_hash_table {
    void *malloc_ptr;
    size_t num_slots;
    size_t key_length;
    size_t num_items;
    hash_kv items[];
};


/* 32-bit FNV-1a hash function */
inline uint32_t fnv_1a(const void *buffer,
                       size_t length)
{
    const uint32_t FNV_PRIME = 16777619UL;
    const uint32_t FNV_OFFSET_32 = 2166136261UL;

    const uint8_t *restrict octets = buffer;

    uint32_t hash = FNV_OFFSET_32;

    LIKELY(length == 6);
    #pragma unroll 6
    for (size_t i = 0; i < length; i++) {
        hash ^= octets[i];
        hash *= FNV_PRIME;
    }

    return hash;
}


inline size_t nearest_pow2(size_t number)
{
    size_t value = 1;

    while (value < number) {
        value *= 2;
    }

    return value;
}


bc_hash_table *init_hash_table(size_t num_items)
{
    // Maximum 0.67 load factor
    size_t num_slots_non_pow2 = num_items * 3 / 2;
    size_t num_slots = nearest_pow2(num_slots_non_pow2);

    // Allocate memory for flat hash table with an extra margin for alignment
    bc_hash_table *hash_table;
    void *malloc_ptr = malloc(sizeof(*hash_table) + num_slots * sizeof(*(hash_table->items)) + 127);

    if (malloc_ptr == NULL) {
        perror("Error: memory allocation failed for hash table");
        exit(EXIT_FAILURE);
    }

    // Align key-value struct array member to typical cache line
    uintptr_t aligned_address = ((uintptr_t) malloc_ptr + 63) & ~63;
    aligned_address += 64;

    hash_table = (bc_hash_table*) (aligned_address - offsetof(bc_hash_table, items));

    assert(((uintptr_t) hash_table->items & 63) == 0);

    hash_table->malloc_ptr = malloc_ptr;
    hash_table->num_slots = num_slots;
    hash_table->key_length = sizeof((*(hash_table->items)).key);
    hash_table->num_items = 0;

    memset(hash_table->items, 0, num_slots * sizeof(*(hash_table->items)));

    return hash_table;
}


/* Inserts an entry into a hash table if
   the key is not already present, otherwise,
   flag the existing entry as a duplicate. */
void hash_table_insert(bc_hash_table *ht,
                       const char *key,
                       int8_t value,
                       size_t bc_index,
                       bool overwrite)
{
    uint32_t hash = fnv_1a(key, ht->key_length);
    size_t index = hash & (ht->num_slots - 1);

    while (true) {
        if (index >= ht->num_slots) {
            index = 0;
        }

        hash_kv *kv_slot = &(ht->items[index]);

        // Slot is empty
        if (kv_slot->value[bc_index] == 0) {
            // Copy new key/value pair into hash table
            strncpy(kv_slot->key, key, ht->key_length);

            kv_slot->value[bc_index] = value;
            ht->num_items += 1;

            return;
        }
        // Key to be inserted matches that of an existing entry,
        // i.e., current key is a repeat and not a collision
        else if (strncmp(key, kv_slot->key, ht->key_length) == 0) {
            if (overwrite) {
                kv_slot->value[bc_index] = value;
                return;
            }
            else {
                kv_slot->value[bc_index] = -1;
                return;
            }
        }
        // Collision; probe next slot
        else {
            index++;
        }
    }
}


inline int hash_table_lookup(const bc_hash_table *ht,
                             const char *key,
                             size_t bc_index)
{
    uint32_t hash = fnv_1a(key, ht->key_length);
    size_t index = hash & (ht->num_slots - 1);

    int8_t value;

    do {
        if (index >= ht->num_slots) {
            index = 0;
        }

        value = ht->items[index].value[bc_index];

        if (value != 0 && strncmp(key, ht->items[index].key, ht->key_length) == 0) {
            break;
        }

        index++;
    } while (value != 0);

    return value;
}


/* Remove entries from the hash table that
   refer to ambiguous barcode mismatches. */
void prune_hash_table(bc_hash_table **ht_double_ptr)
{
    if (*ht_double_ptr == NULL) {
        return;
    }

    bc_hash_table *original_table = *ht_double_ptr;
    size_t num_unique_keys = 0;
    size_t num_slots_with_duplicates = 0;

    for (size_t i = 0; i < original_table->num_slots; i++) {
        if (original_table->items[i].value[0] > 0 ||
            original_table->items[i].value[1] > 0) {
            num_unique_keys++;
        }
        else if (original_table->items[i].value[0] < 0 ||
                 original_table->items[i].value[1] < 0) {
            num_slots_with_duplicates++;
        }
    }

    if (num_slots_with_duplicates == 0) {
        return;
    }

    size_t pruned_num_slots = nearest_pow2(num_unique_keys * 3 / 2);
    bc_hash_table *pruned_table = init_hash_table(pruned_num_slots);

    if (pruned_table == NULL) {
        perror("Error: memory allocation failed for pruned hash table");
        exit(EXIT_FAILURE);
    }

    for (size_t i = 0; i < original_table->num_slots; i++) {
        for (size_t j = 0; j < 2; j++) {
            if (original_table->items[i].value[j] > 0) {
                hash_table_insert(pruned_table, original_table->items[i].key,
                                  original_table->items[i].value[j], j, false);
            }
        }
    }

    destroy_hash_table(ht_double_ptr);

    *ht_double_ptr = pruned_table;
}


void destroy_hash_table(bc_hash_table **ht_double_ptr)
{
    bc_hash_table *hash_table = *ht_double_ptr;

    if (hash_table->malloc_ptr != NULL) {
        free(hash_table->malloc_ptr);
    }

    *ht_double_ptr = NULL;
}
