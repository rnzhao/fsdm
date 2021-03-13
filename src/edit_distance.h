/*
   Copyright (c) 2020 Roy Zhao <roy.zhao@uci.edu>
   Mozilla Public License Version 2.0
*/

#ifndef EDIT_DISTANCE_H
#define EDIT_DISTANCE_H

#include <stddef.h>
#include <stdint.h>

extern uint64_t uint_pow(unsigned long base, unsigned long exponent);

extern unsigned long calc_num_combos(unsigned int bc_len,
                                     unsigned int num_barcodes,
                                     unsigned int max_mismatches);

extern void generate_all_bc_combos(size_t len,
                                   char (*output)[len + 1]);

extern unsigned int hamming_distance(const char *restrict seq_1,
                                     const char *restrict seq_2,
                                     size_t length);

extern int damerau_levenshtein(const char *restrict seq_1,
                               const char *restrict seq_2,
                               const int len_1);

extern int nw_offset(const char *restrict seq_1,
                     const char *restrict seq_2,
                     size_t len);

#endif
