/*
   Copyright (c) 2020 Roy Zhao <roy.zhao@uci.edu>
   Mozilla Public License Version 2.0
*/

#include "edit_distance.h"

#include <stdarg.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>


static inline int min(const size_t num_items, ...)
{
    va_list items;
    va_start(items, num_items);

    int min_value = va_arg(items, int);

    for (size_t i = 1; i < num_items; i++) {
        int current_item = va_arg(items, int);

        if (current_item < min_value) {
            min_value = current_item;
        }
    }

    va_end(items);

    return min_value;
}


static inline uint64_t factorial(unsigned int n)
{
    static const uint64_t n_factorial[21] = {
        1, 1, 2, 6, 24,
        120, 720, 5040, 40320,
        362880, 3628800, 39916800,
        479001600, 6227020800, 87178291200,
        1307674368000, 20922789888000, 355687428096000,
        6402373705728000, 121645100408832000, 2432902008176640000
    };

    if (n >= sizeof(n_factorial) / sizeof(n_factorial[0])) {
        return 0;
    }

    return n_factorial[n];
}


static inline uint64_t n_choose_k(unsigned int n,
                                  unsigned int k)
{
    return factorial(n) / (factorial(k) * factorial(n - k));
}


uint64_t uint_pow(unsigned long base,
                  unsigned long exponent)
{
    uint64_t value = 1;

    while (exponent > 0) {
        if (exponent & 1) {
            value *= base;
        }

        base *= base;
        exponent /= 2;
    }

    return value;
}


unsigned long calc_num_combos(unsigned int bc_len,
                              unsigned int num_barcodes,
                              unsigned int max_mismatches)
{
    if (max_mismatches == 0 || max_mismatches >= bc_len) {
        return num_barcodes;
    }

    unsigned long total_possible_combos = 0;

    for (size_t k = 1; k <= max_mismatches; k++) {
        unsigned long k_mismatch_combos = uint_pow(4, k);
        k_mismatch_combos *= n_choose_k(bc_len, k);
        total_possible_combos += k_mismatch_combos;
    }

    total_possible_combos *= num_barcodes;

    unsigned long max_possible_permutations = uint_pow(4, bc_len);

    if (total_possible_combos > max_possible_permutations) {
        total_possible_combos = max_possible_permutations;
    }

    return total_possible_combos;
}


void generate_all_bc_combos(size_t len,
                            char (*output)[len + 1])
{
    const char bases[4] = {'A', 'T', 'G', 'C'};

    for (size_t i = 0; i < len; i++) {
        unsigned long prev_num_combos = uint_pow(4, i);

        for (size_t j = 1; j < 4; j++) {
            for (size_t k = 0; k < prev_num_combos; k++) {
                strncpy(output[j * prev_num_combos + k] + len - i, output[k] + len - i, i);
            }
        }

        for (size_t j = 0; j < 4; j++) {
            for (size_t k = 0; k < prev_num_combos; k++) {
                output[j * prev_num_combos + k][len - i - 1] = bases[j];
            }
        }
    }
}


unsigned int hamming_distance(const char *restrict seq_1,
                              const char *restrict seq_2,
                              size_t length)
{
    unsigned int num_mismatches = 0;

    #pragma unroll 6
    #pragma clang loop vectorize(enable)
    #pragma GCC ivdep
    for (size_t i = 0; i < length; i++) {
        num_mismatches += (seq_1[i] != seq_2[i]);
    }

    return num_mismatches;
}


int damerau_levenshtein(const char *restrict seq_1,
                        const char *restrict seq_2,
                        const int len)
{
    int da[256] = {0};
    int max_dist = len + len;

    // Populate initial scores
    int dpm[len + 2][len + 2];
    dpm[0][0] = max_dist;

    for (int i = 0; i < len + 1; i++) {
        dpm[i + 1][0] = max_dist;
        dpm[i + 1][1] = i;
    }
    for (int j = 0; j < len + 1; j++) {
        dpm[0][j + 1] = max_dist;
        dpm[1][j + 1] = j;
    }

    // Fill in distance matrix
    for (int i = 1; i < len + 1; i++) {
        int db = 0;

        uint8_t da_i;
        int k, l, cost;

        for (int j = 1; j < len + 1; j++) {
            da_i = (uint8_t) seq_2[j - 1];
            k = da[da_i];
            l = db;

            if (seq_1[i - 1] == seq_2[j - 1]) {
                cost = 0;
                db = j;
            }
            else {
                cost = 1;
            }

            dpm[i + 1][j + 1] = min(4,
                                    dpm[i][j] + cost,                    // substitution
                                    dpm[i + 1][j] + 1,                   // insertion
                                    dpm[i][j + 1] + 1,                   // deletion
                                    dpm[k][l] + (i-k-1) + 1 + (j-l-1));  // transposition
        }

        da_i = (uint8_t) seq_1[i - 1];
        da[da_i] = i;
    }

    return dpm[len + 1][len + 1];
}


int nw_offset(const char *restrict seq_1,
              const char *restrict seq_2,
              size_t len)
{
    enum {
        MATCH = 1,
        MISMATCH = -1,
        INDEL = -1
    };

    enum {
        DIAG = 0,
        UP = 1,
        LEFT = 2
    };

    int scores[len + 1][len + 1];
    int traceback[len + 1][len + 1];

    for (size_t i = 0; i < len + 1; i++) {
        for (size_t j = 0; j < len + 1; j++) {
            scores[i][j] = 0;
            traceback[i][j] = 0;
        }
    }

    for (size_t i = 1; i < len + 1; i++) {
        scores[i][0] = INDEL * i;
        traceback[i][0] = UP;
    }
    for (size_t j = 1; j < len + 1; j++) {
        scores[0][j] = INDEL * j;
        traceback[0][j] = LEFT;
    }

    for (size_t i = 1; i < len + 1; i++) {
        for (size_t j = 1; j < len + 1; j++) {
            int s_m = (seq_1[i - 1] == seq_2[j - 1]) ? MATCH : MISMATCH;
            int overlap_score = scores[i - 1][j - 1] + s_m;

            int s1_gap = scores[i - 1][j] + INDEL;
            int s2_gap = scores[i][j - 1] + INDEL;
            int gap_score, gap_move;

            if (s1_gap > s2_gap) {
                gap_score = s1_gap;
                gap_move = UP;
            }
            else {
                gap_score = s2_gap;
                gap_move = LEFT;
            }

            if (overlap_score > gap_score) {
                scores[i][j] = overlap_score;
                traceback[i][j] = DIAG;
            }
            else {
                scores[i][j] = gap_score;
                traceback[i][j] = gap_move;
            }
        }
    }

    size_t aln_index = len * 2 - 1;
    size_t i_1 = len, i_2 = len;

    int offset = 0;
    int offset_index[2] = {aln_index};
    bool offset_reached = false;

    while ((! offset_reached) && (i_1 > 0 || i_2 > 0)) {
        int traceback_direction = traceback[i_1][i_2];

        if (traceback_direction == DIAG) {
            offset_reached = true;
            offset = offset_index[1] - offset_index[0];
        }
        else if (traceback_direction == LEFT) {
            offset_index[1] = aln_index - 1;
            i_2 -= 1;
        }
        else if (traceback_direction == UP) {
            offset_index[0] = aln_index - 1;
            i_1 -= 1;
        }

        aln_index--;
    }

    return offset;
}
