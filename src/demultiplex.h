/*
   Copyright (c) 2020 Roy Zhao <roy.zhao@uci.edu>
   Mozilla Public License Version 2.0
*/

#ifndef DEMULTIPLEX_H
#define DEMULTIPLEX_H

#include "bc_hash.h"
#include "parse_seq.h"

#include <stddef.h>

typedef struct bc_counter {
    unsigned int num_bc1;
    unsigned int num_bc2;
    unsigned int counts[][4];
} bc_counter;

void demultiplex_fastq_pair(const char **fastq_pair,
                            library_seqs *fs2_seqs,
                            bc_hash_table *hash_table,
                            bc_counter *bc_combo_counts,
                            int ad_fl_mismatches,
                            int ed_threshold,
                            const bool *valid_alleles);

#endif
