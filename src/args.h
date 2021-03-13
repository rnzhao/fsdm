/*
   Copyright (c) 2020 Roy Zhao <roy.zhao@uci.edu>
   Mozilla Public License Version 2.0
*/

#ifndef FSDM_ARGS_H
#define FSDM_ARGS_H

#include <stdbool.h>
#include <stddef.h>

typedef struct args {
    const char *fasta_file;
    const char **fastq_files;
    char *outfile;
    bool output_all;
    int num_fastq_pairs;
    int bc_mismatches;
    int ad_fl_mismatches;
    int ed_threshold;
} args;

extern args parse_args(int argc, const char **argv);

#endif
