/*
   Copyright (c) 2020 Roy Zhao <roy.zhao@uci.edu>
   Mozilla Public License Version 2.0
*/

#ifndef PARSE_SEQ_H
#define PARSE_SEQ_H

#include <stdbool.h>
#include <stddef.h>

enum { MAX_SEQ_LEN = 304 };

enum {
    ALLELE_A,
    ALLELE_C,
    ALLELE_G,
    ALLELE_T
};

typedef struct barcode_t {
    char seq[7];
    int label;
} barcode_t;

typedef struct read_segment {
    size_t offset;
    size_t length;
    char seq[MAX_SEQ_LEN];
} read_segment;

typedef struct prototype {
    size_t allele_offset;
    read_segment *segments[5];
} prototype;

typedef struct library_seqs {
    barcode_t *barcodes[2];
    unsigned int num_barcodes[2];
    read_segment adapters[2];
    read_segment flanking[2];
    read_segment alleles[4];
    read_segment prototype_strings[2];
    prototype prototypes[2];
} library_seqs;

extern size_t allele_char_to_enum(char allele);
extern bool all_standard_barcodes(const library_seqs *fs2_seqs);
extern library_seqs *load_fasta_sequences(const char *filepath);
extern void parse_prototypes(library_seqs *fs2_seqs);

#endif
