/*
   Copyright (c) 2020 Roy Zhao <roy.zhao@uci.edu>
   Mozilla Public License Version 2.0
*/

#include "fs2_barcodes.h"
#include "parse_seq.h"
#include "kseq.h"

// #include <ctype.h>
#include <errno.h>
#include <stdbool.h>
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <zlib.h>

KSEQ_INIT(gzFile, gzread)


static inline void str_toupper(char *str)
{
    for (size_t i = 0; str[i] != '\0'; i++) {
        str[i] = (char) toupper((unsigned char) str[i]);
    }
}


static bool str_in_list(const char *str,
                        const char **list,
                        size_t length)
{
    size_t i, j;

    for (i = 0, j = 0; i < length; i++) {
        j++;

        if (strcmp(str, list[i]) == 0) {
            break;
        }
    }

    return i < j;
}


static int num_fasta_seqs_read(const read_segment *arr,
                               size_t arr_length)
{
    int n = 0;

    for (size_t i = 0; i < arr_length; i++) {
        n += (arr[i].seq[0] != '\0');
    }

    return n;
}


static bool valid_fasta_seq(const kseq_t *seq)
{
    const char *valid_names[] = {"bc1", "bc2", "allele",
                                 "adapter1", "adapter2",
                                 "flanking1", "flanking2",
                                 "prototype1", "prototype2"};
    const size_t num_valid_names = sizeof(valid_names) / sizeof(valid_names[0]);

    if (! str_in_list(seq->name.s, valid_names, num_valid_names)) {
        fprintf(stderr, "Error: unrecognized sequence type '%s'\n", seq->name.s);
        return false;
    }

    bool seq_is_valid = true;

    if (seq->seq.l > MAX_SEQ_LEN) {
        fprintf(stderr, "Error: sequence '%s' is too long\n", seq->name.s);
        seq_is_valid = false;
    }
    else if (seq->comment.l > MAX_SEQ_LEN) {
        fprintf(stderr, "Error: comment for sequence '%s' is too long\n", seq->name.s);
        seq_is_valid = false;
    }
    else if (strncmp(seq->name.s, "bc", 2) == 0) {
        if (seq->seq.l != 6) {
            fprintf(stderr, "Error: invalid barcode length: '%s'\n", seq->seq.s);
            seq_is_valid = false;
        }

        char *end = NULL;
        strtol(seq->comment.s, &end, 10);

        if (*end != '\0') {
            fprintf(stderr, "Error: barcode label must be an integer: '%s'\n", seq->comment.s);
            seq_is_valid = false;
        }
    }
    else if (strcmp(seq->name.s, "allele") == 0) {
        char allele[2] = {(char) toupper((unsigned char) seq->seq.s[0]), '\0'};

        if (seq->seq.l != 1 || strpbrk(allele, "ACGT") == NULL) {
            fprintf(stderr, "Error: invalid allele '%s'\n", seq->seq.s);
            seq_is_valid = false;
        }
    }
    else if (strncmp(seq->name.s, "prototype", 9) == 0) {
        char prototype_str[strlen(seq->seq.s)];
        strcpy(prototype_str, seq->seq.s);
        char *segment = strtok(prototype_str, "|");
        int num_segments = 0;

        if (strncmp(seq->seq.s, "bc", 2) != 0) {
            fprintf(stderr, "Error: invalid prototype '%s' "
                    "(must start with 'bc1' or 'bc2')\n", seq->seq.s);
            seq_is_valid = false;
        }

        while (segment != NULL) {
            if (! str_in_list(segment, valid_names, num_valid_names) || num_segments > 5) {
                fprintf(stderr, "Error: unrecognized segment '%s' in prototype '%s'\n",
                        segment, seq->seq.s);
                seq_is_valid = false;
                break;
            }

            segment = strtok(NULL, "|");
            num_segments++;
        }
    }

    return seq_is_valid;
}


size_t allele_char_to_enum(char allele)
{
    size_t allele_i = 0;

    switch (allele) {
        case 'A':
            allele_i = ALLELE_A; break;
        case 'C':
            allele_i = ALLELE_C; break;
        case 'G':
            allele_i = ALLELE_G; break;
        case 'T':
            allele_i = ALLELE_T; break;
    }

    return allele_i;
}


bool all_standard_barcodes(const library_seqs *fs2_seqs) {
    size_t num_standard_bc = sizeof(FS2_BARCODES) / sizeof(FS2_BARCODES[0]);

    for (size_t s = 0; s < 2; s++) {
        for (size_t i = 0; i < fs2_seqs->num_barcodes[s]; i++) {
            bool bc_is_standard = false;

            for (size_t j = 0; j < num_standard_bc; j++) {
                if (strcmp(fs2_seqs->barcodes[s][i].seq, FS2_BARCODES[j]) == 0) {
                    bc_is_standard = true;
                    break;
                }
            }

            if (! bc_is_standard) {
                return false;
            }
        }
    }

    return true;
}


library_seqs *load_fasta_sequences(const char *filepath)
{
    gzFile fp = gzopen(filepath, "r");

    if (fp == NULL) {
        fprintf(stderr, "Error: unable to read file '%s': %s\n",
                filepath, strerror(errno));
        return NULL;
    }

    kseq_t *seq = kseq_init(fp);
    bool error_occurred = false;
    bool allocation_error = false;
    size_t num_bc_slots[2] = {0};

    library_seqs *fs2_seqs = calloc(1, sizeof(*fs2_seqs));

    if (fs2_seqs == NULL) {
        allocation_error = true;
        goto CLEANUP;
    }

    while (kseq_read(seq) >= 0) {
        if (! valid_fasta_seq(seq)) {
            error_occurred = true;
            break;
        }

        if (seq->name.s[0] != 'p') {
            str_toupper(seq->seq.s);
        }

        char *copy_str = NULL;
        char *dest_str = NULL;

        if (strncmp(seq->name.s, "bc", 2) == 0) {
            int BC_1_OR_2 = atoi(seq->name.s + 2) - 1;

            if (fs2_seqs->num_barcodes[BC_1_OR_2] == num_bc_slots[BC_1_OR_2]) {
                num_bc_slots[BC_1_OR_2] += 48;

                void *alloc_tmp = realloc(fs2_seqs->barcodes[BC_1_OR_2],
                                          num_bc_slots[BC_1_OR_2] * sizeof(barcode_t));

                if (alloc_tmp == NULL) {
                    error_occurred = true;
                    allocation_error = true;
                    break;
                }

                fs2_seqs->barcodes[BC_1_OR_2] = alloc_tmp;
            }

            size_t bc_i = fs2_seqs->num_barcodes[BC_1_OR_2];
            copy_str = seq->seq.s;
            dest_str = fs2_seqs->barcodes[BC_1_OR_2][bc_i].seq;
            fs2_seqs->barcodes[BC_1_OR_2][bc_i].label = atoi(seq->comment.s);

            fs2_seqs->num_barcodes[BC_1_OR_2]++;
        }
        else if (strcmp(seq->name.s, "allele") == 0) {
            size_t allele_i = allele_char_to_enum(seq->seq.s[0]);
            copy_str = seq->comment.s;
            dest_str = fs2_seqs->alleles[allele_i].seq;
        }
        else {
            size_t last_char_index = strlen(seq->name.s) - 1;
            int seq_1_or_2 = atoi(seq->name.s + last_char_index) - 1;
            copy_str = seq->seq.s;

            switch (seq->name.s[0]) {
                case 'a':
                    dest_str = fs2_seqs->adapters[seq_1_or_2].seq; break;
                case 'f':
                    dest_str = fs2_seqs->flanking[seq_1_or_2].seq; break;
                case 'p':
                    dest_str = fs2_seqs->prototype_strings[seq_1_or_2].seq; break;
            }
        }

        strcpy(dest_str, copy_str);
    }

    if (num_fasta_seqs_read(fs2_seqs->adapters, 2) != 2 ||
        num_fasta_seqs_read(fs2_seqs->flanking, 2) != 2 ||
        num_fasta_seqs_read(fs2_seqs->prototype_strings, 2) != 2 ||
        num_fasta_seqs_read(fs2_seqs->alleles, 4) < 1) {

        fprintf(stderr, "Error: FASTA file must contain two adapter sequences, "
               "flanking sequences, and prototypes and at least one allele\n");
        error_occurred = true;
    }

  CLEANUP:
    if (allocation_error) {
        perror("Error: memory allocation failed while reading FASTA");
    }
    if (error_occurred) {
        free(fs2_seqs->barcodes[0]);
        free(fs2_seqs->barcodes[1]);
        free(fs2_seqs);
        fs2_seqs = NULL;
    }

    kseq_destroy(seq);
    gzclose(fp);

    return fs2_seqs;
}


void parse_prototypes(library_seqs *fs2_seqs)
{
    for (size_t i = 0; i < 2; i++) {
        char *prototype_str = fs2_seqs->prototype_strings[i].seq;
        char *segment = strtok(prototype_str, "|");
        size_t segment_length = 0;
        size_t segment_index = 0;
        size_t offset_counter = 0;

        while (segment != NULL) {
            if (strncmp(segment, "bc", 2) == 0) {
                segment_length = 6;
            }
            else if (strcmp(segment, "allele") == 0) {
                fs2_seqs->prototypes[i].allele_offset = offset_counter;
                segment_length = 1;
            }
            else {
                size_t last_char_index = strlen(segment) - 1;
                int seq_1_or_2 = atoi(segment + last_char_index) - 1;
                read_segment *segment_ptr = NULL;

                if (segment[0] == 'a') {
                    segment_ptr = &(fs2_seqs->adapters[seq_1_or_2]);
                }
                else if (segment[0] == 'f') {
                    segment_ptr = &(fs2_seqs->flanking[seq_1_or_2]);
                }

                segment_length = strlen(segment_ptr->seq);
                segment_ptr->offset = offset_counter;
                segment_ptr->length = segment_length;

                fs2_seqs->prototypes[i].segments[segment_index] = segment_ptr;

                segment_index++;
            }

            offset_counter += segment_length;
            segment = strtok(NULL, "|");
        }
    }
}
