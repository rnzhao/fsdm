/*
   Copyright (c) 2020 Roy Zhao <roy.zhao@uci.edu>
   Mozilla Public License Version 2.0
*/

#include "args.h"
#include "bc_hash.h"
#include "demultiplex.h"
#include "edit_distance.h"
#include "fs2_barcodes.h"
#include "parse_seq.h"

#include <errno.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>


int main(int argc, const char **argv)
{
    args args = parse_args(argc, argv);
    library_seqs *fasta_seqs = load_fasta_sequences(args.fasta_file);

    if (fasta_seqs == NULL) {
        return EXIT_FAILURE;
    }

    parse_prototypes(fasta_seqs);

    unsigned int num_bc[2];
    unsigned int total_num_unique_barcodes;
    unsigned int num_standard_barcodes = sizeof(FS2_BARCODES) / sizeof(FS2_BARCODES[0]);

    if (args.output_all && ! all_standard_barcodes(fasta_seqs)) {
        args.output_all = false;
    }

    if (args.output_all) {
        num_bc[0] = num_standard_barcodes;
        num_bc[1] = num_standard_barcodes;
        total_num_unique_barcodes = num_standard_barcodes;

        puts("Outputting all barcode combinations ('-a' option). "
             "Refer to the standard barcode number labels from 1-48.");
    }
    else {
        num_bc[0] = fasta_seqs->num_barcodes[0];
        num_bc[1] = fasta_seqs->num_barcodes[1];

        if (all_standard_barcodes(fasta_seqs)) {
            total_num_unique_barcodes = num_standard_barcodes;
        }
        else {
            total_num_unique_barcodes = num_bc[0] + num_bc[1];
        }
    }

    bc_counter *counter = calloc(1, sizeof(*counter) + num_bc[0] *
                                    num_bc[1] * sizeof(*counter->counts));
    counter->num_bc1 = num_bc[0];
    counter->num_bc2 = num_bc[1];

    size_t num_slots = calc_num_combos(6, total_num_unique_barcodes, args.bc_mismatches);
    bc_hash_table *hash_table = init_hash_table(num_slots);

    if (args.bc_mismatches > 0) {
        unsigned long num_possible_perms = uint_pow(4, 6);
        char (*all_permutations)[7] = calloc(num_possible_perms, 7);

        generate_all_bc_combos(6, all_permutations);

        for (size_t i = 0; i < num_possible_perms; i++) {
            for (size_t j = 0; j < 2; j++) {
                if (args.output_all && j > 0) {
                    break;
                }

                const char *mm_bc = NULL;
                size_t mm_bc_index = 0;

                for (size_t k = 0; k < num_bc[j]; k++) {
                    const char *current_bc;

                    if (args.output_all) {
                        current_bc = FS2_BARCODES[k];
                    }
                    else {
                        current_bc = fasta_seqs->barcodes[j][k].seq;
                    }

                    int mismatches = hamming_distance(current_bc, all_permutations[i], 6);

                    if (mismatches <= args.bc_mismatches) {
                        if (mm_bc) {
                            mm_bc = NULL;
                            break;
                        }
                        else {
                            mm_bc = current_bc;
                            mm_bc_index = k;
                        }
                    }
                }

                if (mm_bc) {
                    if (args.output_all) {
                        hash_table_insert(hash_table, mm_bc, mm_bc_index + 1, 0, false);
                        hash_table_insert(hash_table, mm_bc, mm_bc_index + 1, 1, false);
                    }
                    else {
                        hash_table_insert(hash_table, mm_bc, mm_bc_index + 1, j, false);
                    }
                }
            }
        }
    }

    if (args.output_all) {
        for (size_t i = 0; i < num_standard_barcodes; i++) {
            hash_table_insert(hash_table, FS2_BARCODES[i], i + 1, 0, true);
            hash_table_insert(hash_table, FS2_BARCODES[i], i + 1, 1, true);
        }
    }
    else {
        for (size_t i = 0; i < 2; i++) {
            for (size_t j = 0; j < num_bc[i]; j++) {
                hash_table_insert(hash_table, fasta_seqs->barcodes[i][j].seq, j + 1, i, true);
            }
        }
    }

    prune_hash_table(&hash_table);

    bool valid_alleles[4] = {false};

    for (size_t i = 0; i < 4; i++) {
        char *allele_name = fasta_seqs->alleles[i].seq;

        if (*allele_name != '\0') {
            valid_alleles[i] = true;
        }
    }

    for (int i = 0; i < args.num_fastq_pairs; i++) {
        demultiplex_fastq_pair(args.fastq_files + 2 * i, fasta_seqs, hash_table, counter,
                               args.ad_fl_mismatches, args.ed_threshold, valid_alleles);
    }

    FILE *output_fp = NULL;

    if (args.outfile) {
        output_fp = fopen(args.outfile, "w");
    }
    else {
        output_fp = stdout;
    }

    fprintf(output_fp, "bc1\tbc2");

    for (size_t i = 0; i < 4; i++) {
        if (valid_alleles[i]) {
            fprintf(output_fp, "\t%s", fasta_seqs->alleles[i].seq);
        }
    }
    fprintf(output_fp, "\n");

    for (size_t i = 0; i < num_bc[0]; i++) {
        for (size_t j = 0; j < num_bc[1]; j++) {
            int bc1_label, bc2_label;

            if (args.output_all) {
                bc1_label = i + 1;
                bc2_label = j + 1;
            }
            else {
                bc1_label = fasta_seqs->barcodes[0][i].label;
                bc2_label = fasta_seqs->barcodes[1][j].label;
            }

            fprintf(output_fp, "%d\t%d", bc1_label, bc2_label);

            for (size_t a = 0; a < 4; a++) {
                if (valid_alleles[a]) {
                    fprintf(output_fp, "\t%d", counter->counts[num_bc[1] * i + j][a]);
                }
            }
            fprintf(output_fp, "\n");
        }
    }

    if (args.outfile) {
        fclose(output_fp);
    }

    return 0;
}
