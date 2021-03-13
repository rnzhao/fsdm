/*
   Copyright (c) 2020 Roy Zhao <roy.zhao@uci.edu>
   Mozilla Public License Version 2.0
*/

#include "demultiplex.h"

#include "bc_hash.h"
#include "edit_distance.h"
#include "parse_seq.h"
#include "kseq.h"

#include <errno.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <zlib.h>

KSEQ_INIT(gzFile, gzread)


static inline bool read_fastq_pair(kseq_t *fq_1,
                                   kseq_t *fq_2,
                                   int status[2])
{
    status[0] = kseq_read(fq_1);
    status[1] = kseq_read(fq_2);

    return (status[0] >= 0) && (status[1] >= 0);
}


void demultiplex_fastq_pair(const char **fastq_pair,
                            library_seqs *fs2_seqs,
                            bc_hash_table *hash_table,
                            bc_counter *bc_combo_counts,
                            int ad_fl_mismatches,
                            int ed_threshold,
                            const bool *valid_alleles)
{
    gzFile fastq_fp[2] = {NULL};

    for (size_t i = 0; i < 2; i++) {
        fastq_fp[i] = gzopen(fastq_pair[i], "r");

        if (fastq_fp[i] == NULL) {
            fprintf(stderr, "Error: unable to read file '%s': %s\n",
                    fastq_pair[i], strerror(errno));
            exit(EXIT_FAILURE);
        }
    }

    kseq_t *fq[2] = {kseq_init(fastq_fp[0]), kseq_init(fastq_fp[1])};
    int read_status[2] = {0};

    while (read_fastq_pair(fq[0], fq[1], read_status)) {
        int bc1 = hash_table_lookup(hash_table, fq[0]->seq.s, 0) - 1;
        int bc2 = hash_table_lookup(hash_table, fq[1]->seq.s, 1) - 1;

        if ((bc1 | bc2) < 0) {
            continue;
        }

        int edit_distance = 0;

        for (size_t i = 0; i < 2; i++) {
            size_t segment_index = 0;
            read_segment **segments = fs2_seqs->prototypes[i].segments;

            while (segments[segment_index] && edit_distance <= ed_threshold) {
                read_segment *segment = segments[segment_index];

                int segment_ed = damerau_levenshtein(segment->seq,
                                                     fq[i]->seq.s + segment->offset,
                                                     segment->length);

                if (segment_ed <= ad_fl_mismatches) {
                    edit_distance += segment_ed;
                }
                else {
                    edit_distance = ed_threshold + 1;
                }

                segment_index++;
            }
        }

        if (edit_distance > ed_threshold) {
            continue;
        }

        char allele = fq[0]->seq.s[fs2_seqs->prototypes[0].allele_offset];
        size_t allele_i = allele_char_to_enum(allele);

        if (! valid_alleles[allele_i]) {
            read_segment *left_flanking = &(fs2_seqs->flanking[0]);

            int allele_offset = nw_offset(left_flanking->seq,
                                          fq[0]->seq.s + left_flanking->offset,
                                          left_flanking->length);
            allele_offset += (int) fs2_seqs->prototypes[0].allele_offset;

            if (0 < allele_offset < fq[0]->seq.l) {
                allele_i = allele_char_to_enum(fq[0]->seq.s[allele_offset]);
            }
        }

        bc_combo_counts->counts[bc_combo_counts->num_bc2 * bc1 + bc2][allele_i] += 1;
    }

    kseq_destroy(fq[0]);
    kseq_destroy(fq[1]);

    if ((read_status[0] & read_status[1]) != -1) {
        fprintf(stderr, "Warning: Files in FASTQ pair have different number "
                "of reads: '%s', '%s'\n", fastq_pair[0], fastq_pair[1]);
    }
}
