/*
   Copyright (c) 2020 Roy Zhao <roy.zhao@uci.edu>
   Mozilla Public License Version 2.0
*/

#include "args.h"
#include "argparse.h"

#include <errno.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>


args parse_args(int argc, const char **argv)
{
    static const char *usage[] = {
        "fsdm [options] <sequences.fa> <reads_1.fq> <reads_2.fq>",
        "(FASTQ files can be gzipped or uncompressed, and multiple pairs can be provided at once.)",
        NULL
    };

    static const char *description = "fsdm (FREQ-Seq^2 Demultiplexer) v1.0.1";
    static const bool SHOW_HELP_OPT = false;

    args parsed_args = {
        .num_fastq_pairs = 0,
        .outfile = NULL,
        .output_all = false,
        .bc_mismatches = 0,
        .ad_fl_mismatches = 1,
        .ed_threshold = 4
    };

    struct argparse_option arguments[] = {
        OPT_HELP(SHOW_HELP_OPT),

        OPT_GROUP("Options"),
        OPT_STRING('o', NULL, &parsed_args.outfile,
                   "Output file (results are printed to stdout if unspecified)",
                   NULL, 0, 0),
        OPT_BOOLEAN('a', NULL, &parsed_args.output_all,
                    "Output all possible barcode combinations",
                    NULL, 0, 0),
        OPT_INTEGER(0, "bm", &parsed_args.bc_mismatches,
                    "Number of mismatches allowed in a barcode sequence (default 0)",
                    NULL, 0, 0),
        OPT_INTEGER(0, "mm", &parsed_args.ad_fl_mismatches,
                    "Number of mismatches allowed in each adapter or flanking sequence (default 1)",
                    NULL, 0, 0),
        OPT_INTEGER(0, "ed", &parsed_args.ed_threshold,
                    "Maximum edit distance allowed across all adapter and flanking sequences (default 4)",
                    NULL, 0, 0),

        OPT_END()
    };

    struct argparse parser;
    argparse_init(&parser, arguments, usage, 0);
    argparse_describe(&parser, description, NULL);

    if (argc == 1) {
        argparse_usage(&parser, true);
        exit(EXIT_FAILURE);
    }

    argc = argparse_parse(&parser, argc, argv);

    if (argc < 3 || ((argc - 1) & 1) != 0) {
        fprintf(stderr, "Error: invalid number of FASTA/FASTQ files\n\n\n");
        argparse_usage(&parser, false);
        exit(EXIT_FAILURE);
    }

    bool argument_error = false;

    if (parsed_args.bc_mismatches < 0 || parsed_args.ad_fl_mismatches < 0 ||
        parsed_args.ed_threshold < 0) {
        fprintf(stderr, "Error: number of allowed mismatches cannot be negative\n");
        argument_error = true;
    }
    else if (parsed_args.bc_mismatches >= 6) {
        fprintf(stderr, "Error: number of barcode mismatches must be lower than barcode length\n");
        argument_error = true;
    }

    for (int i = 0; i < argc; i++) {
        const char *seq_file = argv[i];

        FILE *fp = fopen(seq_file, "r");

        if (fp == NULL) {
            fprintf(stderr, "Error: unable to read file '%s': %s\n",
                    seq_file, strerror(errno));
            argument_error = true;
        }

        fclose(fp);
    }

    if (parsed_args.outfile) {
        FILE *fp = fopen(parsed_args.outfile, "a");

        if (fp == NULL) {
            fprintf(stderr, "Error: unable to open output file '%s': %s\n",
                    parsed_args.outfile, strerror(errno));
            argument_error = true;
        }

        fclose(fp);
    }

    if (argument_error) {
        exit(EXIT_FAILURE);
    }

    parsed_args.fasta_file = argv[0];
    parsed_args.fastq_files = argv + 1;
    parsed_args.num_fastq_pairs = (argc - 1) / 2;

    return parsed_args;
}
