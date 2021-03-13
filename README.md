# *fsdm* (FREQ-Seq<sup>2</sup> Demultiplexer)

FREQ-Seq<sup>2</sup> is a method for quantifying the relative frequencies of different alleles at loci of interest in mixed population samples using customized double-barcoded libraries. There are 48 different sequences available for each of the two barcodes, enabling 48<sup>2</sup> (2,304) different possible sample barcoding configurations per library.

This repository contains the accompanying software utility, *fsdm*, which processes raw the FASTQ sequencing data of a FREQ-Seq<sup>2</sup> library and outputs the demultiplexed coverage for each sample in the library. *fsdm* is optimized for efficient demultiplexing of massive quantities of NGS sequencing data while using very little memory. As such, it is suitable for use on both HPC systems and personal computers.

## Install

```
git clone https://github.com/rnzhao/fsdm.git
cd fsdm
make
```

*fsdm* requires only a C99 compiler and zlib to build. It should compile on most Mac or Linux systems without needing to install additional software or dependencies.

## Usage

*fsdm* requires a FASTA file containing the library sequences and one or more pairs of FASTQ files. You can optionally specify a maximum permitted number of mismatches in the barcode, adapter, and flanking sequences, as well as a maximum permitted cumulative edit distance for each read pair.

The included FASTA file contains all 48 FREQ-Seq<sup>2</sup> barcodes as well as placeholders for providing *fsdm* with your specific library sequences. The alleles and flanking sequences should be adjusted for each library. Unused barcodes can be removed from the FASTA if you don't want to include them in the program's output.

For an overview of the usage and command line options, run `fsdm -h`.

## License

Mozilla Public License (MPL) 2.0
