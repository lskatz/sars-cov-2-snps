# sars-cov-2-snps
Pipeline for identifying SNPs in SARS-CoV-2

Uses Paden et al 2020 EID for Illumina PE workflow at
https://github.com/CDCgov/SARS-CoV-2_Sequencing/tree/master/protocols/CDC-Comprehensive

## Usage

Options include where to place your temporary directory
and your output directory.  Defaults are available for each
and so these options are optional.

    Usage: sars-cov-2-snps.pl [options] ref.fasta *R1*.fastq.gz
    R2s are found by replacing _1.fastq.gz with _2.fastq.gz
                            or _R1_ with _R2_
    --help   This useful help menu
    --tempdir ''
    --outdir  ''

## Workflow

For each R1 file:

1. The R2 file is found by simple replace from the R1 read
2. Cutadapt to trim adapters
3. Reads mapped to reference with bowtie2
4. SNPs identified with `samtools mpileup`
5. Consensus fasta created with `vcfutils.pl` and seqtk. Sed is used to replace lower quality bases with N.

