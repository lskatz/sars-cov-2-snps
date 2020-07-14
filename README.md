# sars-cov-2-snps
Pipeline for identifying SNPs in SARS-CoV-2

Uses Paden et al 2020 EID for Illumina PE workflow at
https://github.com/CDCgov/SARS-CoV-2_Sequencing/tree/master/protocols/CDC-Comprehensive

## Usage

    Usage: sars-cov-2-snps.pl [options] ref.fasta *R1*.fastq.gz
    R2s are found by replacing _1.fastq.gz with _2.fastq.gz
                            or _R1_ with _R2_
    --help   This useful help menu
    --tempdir ''
    --outdir  ''

