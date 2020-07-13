#!/usr/bin/env perl 

use warnings;
use strict;
use Data::Dumper;
use Getopt::Long;
use File::Basename qw/basename/;
use File::Copy qw/mv cp/;
use File::Temp qw/tempdir/;

local $0 = basename $0;
sub logmsg{local $0=basename $0; print STDERR "$0: @_\n";}
exit(main());

sub main{
  my $settings={};
  GetOptions($settings,qw(help tempdir=s)) or die $!;
  die usage() if($$settings{help} || @ARGV < 3);
  $$settings{tempdir} //= tempdir("$0.XXXXXX", CLEANUP=>1, TMPDIR=>1);
  $$settings{outdir} //= "./$0.out";

  my $ref = shift(@ARGV);
  my %fasta;
  while(@ARGV){
    my($R1,$R2) = splice(@ARGV,0,2);
    my $fasta = singleSampleWorkflow($R1, $R2, $ref, $settings);
    $fasta{$R1} = $fasta;
  }

  die Dumper \%fasta;

  return 0;
}

sub singleSampleWorkflow{
  my($R1, $R2, $ref, $settings) = @_;

  my ($trimmedR1, $trimmedR2) = adapterTrim($R1, $R2, $settings);
  my $bam = mapReads($trimmedR1, $trimmedR2, $ref, $settings);
  my $fasta = snps($bam, $ref, $settings);

  my $outfasta = "$$settings{outdir}/".basename($fasta);
  cp($fasta, $outfasta);

  return $outfasta;
}

sub adapterTrim{
  my($read1, $read2, $settings) = @_;
  my $threads = 1;

  my $R1out = "$$settings{tempdir}/".basename($read1);
  my $R2out = "$$settings{tempdir}/".basename($read2);

  my $stdout = qx(cutadapt -j $threads -g GTTTCCCAGTCACGATA -G GTTTCCCAGTCACGATA -a TATCGTGACTGGGAAAC -A TATCGTGACTGGGAAAC -g ACACTCTTTCCCTACACGACGCTCTTCCGATCT -G ACACTCTTTCCCTACACGACGCTCTTCCGATCT -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT -n 3 -m 75 -q 25 --interleaved $read1 $read2 | cutadapt -j $threads --interleaved -m 75 -u 30 -u -30 -U 30 -U -30 -o $R1out -p $R2out -);

  return($R1out, $R2out);
}

sub mapReads{
  my($R1, $R2, $ref, $settings) = @_;
  my $threads = 1;
  my $sam = "$$settings{tempdir}/".basename($R1).".sam";
  my $bam = $sam;
     $bam =~ s/\.sam$/.bam/;

  my $index = "$$settings{tempdir}/bowtieindex";
  if(! -e "$index.1.bt2"){
    system("bowtie2-build $ref $index");
    die if $?;
  }

  system("bowtie2 --sensitive-local -p $threads -x $index  -1 $R1 -2 $R2 -S $sam"); die if $?;
  system("samtools view -b $sam | samtools sort - -o $bam"); die if $?;
  system("samtools index $bam"); die if $?;

  unlink($sam);
  return $bam;
}

sub snps{
  my($bam, $ref, $settings) = @_;

  my $vcf = "$$settings{tempdir}/".basename($bam,".bam").".vcf";
  my $consensusfasta = "$$settings{tempdir}/".basename($bam,".bam").".fasta";

  system("samtools mpileup -aa -d 8000 -uf $ref $bam | bcftools call -Mc |tee -a $vcf | vcfutils.pl vcf2fq -d 100 -D 100000000 | seqtk seq -A - | sed '2~2s/[actg]/N/g' > $consensusfasta");

  return $consensusfasta;
}

sub usage{
  "$0: runs Clint's SNP workflow on raw illumina reads
  Usage: $0 [options] ref.fasta R1.fastq.gz R2.fastq.gz [R1.fastq.gz R2.fastq.gz...]
  --help   This useful help menu
  --tempdir ''
  --outdir  ''
  "
}
