#!/usr/bin/env perl 

use warnings;
use strict;
use Data::Dumper;
use Getopt::Long;
use File::Basename qw/basename/;
use File::Copy qw/mv cp/;
use File::Temp qw/tempdir/;
use File::Which qw/which/;

local $0 = basename $0;
sub logmsg{local $0=basename $0; print STDERR "$0: @_\n";}
exit(main());

sub main{
  my $settings={};
  GetOptions($settings,qw(help outdir=s tempdir=s)) or die $!;
  die usage() if($$settings{help} || @ARGV < 2);
  $$settings{tempdir} //= tempdir("$0.XXXXXX", CLEANUP=>1, TMPDIR=>1);
  $$settings{outdir} //= "./$0.out";

  for my $exec(qw(cutadapt bowtie2-build samtools bcftools tabix bgzip)){
    my $path = which($exec);
    if(!$path){
      die "ERROR: could not find $exec in PATH";
    }
  }

  mkdir $$settings{outdir};
  my $ref = shift(@ARGV);
  my %result;
  while(@ARGV){
    my $R1 = shift(@ARGV);
    my $R2 = $R1;
       $R2 =~ s/_1\.fastq.gz$/_2.fastq.gz/;
       # substitute R2
    if(!-e $R2){
      $R2 = $R1;
      $R2 =~ s/_R1_/_R2_/;
    }
    if(!-e $R2 || $R1 eq $R2){
      logmsg "SKIP: could not find R2 for $R1";
      next;
    }

    logmsg "Workflow on $R1 / $R2";
    if(-e "$$settings{outdir}/".basename($R1).".fasta"){
      logmsg "SKIP: found $$settings{outdir}/".basename($R1).".fasta";
      next;
    }
    my $res = singleSampleWorkflow($R1, $R2, $ref, $settings);
    $result{$R1} = $res;

    logmsg "vcf:   $$res{vcf}";
    logmsg "fasta: $$res{fasta}";
  }

  return 0;
}

sub singleSampleWorkflow{
  my($R1, $R2, $ref, $settings) = @_;

  my ($trimmedR1, $trimmedR2) = adapterTrim($R1, $R2, $settings);
  my $bam   = mapReads($trimmedR1, $trimmedR2, $ref, $settings);
  my $vcf   = snps($bam, $ref, $settings);
  my $fasta = consensus($vcf, $ref, $settings);

  my $outfasta = "$$settings{outdir}/".basename($fasta);
  my $outvcf   = "$$settings{outdir}/".basename($vcf);
  cp($fasta, $outfasta) or die "ERROR copying $fasta => $outfasta: $!";
  cp($vcf  , $outvcf) or die "ERROR copying $vcf => $outvcf: $!";

  return {vcf=>$outvcf, fasta=>$outfasta};
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

  if(-e "$bam.bai"){
    logmsg "SKIP: found $bam.bai";
    return $bam;
  }

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
  my $gcf = $vcf; # vcf.gz
     $gcf =~ s/vcf$/vcf.gz/;

  if(-e "$gcf.tbi"){
    logmsg "SKIP: found $gcf.tbi";
    return $gcf;
  }

  system("samtools mpileup -aa -d 8000 -uf $ref $bam | bcftools call -Mc > $vcf");
  die if $?;

  system("bgzip $vcf"); # creates vcf.gz
  die if $?;
  system("tabix -f $gcf"); # Creates vcf.gz.tbi
  die if $?;

  return $gcf;
}

sub consensus{
  my($vcf, $ref, $settings) = @_;

  my $consensusfasta = "$$settings{tempdir}/".basename($vcf,".vcf.gz").".fasta";

  if(-e $consensusfasta){
    logmsg "SKIP: found $consensusfasta";
    return $consensusfasta;
  }

  system("zcat $vcf | vcfutils.pl vcf2fq -d 100 -D 100000000 | seqtk seq -A - | sed '2~2s/[actg]/N/g' > $consensusfasta.tmp");

  mv("$consensusfasta.tmp", $consensusfasta);

  return $consensusfasta;
}

sub usage{
  "$0: runs Clint's SNP workflow on raw illumina reads
  Usage: $0 [options] ref.fasta *R1*.fastq.gz
  R2s are found by replacing _1.fastq.gz with _2.fastq.gz
                          or _R1_ with _R2_
  --help   This useful help menu
  --tempdir ''
  --outdir  ''
  "
}
