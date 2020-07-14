#!/bin/bash -l

set -e

# Read ARGV
OUTDIR=$1
REF=$2
R1s=$3
SLOTS_PER_JOB=1 # manually change this as needed

if [ "$R1s" == "" ]; then
  scriptName=$(basename $0)
  echo "Usage: $scriptName outdir reference.fasta fofn.txt"
  echo "  Only supply the R1s for this workflow in the file of filenames."
  exit 1;
fi

if [ -f $OUTDIR ]; then
  echo "ERROR: $OUTDIR is not a directory";
  exit 1;
fi
mkdir -pv $OUTDIR

TMP=$(mktemp --tmpdir='.' --directory sars-snps.XXXXXXXX)
mkdir -p $TMP/log
echo "tmp dir is $TMP "

# CTRL file will have one SRA run ID per line
CTRL_FILE="$TMP/array.txt"
cat $R1s | perl -lane '
  for my $file(@F){
    print $file;
  }
' > $CTRL_FILE
echo "CTRL_FILE is $CTRL_FILE"


qsub -N sarsSNPs -o $TMP/log -j y -pe smp $SLOTS_PER_JOB -V -cwd -t 1-$(cat $CTRL_FILE | wc -l) \
  -v "REF=$REF" -v "OUTDIR=$OUTDIR" -v "CTRL_FILE=$CTRL_FILE" <<- "END_OF_SCRIPT"
  #!/bin/bash -l
  set -e
  source /etc/profile.d/modules.sh
  module purge
  module load bcftools samtools Python/3.7 cutadapt bowtie2
  
  # Set up filenames
  INFILE=$(sed -n ${SGE_TASK_ID}p $CTRL_FILE | awk '{print $1}')
  tmpdir=/scratch/$USER
  mkdir -p $tmpdir
  outdir=$(mktemp --tmpdir=$tmpdir --directory sars-snps.XXXXXX);
  trap "rm -rf $outdir" EXIT

  sars-cov-2-snps.pl --tempdir $outdir -o $OUTDIR $REF $INFILE

END_OF_SCRIPT

