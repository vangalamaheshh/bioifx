#!/bin/bash
# vim: syntax=sh tabstop=2 expandtab

genome_dir=$1
fasta_file=$2
gtf_file=$3
num_threads=$4

if [[ -z $genome_dir || -z $fasta_file || -z $gtf_file || -z $num_threads ]]; then
  read -r -d '' usage <<'EOF' 
    Usage: $0 <genome_dir> <fasta_file> <gtf_file> <num_threads>
EOF
  echo "$usage"
  exit 1
fi

echo -n "Start time: "
date
STAR --runThreadN $num_threads \
  --runMode genomeGenerate \
  --genomeDir $genome_dir \
  --genomeFastaFiles $fasta_file \
  --sjdbGTFfile $gtf_file \
  --sjdbOverhang 100
echo -n "End time: "
date

