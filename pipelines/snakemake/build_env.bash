#!/bin/bash

echo y | conda update conda

echo y | conda install -c https://conda.anaconda.org/bioconda \
	python=3.5 \
	java-jdk=8 \
	pandas=0.17 \
	snakemake=3.4 \
	star=2.5 \
	star-fusion=0.5 \
	cufflinks=2.2 \
	samtools=1.2

echo y | conda install -c r r-essentials

