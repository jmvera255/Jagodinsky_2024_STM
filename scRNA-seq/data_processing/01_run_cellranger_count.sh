#! /usr/bin/env bash
# 01_run_cellranger_count.sh

# terminate if error occurs
set -e

# get extract software and ref data files
tar xzf refdata-gex-mm10-2020-A.tar.gz
tar xzf cellranger-6.1.2.tar.gz
export PATH=$PWD/cellranger-6.1.2:$PATH

cellranger count --disable-ui --no-bam --nosecondary --localcores=4 \
--id=<sample> --fastqs=<sample>_fastq --transcriptome=refdata-gex-mm10-2020-A