#!/bin/bash

module load SRA-Toolkit/2.9.2-ubuntu64

parallel -j 30 fastq-dump --skip-technical -F --gzip --split-files -O /data/RNA-Seq/2018-MAVO-MCF/Schmidt-29730990-thp-th0-iTreg/GSE94396/Unaligned/Project_Schmidt-29730990-GSE94396-thp-th0-iTreg/ {} ::: $(cat SRA-manifest)
## A nice explanation of other fastq-dump options are provided by Rob Edward's group: https://edwards.sdsu.edu/research/fastq-dump/
