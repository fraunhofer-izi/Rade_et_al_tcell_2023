#!/bin/bash

module load SRA-Toolkit/2.9.2-ubuntu64

parallel -j 30 fastq-dump --skip-technical -F --gzip --split-files -O /data/RNA-Seq/2018-MAVO-MCF/Tuomela-26967054-GSE52260-thp-th0-th17/Unaligned/Project_Tuomela-26967054-GSE52260-thp-th0-th17/ {} ::: $(cat SRA-manifest)
## A nice explanation of other fastq-dump options are provided by Rob Edward's group: https://edwards.sdsu.edu/research/fastq-dump/
