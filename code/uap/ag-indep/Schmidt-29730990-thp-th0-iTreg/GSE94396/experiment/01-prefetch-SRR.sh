#!/bin/bash

module load SRA-Toolkit/2.9.2-ubuntu64

## change in line 38 the path where files should be downloaded
export VDB_CONFIG=custom.kfg

## test: vdb-config /repository/user/main

parallel -j 10 prefetch {} ::: $(cat SRR-accession-list)

## The -j 1 specifies the number of threads to use. Using 1 limits to downloading one file at a time
## (simultaneous downloads may be faster, depending on your computer and network).

## The ::: $(cat SRR_Acc_List.txt) passes the contents of SRR_Acc_List.txt as arguments to the parallel
