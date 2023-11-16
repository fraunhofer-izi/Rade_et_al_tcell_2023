# Rade_et_al_tcell_2023

This repository contains code used to produce the results in: Michael Rade, Sebastian Böhlen, Vanessa Neuhaus, Dennis Löffler, Conny Blumert, Ulrike Köhl, Susann Dehmel, Katherina Sewald, and Kristin Reiche, A time-resolved meta-analysis of consensus gene expression profiles during human T-cell activation, May 2023, PREPRINT (Version 1) available at [biorxiv](https://www.biorxiv.org/content/10.1101/2023.05.03.538418v1)

``` sh
## >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
## Prepare microarray data (download from NCBI GEO, normalize, custom CDF etc.)
## >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# thp, th1, Aijö et al.
Rscript code/eda/microarray/Aijö-23110343-GSE32959-thp-th0-th1-th2/01-download-rawdata-aijö.R
Rscript code/eda/microarray/Aijö-23110343-GSE32959-thp-th0-th1-th2/02-phenoData-aijö.R
Rscript code/eda/microarray/Aijö-23110343-GSE32959-thp-th0-th1-th2/03-eda-prep-aijö.R
Rscript code/eda/microarray/Aijö-23110343-GSE32959-thp-th0-th1-th2/04-plm-aijö.R

# thp, th0, th2, Elo et al.
Rscript code/eda/microarray/Elo-20620947-GSE17974-thp-th0-th2/01-download-rawdata-elo.R
Rscript code/eda/microarray/Elo-20620947-GSE17974-thp-th0-th2/02-phenoData-elo.R
Rscript code/eda/microarray/Elo-20620947-GSE17974-thp-th0-th2/03-eda-prep-elo.R
Rscript code/eda/microarray/Elo-20620947-GSE17974-thp-th0-th2/04-plm-elo.R

## >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
## Prepare RNA-Seq data
## >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# thp, th0, treg (Schmidt)
# download SRA: code/uap/ag-indep/Schmidt-29730990-thp-th0-iTreg/GSE94396/experiment
# download SRA: code/uap/ag-indep/Schmidt-29730990-thp-th0-iTreg/GSE96538/experiment
# uap: code/uap/ag-indep/Schmidt-29730990-thp-th0-iTreg/GSE94396/pipeline
# uap: code/uap/ag-indepSchmidt-29730990-thp-th0-iTreg/GSE96538/pipeline

# thp, th0, itreg (Ullah)
# download SRA: code/uap/ag-indep/Ullah-29466736-GSE90569-thp-th0-iTreg/experiment
# uap: code/uap/ag-indep/Ullah-29466736-GSE90569-thp-th0-iTreg/pipeline

# thp, th17 (Tuomela)
# download SRA: code/uap/ag-indep/Tuomela-26967054-GSE52260-thp-th0-th17/experiment
# uap: code/uap/ag-indep/Tuomela-26967054-GSE52260-thp-th0-th17/pipeline

# thp, th0, (Verificaton 1)
# download gene counts: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi

# thp, th0, medium, (Verificaton 2)
# Download data: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE197067
Rscript code/eda/rna-seq/intern-datasets/thp-th0/01_counts_vst.R

# Discovery (RNA-Seq): One object for Schmidt, Ullah and Tuomela
Rscript code/eda/rna-seq/public-datasets/01_HTSeqCounts.R
Rscript code/eda/rna-seq/public-datasets/02_phenoData.R
Rscript code/eda/rna-seq/public-datasets/03_eda-prep.R

## >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
## DGEA
## >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# thp-th2 (Elo)
Rscript code/eda/microarray/Elo-20620947-GSE17974-thp-th0-th2/05-limma-elo.R
# thp-th1 (Aijö)
Rscript code/eda/microarray/Aijö-23110343-GSE32959-thp-th0-th1-th2/05-limma-aijö.R
# thp-th0 ()
Rscript code/eda/rna-seq/public-datasets/dgea_counts/thp-th0/eda-limma-thp-th0.R
# thp-itreg (Ullah)
Rscript code/eda/rna-seq/public-datasets/dgea_counts/thp-itreg/eda-limma-thp-itreg.R
# thp-th17 (Tuomela)
Rscript code/eda/rna-seq/public-datasets/dgea_counts/thp-th17/eda-limma-thp-th17.R
# thp-th0 (Arcelus, Verificaton 1)
Rscript code/eda/rna-seq/public-datasets/dgea_counts/cd4Mem-th0/limma_arcelus.R
# thp-th0, (Verificaton 2)
Rscript code/eda/rna-seq/intern-datasets/thp-th0/03_limma_act_vs_thp.R
Rscript code/eda/rna-seq/intern-datasets/thp-th0/04_limma_ctrl.R

## >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
## Meta-Analysis
## >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# One object for all datasets
Rscript code/eda/meta-effect/01_prep-meta.R FALSE FALSE
Rscript code/edapc/meta-effect/01_prep-meta.R TRUE FALSE
Rscript code/eda/meta-effect/01_prep-meta.R TRUE TRUE

# Combined effect size for data from the Discovery Set
Rscript code/eda/meta-effect/02_effect-size.R
# Combined effect sizes for all data (Discovery and Verification Sets)
Rscript code/eda/meta-effect/02_effect-size.R TRUE

## >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
## NMF
## >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# NMF; Discovery
Rscript code/nmf/nmf_array_aijoe.R 'vp18' 200 2:10 1
Rscript code/nmf/nmf_array_elo.R 'vp18' 200 2:10 1
Rscript code/nmf/nmf_rnaseq_meta.R 'vp18' 200 2:10 1

# Metagenes from discovery
Rscript code/nmf/pipeline_discovery_01.R

# NMF; Verification
Rscript code/nmf/nmf_rnaseq_izi.R 'vp18' 200 2:10 1
Rscript code/nmf/nmf_rnaseq_izi_all.R 'vp15' 200 2:10 1
Rscript code/nmf/nmf_rnaseq_arcelus.R 'vp18' 200 2:10 1
Rscript code/nmf/pipeline_verification_02.R

```
