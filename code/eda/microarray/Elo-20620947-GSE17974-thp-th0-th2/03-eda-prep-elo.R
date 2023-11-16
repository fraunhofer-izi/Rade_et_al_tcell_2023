.cran_packages = c("devtools", "dplyr", "yaml", "ff")
.bioc_packages = c("affy", "simpleaffy", "AnnotationDbi", "vsn")

# http://brainarray.mbni.med.umich.edu/Brainarray/Database/CustomCDF/23.0.0/gencodeg.asp
.url_packages = c("hgu133plus2hsgencodegcdf", "hgu133plus2hsgencodegprobe",
                  "pd.hgu133plus2.hs.gencodeg")
.inst = .url_packages %in% installed.packages()
if (any(!.inst)) {
  devtools::install_url("http://mbni.org/customcdf/23.0.0/gencodeg.download/hgu133plus2hsgencodegcdf_23.0.0.tar.gz")
  devtools::install_url("http://mbni.org/customcdf/23.0.0/gencodeg.download/hgu133plus2hsgencodegprobe_23.0.0.tar.gz")
  devtools::install_url("http://mbni.org/customcdf/23.0.0/gencodeg.download/pd.hgu133plus2.hs.gencodeg_23.0.0.tar.gz")
}

## Install CRAN packages (if not already installed)
.inst = .cran_packages %in% installed.packages()
if (any(!.inst)) {
  install.packages(.cran_packages[!.inst], repos = "http://cran.rstudio.com/")
}

# Install bioconductor packages (if not already installed)
.inst <- .bioc_packages %in% installed.packages()
if (any(!.inst)) {
  library(BiocManager)
  BiocManager::install(.bioc_packages[!.inst], ask = F)
}

list.of.packages = c(.cran_packages, .bioc_packages, .url_packages)
# list.of.packages = c(.cran_packages, .bioc_packages, .url_packages)

## Loading library
for (pack in list.of.packages) {
  suppressMessages(library(
    pack,
    quietly = TRUE,
    verbose = FALSE,
    character.only = TRUE
  ))
}

# # Wenn folgende Lib "pd.hgu133plus2.hs.gencodeg" nicht laden will:
# # install.packages("https://cran.r-project.org/src/contrib/Archive/ff/ff_2.2-14.tar.gz", repos=NULL)
# # https://community.rstudio.com/t/unable-to-install-bioconductor-package/75223
#
## >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
## MANIFEST
## >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
manifest = yaml.load_file("assets/manifest.yaml")
workdata = manifest$workdata
rawdata = manifest$rawdata

cel.dir = paste0(rawdata, manifest$rawdata_elo_thp_th0_th2)
meta.table.dir = manifest$phenodata_elo_thp_th0_th2
out.file = paste0(workdata, manifest$eda_prep_elo_thp_th2)
out.affybatch = paste0(workdata, manifest$affybatch_elo_thp_th2)

## >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
## Load GC.v29 annotation (GENE_NAME, GENE_TYPE, GENE_TYPE_CLUSTER)
## Clustered subbiotypes based on
## ftp://ftp.ebi.ac.uk/pub/databases/gencode/_README_stats.txt
## >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
gc.v29.anno = readRDS(paste0(workdata, manifest$gencode_v29_features))

## >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
## CEL FILENAMES
## >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
cel.files = list.files(cel.dir)
cel.files.df = data.frame(FILE_NAME = cel.files, SAMPLE_NAME = gsub(".CEL.gz", "", cel.files))

meta.table = readRDS(paste0(workdata, meta.table.dir, "elo-thp-th0-th2.Rds"))[[1]]
meta.table$FILE_NAME = cel.files.df$FILE_NAME[match(meta.table$SAMPLE_NAME, cel.files.df$SAMPLE_NAME)]
rownames(meta.table) = meta.table$FILE_NAME

meta.table = meta.table[!grepl("th0.*", meta.table$SHORT_NAME), ]
# meta.table = meta.table[!grepl("th2_4h$", meta.table$SHORT_NAME), ]
meta.table = droplevels(meta.table)
meta.table = AnnotatedDataFrame(meta.table)

## >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
## Read cel files with affy package & custon CDF (BRAINARRAY)
## http://brainarray.mbni.med.umich.edu/Brainarray/Database/CustomCDF/23.0.0/gencodeg.asp
## >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
affybatch = NULL
affybatch = affy::ReadAffy(verbose=TRUE,
                    filenames = file.path(cel.dir, meta.table$FILE_NAME),
                    phenoData = meta.table,
                    sampleNames = as.character(meta.table$SAMPLE_ID),
                    cdfname="hgu133plus2hsgencodeg",
                    compress = TRUE)

## >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
## simpleaffy
## >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# n = cdfName(affybatch)
# cleancdfname(n)
setQCEnvironment("hgu133plus2hsgencodegcdf", "assets/")
simpleAffy.qc = qc(affybatch)

## >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
## For RNA degradation plot
## >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
rna.deg = affy::AffyRNAdeg(affybatch)

## >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
## RMA
## >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
eset.rma = affy::rma(affybatch)
eset.vsn = vsnrma(affybatch)
eset.rma.raw = affy::rma(affybatch, background = FALSE, normalize = FALSE)

esets = list(rma = eset.rma, vsn = eset.vsn)

## >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
## Annotation of the ENSEMBL IDs -> DF with ENSEMBL ID, gene symbol and biotype
## >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
for (i in names(esets)) {
  eset = NULL
  eset.ftr = NULL
  eset = esets[[i]]

  eset.ftr = data.frame(gsub("\\_.*", "", rownames(exprs(eset))),
                            row.names = rownames(exprs(eset)),
                            stringsAsFactors = F)
  colnames(eset.ftr) = "PROBE_ID"

  rows = rownames(eset.ftr)
  eset.ftr = merge(x = eset.ftr, y = gc.v29.anno,
                       by.x = "PROBE_ID", by.y = "ENSEMBL_ID",
                       all.x = TRUE)

  # nrow(eset.ftr[!duplicated(eset.ftr$PROBE_ID), ]) == nrow(eset.ftr)
  colnames(eset.ftr) = c("ENSEMBL_ID", "CHR", "GENE_SYMBOL", "BIOTYPE", "BIOTYPE_CLUSTER",
                             "ENSEMBL_ID_ABBR", "GENE_SYMBOL_DUPL_MARKED", "ENTREZ", "DESCRIPTION")
  rownames(eset.ftr) = rows
  eset.ftr$PROBE_ID = rownames(eset.ftr)

  print(identical(rownames(eset.ftr), rownames(eset)))

  ## >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  ## Add features
  ## >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  fData(eset)$PROBE_ID = rownames(exprs(eset))
  fData(eset) = merge(x = fData(eset), y = eset.ftr,
                          by.x = "PROBE_ID", by.y = "row.names",
                          all.x = TRUE)
  fData(eset)$PROBE_ID.y = NULL
  # restore rownames after left_join
  rownames(fData(eset)) = fData(eset)$PROBE_ID

  print(identical(rownames(fData(eset)), fData(eset)$PROBE_ID))

  esets[[i]] = eset
}

obj = list(eset.rma.raw = eset.rma.raw,
           esets = esets,
           simpleAffy.qc = simpleAffy.qc,
           rna.deg = rna.deg)
saveRDS(obj, file = out.file)

saveRDS(affybatch, file = out.affybatch)
