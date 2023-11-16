.cran_packages = c("dplyr", "yaml")
.bioc_packages = c("Biobase", "GEOquery")

## Install CRAN packages (if not already installed)
.inst = .cran_packages %in% installed.packages()
if (any(!.inst)) {
  install.packages(.cran_packages[!.inst], repos = "http://cran.rstudio.com/")
}

# Install bioconductor packages (if not already installed)
.inst <- .bioc_packages %in% installed.packages()
if (any(!.inst)) {
  source("http://bioconductor.org/biocLite.R")
  biocLite(.bioc_packages[!.inst], ask = F)
}

list.of.packages = c(.cran_packages, .bioc_packages)

## Loading library
for (pack in list.of.packages) {
  suppressMessages(library(
    pack,
    quietly = TRUE,
    verbose = FALSE,
    character.only = TRUE
  ))
}

## >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
## MANIFEST
## >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
manifest = yaml.load_file("assets/manifest.yaml")
workdata     = manifest$workdata

out.dir = paste0(workdata, manifest$phenodata_aijö_thp_th0_th1_th2)
out.file = paste0(out.dir, "/aijö-thp-th0-th1-th2.Rds")

## >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
##
## >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
gse = getGEO("GSE32959", GSEMatrix =TRUE, destdir =out.dir)
gse = gse$GSE32959_series_matrix.txt.gz

meta.table = pData(gse)
## >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
##
## >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
meta.table = dplyr::mutate(
  meta.table,
  Short_Name = dplyr::case_when(
    grepl("\\s0h", source_name_ch1) ~ "thp",
    grepl("IL-4\\streated\\,\\s12h", source_name_ch1) ~ "th2_12h",
    grepl("IL-4\\streated\\,\\s24h", source_name_ch1) ~ "th2_24h",
    grepl("IL-4\\streated\\,\\s48h", source_name_ch1) ~ "th2_48h",
    grepl("IL-4\\streated\\,\\s72h", source_name_ch1) ~ "th2_72h",
    grepl("IL-12\\streated\\,\\s12h", source_name_ch1) ~ "th1_12h",
    grepl("IL-12\\streated\\,\\s24h", source_name_ch1) ~ "th1_24h",
    grepl("IL-12\\streated\\,\\s48h", source_name_ch1) ~ "th1_48h",
    grepl("IL-12\\streated\\,\\s72h", source_name_ch1) ~ "th1_72h",
    grepl("activated\\,\\s12h", source_name_ch1) ~ "th0_12h",
    grepl("activated\\,\\s24h", source_name_ch1) ~ "th0_24h",
    grepl("activated\\,\\s48h", source_name_ch1) ~ "th0_48h",
    grepl("activated\\,\\s72h", source_name_ch1) ~ "th0_72h",
  )
)
meta.table$tissue = "cord blood"
meta.table$description = gsub("biological\\sreplicate\\s", "", meta.table$description)
meta.table = dplyr::rename(meta.table,
                               Sample_Name = geo_accession,
                               Biological_Replicate = description) %>%
  dplyr::mutate(Short_Name_Rep = paste0(Short_Name, "_", Biological_Replicate))

## >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
## METATABLE
## >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
studies.l = list()

study = meta.table

## add column with # for techn. replicates
c.unique.GSM = plyr::count(study, vars = "Sample_Name")
colnames(c.unique.GSM)[2] = "TECH_REPS"
study =  dplyr::left_join(study, c.unique.GSM, by = 'Sample_Name')
study = distinct(study, Sample_Name, .keep_all = TRUE)

## add column with # for biol. replicates
c.unique.Rep = plyr::count(study, vars = "Short_Name")
colnames(c.unique.Rep)[2] = "BIO_REPS"
study =  dplyr::left_join(study, c.unique.Rep, by = 'Short_Name')

if(!"Study_Name" %in% colnames(study)) study$Study_Name = NA

study$Instrument = "Affymetrix GeneChip Human Genome U133 Plus 2.0"
study = dplyr::select(
  study,
  Sample_Name, Instrument,
  organism_ch1, tissue, Study_Name, Short_Name_Rep, Biological_Replicate,
  Short_Name, TECH_REPS, BIO_REPS
)

study = dplyr::mutate(
  study,
  tissue = case_when(
    grepl("^Peripheral", tissue) ~ "PBMC",
    grepl("^cord", tissue) ~ "Cord_Blood"
  )
)

study$HOURS = gsub("th[1|2|0]_", "", study$Short_Name)
study$HOURS = gsub("thp", "0h", study$HOURS)
study$GSE = "GSE32959"
study$STUDY = "aijö"
study$ORGANIZATION = "Turku_University"

colnames(study) = toupper(colnames(study))

study$HOURS = factor(study$HOURS, levels = c("0h", "12h", "24h", "48h", "72h"))

study = droplevels(study)
study$SAMPLE_ID = study$SHORT_NAME_REP
study = study %>% mutate_if(is.character,as.factor)
rownames(study) = study$SAMPLE_ID

studies.l[["aijö"]] = study

saveRDS(studies.l, file = out.file)


