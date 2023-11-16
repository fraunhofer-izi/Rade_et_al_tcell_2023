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

output = paste0(workdata, manifest$phenodata_elo_thp_th0_th2, "/elo-thp-th0-th2.Rds")

## >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
##
## >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
gse = getGEO("GSE17974", GSEMatrix =TRUE, destdir = paste0(workdata, manifest$phenodata_elo_thp_th0_th2))
gse = gse$GSE17974_series_matrix.txt.gz

## Elo-20620947-GSE17974
meta.table.elo = pData(gse)

## >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
##
## >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
meta.table.elo = dplyr::mutate(
  meta.table.elo,
  Short_Name = dplyr::case_when(
    grepl("\\s0h", source_name_ch1) ~ "thp",
    grepl("treated\\,\\s0.5h", source_name_ch1) ~ "th2_0.5h",
    grepl("treated\\,\\s1h", source_name_ch1) ~ "th2_1h",
    grepl("treated\\,\\s2h", source_name_ch1) ~ "th2_2h",
    grepl("treated\\,\\s4h", source_name_ch1) ~ "th2_4h",
    grepl("treated\\,\\s6h", source_name_ch1) ~ "th2_6h",
    grepl("treated\\,\\s12h", source_name_ch1) ~ "th2_12h",
    grepl("treated\\,\\s24h", source_name_ch1) ~ "th2_24h",
    grepl("treated\\,\\s48h", source_name_ch1) ~ "th2_48h",
    grepl("treated\\,\\s72h", source_name_ch1) ~ "th2_72h",
    grepl("activated\\,\\s0.5h", source_name_ch1) ~ "th0_0.5h",
    grepl("activated\\,\\s1h", source_name_ch1) ~ "th0_1h",
    grepl("activated\\,\\s2h", source_name_ch1) ~ "th0_2h",
    grepl("activated\\,\\s4h", source_name_ch1) ~ "th0_4h",
    grepl("activated\\,\\s6h", source_name_ch1) ~ "th0_6h",
    grepl("activated\\,\\s12h", source_name_ch1) ~ "th0_12h",
    grepl("activated\\,\\s24h", source_name_ch1) ~ "th0_24h",
    grepl("activated\\,\\s48h", source_name_ch1) ~ "th0_48h",
    grepl("activated\\,\\s72h", source_name_ch1) ~ "th0_72h",
  )
)
meta.table.elo$description = gsub("biological\\sreplicate\\s", "", meta.table.elo$description)
meta.table.elo = dplyr::rename(meta.table.elo,
                               Sample_Name = geo_accession,
                               Biological_Replicate = description,
                               tissue = `tissue:ch1`) %>%
  dplyr::mutate(Short_Name_Rep = paste0(Short_Name, "_", Biological_Replicate))

## >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
## METATABLE
## >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
studies.l = list()

study = meta.table.elo

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

study$HOURS = gsub("th[17|2|0]_", "", study$Short_Name)
study$HOURS = gsub("thp", "0h", study$HOURS)
study$GSE = "GSE17974"
study$STUDY = "elo"
study$ORGANIZATION = "Turku_University"

colnames(study) = toupper(colnames(study))

study$HOURS = factor(study$HOURS, levels = c("0h", "0.5h", "1h", "2h", "4h",
                                             "6h", "12h", "24h", "48h", "72h"))

study = droplevels(study)
study$SAMPLE_ID = study$SHORT_NAME_REP
study = study %>% mutate_if(is.character,as.factor)
rownames(study) = study$SAMPLE_ID

studies.l[["elo"]] = study

saveRDS(studies.l, file = output)


