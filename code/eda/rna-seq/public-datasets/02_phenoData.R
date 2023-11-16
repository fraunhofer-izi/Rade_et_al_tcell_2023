## >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
## LIBRARIES
## >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
.cran_packages = c("dplyr", "yaml", "naturalsort")
.bioc_packages = c("GEOquery")

## library(ggplotify)

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

## >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
##
## >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
manifest = yaml.load_file("assets/manifest.yaml")
home = manifest$homes
work = manifest$workdata

## Tuomela-26967054-GSE52260
meta.table.tuomela = read.table( "code/uap/ag-indep/Tuomela-26967054-GSE52260-thp-th0-th17/experiment/human/SraRunTable", header = T, sep = "\t")
# head(meta.table.tuomela)

## Ullah-29466736-GSE90569
meta.table.ullah = read.table("code/uap/ag-indep/Ullah-29466736-GSE90569-thp-th0-iTreg/experiment/SraRunTable", header = T, sep = "\t")
# head(meta.table.ullah)

## Schmidt-29730990-GSE94396
meta.table.schmidt = read.table("code/uap/ag-indep/Schmidt-29730990-thp-th0-iTreg/GSE94396/experiment/SraRunTable", header = T, sep = "\t")
gse = getGEO("GSE94396")
pd.schmidt = gse$GSE94396_series_matrix.txt.gz
pd.schmidt = pData(pd.schmidt)
pd.schmidt = subset(pd.schmidt, `sample.group:ch1` == "G01" | `sample.group:ch1` == "G02" | `sample.group:ch1` == "G04")
# head(meta.table.schmidt)

## Schmidt-29730990-GSE96538
meta.table.schmidt.ind = read.table("code/uap/ag-indep/Schmidt-29730990-thp-th0-iTreg/GSE96538/experiment/SraRunTable", header = T, sep = "\t")
# head(meta.table.schmidt.ind)

## >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
## MAKE UNIFORM COLLUMNS OF INTEREST
## >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
## Tuomela-26967054-GSE52260
meta.table.tuomela = dplyr::mutate(
  meta.table.tuomela,
  Short_Name = case_when(
    grepl("\\swithout", source_name) ~ "thp",
    grepl("\\s0.5h after activation$", source_name) ~ "th0_0.5h",
    grepl("\\s1h after activation$", source_name) ~ "th0_1h",
    grepl("\\s2h after activation$", source_name) ~ "th0_2h",
    grepl("\\s4h after activation$", source_name) ~ "th0_4h",
    grepl("\\s6h after activation$", source_name) ~ "th0_6h",
    grepl("\\s12h after activation$", source_name) ~ "th0_12h",
    grepl("\\s24h after activation$", source_name) ~ "th0_24h",
    grepl("\\s48h after activation$", source_name) ~ "th0_48h",
    grepl("\\s72h after activation$", source_name) ~ "th0_72h",
    grepl("\\s0.5h after activation and polarization", source_name) ~ "th17_0.5h",
    grepl("\\s1h after activation and polarization", source_name) ~ "th17_1h",
    grepl("\\s2h after activation and polarization", source_name) ~ "th17_2h",
    grepl("\\s4h after activation and polarization", source_name) ~ "th17_4h",
    grepl("\\s6h after activation and polarization", source_name) ~ "th17_6h",
    grepl("\\s12h after activation and polarization", source_name) ~ "th17_12h",
    grepl("\\s24h after activation and polarization", source_name) ~ "th17_24h",
    grepl("\\s48h after activation and polarization", source_name) ~ "th17_48h",
    grepl("\\s72h after activation and polarization", source_name) ~ "th17_72h"
  )
)
meta.table.tuomela = dplyr::rename(meta.table.tuomela,
  Biological_Replicate = human_biological_replicate,
  tissue = source) %>%
dplyr::mutate(
  Short_Name_Rep = paste0(Short_Name,"_",Biological_Replicate))

## Ullah-29466736-GSE90569
meta.table.ullah = dplyr::mutate(
  meta.table.ullah,
  Short_Name = case_when(
    grepl("\\swithout", source_name) ~ "thp",
    grepl("\\s0.5h after activation$", source_name) ~ "th0_0.5h",
    grepl("\\s1h after activation$", source_name) ~ "th0_1h",
    grepl("\\s2h after activation$", source_name) ~ "th0_2h",
    grepl("\\s4h after activation$", source_name) ~ "th0_4h",
    grepl("\\s6h after activation$", source_name) ~ "th0_6h",
    grepl("\\s12h after activation$", source_name) ~ "th0_12h",
    grepl("\\s24h after activation$", source_name) ~ "th0_24h",
    grepl("\\s48h after activation$", source_name) ~ "th0_48h",
    grepl("\\s72h after activation$", source_name) ~ "th0_72h",
    grepl("\\s0.5h after activation and polarization", source_name) ~ "itreg_0.5h",
    grepl("\\s1h after activation and polarization", source_name) ~ "itreg_1h",
    grepl("\\s2h after activation and polarization", source_name) ~ "itreg_2h",
    grepl("\\s4h after activation and polarization", source_name) ~ "itreg_4h",
    grepl("\\s6h after activation and polarization", source_name) ~ "itreg_6h",
    grepl("\\s12h after activation and polarization", source_name) ~ "itreg_12h",
    grepl("\\s24h after activation and polarization", source_name) ~ "itreg_24h",
    grepl("\\s48h after activation and polarization", source_name) ~ "itreg_48h",
    grepl("\\s72h after activation and polarization", source_name) ~ "itreg_72h"
  )
)

meta.table.ullah = dplyr::mutate(
  meta.table.ullah,
  Biological_Replicate = case_when(
    grepl("human biological replicate: 1", biological_replicate) ~ "1",
    grepl("human biological replicate: 2", biological_replicate) ~ "2",
    grepl("human biological replicate: 3", biological_replicate) ~ "3"
  )
) %>%
  dplyr::mutate(
    Short_Name_Rep = paste0(Short_Name,"_",Biological_Replicate))

## Schmidt-29730990-GSE94396
meta.table.schmidt$tissue = "Peripheral blood"
pd.schmidt$source_name = paste0(pd.schmidt$`sample.group:ch1`,"_", pd.schmidt$`sample.time:ch1`)
meta.table.schmidt$source_name =pd.schmidt$source_name[match(meta.table.schmidt$Sample_Name, rownames(pd.schmidt))]
meta.table.schmidt$Biological_Replicate = pd.schmidt$`experiment.donor:ch1`[match(meta.table.schmidt$Sample_Name, rownames(pd.schmidt))]

meta.table.schmidt = dplyr::mutate(
  meta.table.schmidt,
  Short_Name = case_when(
    grepl("G01_0", source_name) ~ "thp",
    grepl("G02_2$", source_name) ~ "th0_2h",
    grepl("G02_6", source_name) ~ "th0_6h",
    grepl("G02_24", source_name) ~ "th0_24h",
    grepl("G02_48", source_name) ~ "th0_48h",
    grepl("G02_144", source_name) ~ "itreg_144h",
    grepl("G04_2$", source_name) ~ "itreg_2h",
    grepl("G04_6", source_name) ~ "itreg_6h",
    grepl("G04_24", source_name) ~ "itreg_24h",
    grepl("G04_48", source_name) ~ "itreg_48h",
    grepl("G04_144", source_name) ~ "itreg_144h"
  )
)

meta.table.schmidt = dplyr::mutate(
  meta.table.schmidt,
  Biological_Replicate = case_when(
    grepl("A", Biological_Replicate) ~ "1",
    grepl("B", Biological_Replicate) ~ "2",
    grepl("D", Biological_Replicate) ~ "3"
  )
) %>%
  dplyr::mutate(
    Short_Name_Rep = paste0(Short_Name,"_",Biological_Replicate)
    )

## Schmidt-29730990-GSE96538
meta.table.schmidt.ind = dplyr::mutate(
  meta.table.schmidt.ind,
  Short_Name = case_when(
    grepl("\\swithout", source_name) ~ "thp",
    grepl("\\s0.5h after activation$", source_name) ~ "th0_0.5h",
    grepl("\\s1h after activation$", source_name) ~ "th0_1h",
    grepl("\\s2h after activation$", source_name) ~ "th0_2h",
    grepl("\\s4h after activation$", source_name) ~ "th0_4h",
    grepl("\\s6h after activation$", source_name) ~ "th0_6h",
    grepl("\\s12h after activation$", source_name) ~ "th0_12h",
    grepl("\\s24h after activation$", source_name) ~ "th0_24h",
    grepl("\\s48h after activation$", source_name) ~ "th0_48h",
    grepl("\\s72h after activation$", source_name) ~ "th0_72h",
    grepl("\\s96h after activation$", source_name) ~ "th0_96h",
    grepl("\\s120h after activation$", source_name) ~ "th0_120h",
    grepl("\\s0.5h after activation and polarization", source_name) ~ "itreg_0.5h",
    grepl("\\s1h after activation and polarization", source_name) ~ "itreg_1h",
    grepl("\\s2h after activation and polarization", source_name) ~ "itreg_2h",
    grepl("\\s4h after activation and polarization", source_name) ~ "itreg_4h",
    grepl("\\s6h after activation and polarization", source_name) ~ "itreg_6h",
    grepl("\\s12h after activation and polarization", source_name) ~ "itreg_12h",
    grepl("\\s24h after activation and polarization", source_name) ~ "itreg_24h",
    grepl("\\s48h after activation and polarization", source_name) ~ "itreg_48h",
    grepl("\\s72h after activation and polarization", source_name) ~ "itreg_72h",
    grepl("\\s96h after activation and polarization", source_name) ~ "itreg_96h",
    grepl("\\s120h after activation and polarization", source_name) ~ "itreg_120h"
  )
)

meta.table.schmidt.ind = dplyr::mutate(
  meta.table.schmidt.ind,
  Biological_Replicate = case_when(
    grepl("human biological replicate: 1", biological_replicate) ~ "1",
    grepl("human biological replicate: 2", biological_replicate) ~ "2",
    grepl("human biological replicate: 3", biological_replicate) ~ "3"
  )
) %>%
  dplyr::mutate(
    Short_Name_Rep = paste0(Short_Name,"_",Biological_Replicate))

## >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
## ONE METATABLE
## >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
studies = c("meta.table.tuomela", "meta.table.ullah", "meta.table.schmidt", "meta.table.schmidt.ind")
gse.nbr = c("GSE52260", "GSE90569", "GSE94396", "GSE96538")
organization = c("Aalto_University", "Aalto_University", "Karolinska_Institute", "Aalto_University")

studies.l = list()
for (i in 1:length(studies)) {

  study = eval(parse(text = studies[i]))

  ## add column with # for techn. replicates
  c.unique.GSM = plyr::count(study, vars = "Sample_Name")
  colnames(c.unique.GSM)[2] = "TECH_REPS"
  study =  dplyr::left_join(study, c.unique.GSM, by = 'Sample_Name')
  study = distinct(study, Sample_Name, .keep_all = TRUE)

  ## add column with # for biol. replicates
  c.unique.Rep = plyr::count(study, vars = "Short_Name")
  colnames(c.unique.Rep)[2] = "BIO_REPS"
  study =  dplyr::left_join(study, c.unique.Rep, by = 'Short_Name')

  study = dplyr::select(
    study,
    Sample_Name, AvgSpotLen, Instrument, LibraryLayout,
    Organism, tissue, Short_Name_Rep, Biological_Replicate,
    Short_Name, TECH_REPS, BIO_REPS
  )

  study = dplyr::mutate(
    study,
    tissue = case_when(
      grepl("^Peripheral", tissue) ~ "PBMC",
      grepl("^Cord", tissue) ~ "Cord_Blood"
    )
  )

  study$HOURS = gsub("th17_", "", study$Short_Name)
  study$HOURS = gsub("itreg_", "", study$HOURS)
  study$HOURS = gsub("th0_", "", study$HOURS)
  study$HOURS = gsub("thp", "0h", study$HOURS)
  study$GSE = gse.nbr[i]
  study$STUDY = gsub("meta.table.", "", studies[i])
  study$ORGANIZATION = organization[i]
  study$CONDITION = gsub("_.*", "", study$Short_Name)

  colnames(study) = toupper(colnames(study))
  study = droplevels(study)
  study = study %>% mutate_if(is.character,as.factor)
  rownames(study) = study$SAMPLE_NAME


  study$BIOLOGICAL_REPLICATE = paste0("R", study$BIOLOGICAL_REPLICATE, "_", study$STUDY)
  lvls = naturalsort(unique(study$BIOLOGICAL_REPLICATE))
  study$BIOLOGICAL_REPLICATE = factor(study$BIOLOGICAL_REPLICATE, levels = lvls)

  studies.l[[ gsub("meta.table.", "", studies[i]) ]] = study
  # print(study)
  # cat("------------------------------------\n")
}

metaData.rna.seq.kintec.l = studies.l
saveRDS(metaData.rna.seq.kintec.l, file = paste0(work, manifest$phenodata_rnaseq))

