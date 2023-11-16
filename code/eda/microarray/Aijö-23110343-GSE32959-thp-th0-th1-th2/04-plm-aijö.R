.cran_packages = c("yaml")
.bioc_packages = c("Biobase", "affy", "affyPLM")
.url_packages = c("hgu133plus2hsgencodegcdf", "hgu133plus2hsgencodegprobe",
                  "pd.hgu133plus2.hs.gencodeg")
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

list.of.packages = c(.cran_packages, .bioc_packages, .url_packages)

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

# obj with rawdata
affybatch = readRDS(paste0(workdata, manifest$affybatch_aijö_thp_th1))

# Output paths for array/pseudo images
array.images = manifest$array_images_aijö_thp_th1

plm.obj.out = paste0(workdata, manifest$plm_obj_aijö_thp_th1)
plm.nuse.rle.out = paste0(workdata, manifest$plm_nuse_rle_aijö_thp_th1)

## >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
## PLMset
## >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
plm.set = affyPLM::fitPLM(affybatch)
nuse = affyPLM::NUSE(plm.set, type = "values")
rle = affyPLM::RLE(plm.set, type = "values")

saveRDS(plm.set, file = plm.obj.out)
plm.nuse.rle = list(nuse = nuse, rle = rle)
saveRDS(plm.nuse.rle, file = plm.nuse.rle.out)

## >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
## Array/pseudo images
## >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# plm.set = readRDS(plm.obj.out)

for (i in 1:nrow(pData(affybatch))){
  name = paste0(array.images, "/array/arrayimage_", i, ".png")
  png(name, width = 350, height = 350, units = "px", pointsize = 18, res = NA)
  par(mar=c(2.0, 1.1, 1.6, 1.1), oma=c(1, 1, 0, 0))
  image(affybatch[,i], main=pData(affybatch)$SAMPLE_ID[i], cex = 2, col=pseudoPalette(low="yellow",high="red"))
  dev.off()
}
# Weights
for (i in 1:nrow(pData(plm.set))){
  name = paste0(array.images, "/pseudo-weights/pseudoimage_weights_", i, ".png")
  png(name, width = 350, height = 350, units = "px", pointsize = 18, res = NA)
  # png(name, height=4,width=4,res=200,units="in")
  par(mar=c(2.0, 1.1, 1.6, 1.1), oma=c(1, 1, 0, 0))
  image(plm.set, which=i, main=pData(plm.set)$SAMPLE_ID[i])
  dev.off()
}
# Residuals
for (i in 1:nrow(pData(plm.set))){
  name = paste0(array.images, "/pseudo-residuals/pseudoimage_residuals_", i, ".png")
  png(name, width = 350, height = 350, units = "px", pointsize = 18, res = NA)
  # png(name, height=4,width=4,res=200,units="in")
  par(mar=c(2.0, 1.1, 1.6, 1.1), oma=c(1, 1, 0, 0))
  image(plm.set, which=i, type='resids', main=pData(plm.set)$SAMPLE_ID[i])
  dev.off()
}
