library(edgeR)

## >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
## Count data from HiSeqfc3f43 (Michael)
## >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
count_dirs<-"/mnt/workdata1/user_worktmp/michael.rade/sphere/work/uap_out/2018-MAVO/2018-Bonnal-E-MTAB-2319/htseq_hisat2/"

design.table = NULL
design.table$fileName = list.files(count_dirs, pattern = "_counts.txt$", recursive = TRUE)
length(design.table$fileName)

setwd(count_dirs)
x <- readDGE(design.table$fileName, columns=c(1,2), header=FALSE)
x$counts = head(x$counts, -5)

raw.counts = as.data.frame(x$counts)
colnames(raw.counts) = gsub("-htseq_counts","", basename(colnames(raw.counts)) )

dim(raw.counts)
head(raw.counts)
tail(raw.counts)

## >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
## 
## >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
tmp = data.frame(EnsemblID = rownames(raw.counts))
raw.counts.fin = cbind(tmp, raw.counts)

dim(raw.counts.fin)

write.csv(raw.counts.fin, "/mnt/workdata1/user_worktmp/michael.rade/sphere/work/2018-MAVO/analysis_tables/E-MTAB-2319/counts/htseq/counts.csv", quote = FALSE, row.names = FALSE)
save(raw.counts.fin, file="/mnt/workdata1/user_worktmp/michael.rade/sphere/work/2018-MAVO/analysis_tables/E-MTAB-2319/counts/htseq/counts.Rdata")



