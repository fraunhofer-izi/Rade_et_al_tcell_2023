library(ShortRead, quietly = T)
library(ggplot2, quietly = T)
# .cran_packages <- c("ggplot2", "gridExtra")
# .bioc_packages <- c("dada2")
# .inst <- .cran_packages %in% installed.packages()
# if(any(!.inst)) {
#   install.packages(.cran_packages[!.inst])
# }
# .inst <- .bioc_packages %in% installed.packages()
# if(any(!.inst)) {
#   source("http://bioconductor.org/biocLite.R")
#   biocLite(.bioc_packages[!.inst], ask = F)
# }
# # Load packages into session, and print package version
# sapply(c(.cran_packages, .bioc_packages), require, character.only = TRUE)

#' Plot quality profile of a fastq file.
#' 
#' This function plots a visual summary of the distribution of quality scores
#' as a function of sequence position for the input fastq file(s).
#' 
#' The distribution of quality scores at each position is shown as a grey-scale
#' heat map, with dark colors corresponding to higher frequency. The plotted lines
#' show positional summary statistics: green is the mean, orange is the median, and
#' the dashed orange lines are the 25th and 75th quantiles.
#' 
#' If the sequences vary in length, a red line will be plotted showing the percentage
#' of reads that extend to at least that position.
#' 
#' @param fl (Required). \code{character}.
#'  File path(s) to fastq or fastq.gz file(s).
#' 
#' @param n (Optional). Default 500,000.
#'  The number of records to sample from the fastq file.
#' 
#' @param aggregate (Optional). Default FALSE.
#'  If TRUE, compute an aggregate quality profile for all fastq files provided.
#'
#' @return A \code{\link{ggplot}2} object.
#'  Will be rendered to default device if \code{\link{print}ed},
#'  or can be stored and further modified.
#'  See \code{\link{ggsave}} for additional options.
#'  
#' @importFrom ShortRead qa
#' @import ggplot2
#' 
#' @export
#' 
#' @examples
#' plotQualityProfile(system.file("extdata", "sam1F.fastq.gz", package="dada2"))
#' 

plotQualityProfile <- function(fl, n=100000, aggregate=FALSE) {
  # fl = paste0("~/ribo-root", sort(manifest$absolute.filepath)[1:1])
  # fl = gsub("homes/olymp/michael.rade/","",fl)
  
  statdf <- data.frame(Cycle=integer(0), Mean=numeric(0), Q25=numeric(0), Q50=numeric(0), Q75=numeric(0), Cum=numeric(0), file=character(0))
  anndf <- data.frame(minScore=numeric(0), label=character(0), rclabel=character(0), rc=numeric(0), file=character(0))
  
  FIRST <- TRUE
  for(f in fl[!is.na(fl)]) {
    srqa <- qa(f, n=n)
    df <- srqa[["perCycle"]]$quality
    rc <- sum(srqa[["readCounts"]]$read) # Handle aggregate form from qa of a directory
    if (rc >= n) { 
      rclabel <- paste("Reads >= ", n)
    } else {
      rclabel <- paste("Reads: ", rc)
    }
    # Calculate summary statistics at each position
    means <- rowsum(df$Score*df$Count, df$Cycle)/rowsum(df$Count, df$Cycle)
    get_quant <- function(xx, yy, q) { xx[which(cumsum(yy)/sum(yy) >=q)][[1]] }
    q25s <- by(df, df$Cycle, function(foo) get_quant(foo$Score, foo$Count, 0.25), simplify=TRUE)
    q50s <- by(df, df$Cycle, function(foo) get_quant(foo$Score, foo$Count, 0.5), simplify=TRUE)
    q75s <- by(df, df$Cycle, function(foo) get_quant(foo$Score, foo$Count, 0.75), simplify=TRUE)
    cums <- by(df, df$Cycle, function(foo) sum(foo$Count), simplify=TRUE)
    if(!all(sapply(list(names(q25s), names(q50s), names(q75s), names(cums)), identical, rownames(means)))) {
      stop("Calculated quantiles/means weren't compatible.")
    }
    if(FIRST) {
      plotdf <- cbind(df, file=basename(f))
      FIRST <- FALSE
    } else { plotdf <- rbind(plotdf, cbind(df, file=basename(f))) }
    statdf <- rbind(statdf, data.frame(Cycle=as.integer(rownames(means)), Mean=means, 
                                       Q25=as.vector(q25s), Q50=as.vector(q50s), Q75=as.vector(q75s), Cum=10*as.vector(cums)/rc, file=basename(f)))
    anndf <- rbind(anndf, data.frame(minScore=min(df$Score), label=basename(f), rclabel=rclabel, rc=rc, file=basename(f)))
  }
  anndf$minScore <- min(anndf$minScore)
  # Create plot
  if (aggregate) {
    plotdf.summary <- aggregate(Count ~ Cycle + Score, plotdf, sum)
    plotdf.summary$label <- paste(nrow(anndf), "files (aggregated)")
    means <- rowsum(plotdf.summary$Score*plotdf.summary$Count, plotdf.summary$Cycle)/rowsum(plotdf.summary$Count, plotdf.summary$Cycle)
    q25s <- by(plotdf.summary, plotdf.summary$Cycle, function(foo) get_quant(foo$Score, foo$Count, 0.25), simplify=TRUE)
    q50s <- by(plotdf.summary, plotdf.summary$Cycle, function(foo) get_quant(foo$Score, foo$Count, 0.5), simplify=TRUE)
    q75s <- by(plotdf.summary, plotdf.summary$Cycle, function(foo) get_quant(foo$Score, foo$Count, 0.75), simplify=TRUE)
    cums <- by(plotdf.summary, plotdf.summary$Cycle, function(foo) sum(foo$Count), simplify=TRUE)
    statdf.summary <- data.frame(Cycle=as.integer(rownames(means)), Mean=means, Q25=as.vector(q25s), Q50=as.vector(q50s), Q75=as.vector(q75s), Cum=10*as.vector(cums)/rc)
    p <- ggplot(data=plotdf.summary, aes(x=Cycle, y=Score)) + geom_tile(aes(fill=Count)) + 
      scale_fill_gradient(low="#F5F5F5", high="black") + 
      geom_line(data=statdf.summary, aes(y=Mean), color="#66C2A5") +
      geom_line(data=statdf.summary, aes(y=Q25), color="#FC8D62", size=0.25, linetype="dashed") +
      geom_line(data=statdf.summary, aes(y=Q50), color="#FC8D62", size=0.25) +
      geom_line(data=statdf.summary, aes(y=Q75), color="#FC8D62", size=0.25, linetype="dashed") +
      ylab("Quality Score") + xlab("Cycle") +
      annotate("text", x=0, y=2, label=sprintf("Total reads: %d", sum(anndf$rc)), color="black", hjust=0) + 
      theme_bw() + theme(panel.grid=element_blank(), strip.text=element_text(size=8)) + guides(fill=FALSE) + 
      facet_wrap(~label, ncol=6, nrow=6) + ylim(c(0,NA))
  } else {
    p <- ggplot(data=plotdf, aes(x=Cycle, y=Score)) + geom_tile(aes(fill=Count)) + 
      scale_fill_gradient(low="#F5F5F5", high="black") + 
      geom_line(data=statdf, aes(y=Mean), color="#66C2A5") +
      geom_line(data=statdf, aes(y=Q25), color="#FC8D62", size=0.25, linetype="dashed") +
      geom_line(data=statdf, aes(y=Q50), color="#FC8D62", size=0.25) +
      geom_line(data=statdf, aes(y=Q75), color="#FC8D62", size=0.25, linetype="dashed") +
      ylab("Quality Score") + xlab("Cycle") +
      theme_bw() + theme(panel.grid=element_blank(), strip.text=element_text(size=8)) + guides(fill=FALSE) +
      geom_text(data=anndf, aes(x=0, label=rclabel, y=2), color="black", hjust=0) + 
      facet_wrap(~file, ncol=6, nrow=6) + ylim(c(0,NA))
  }
  p
}


manifest = read.table("~/projects/2018-MAVO/studies/kinetic/2018-Schmidt-29730990/GSE94396/scripts/rawdata-manifest",
                      header = F, stringsAsFactors = F)
manifest = manifest$V1
length(manifest)

png("demux_qcplot_1.png", width = 800, height = 800)
plotQualityProfile(sort(manifest[1:36]), n=100000)
dev.off()

# manifest = read.table("~/projects/2018-Mikrobiom-Rostock/scripts/06_maxEEfilter_manifest",
#                       sep = ",", header = T, stringsAsFactors = F)
# png("qcplot-final.png", width = 800, height = 800)
# plotQualityProfile(sort(manifest$absolute.filepath)[1:36], n=650000)
# dev.off()




