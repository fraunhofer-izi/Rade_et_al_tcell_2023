# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# NORMALIZATION METHODS
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
z_score = function(z) {
  rowmean = apply(z, 1, mean, na.rm=TRUE)
  rowsd = apply(z, 1, sd, na.rm=TRUE)
  rv = sweep(z, 1, rowmean,"-")
  rv = sweep(rv, 1, rowsd, "/")
  return(rv)
}
centering = function(z) {
  rowmean = apply(z, 1, mean)
  rv = sweep(z, 1, rowmean,"-")
  return(rv)
}

unit_var = function(z) {
  rowsd = apply(z, 1, sd)
  rv = sweep(z, 1, rowsd,"/")
  return(rv)
}

minmax = function(x) {
  return((x- min(x)) /(max(x)-min(x)))
}

minmax_colwise <- function(x) {
  x <- as.matrix(x)
  minAttr=apply(x, 2, min)
  maxAttr=apply(x, 2, max)
  x <- sweep(x, 2, minAttr, FUN="-")
  x=sweep(x, 2,  maxAttr-minAttr, "/")
  attr(x, 'normalized:min') = minAttr
  attr(x, 'normalized:max') = maxAttr
  return (x)
}

rank_norm = function(mat) {
  return(apply(mat, 2, function(y) 10000*(rank(y) -1)/ length(y)))
}

colnorm = function(z) {
  colsum = apply(z, 2, sum)
  coln = sweep(z, 2, colsum, "/")
  return(coln)
}

quantile_norm = function(data)
{
  quantiles = c(rep(0, nrow(data)))
  na.cols = c(rep(0, nrow(data)))

  for (i in 1:ncol(data))
  {
    data.col = data [order(data [,i], decreasing=TRUE) ,i]
    not.na = which(!is.na(data.col))

    quantiles[not.na] = quantiles[not.na] + data.col[not.na]
    na.cols[not.na] = na.cols[not.na] + 1
  }
  quantiles = quantiles / na.cols

  for (i in 1:ncol(data))
  {
    data [order(data [,i], decreasing=TRUE) ,i] = quantiles
  }

  return(data)
}

## >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
## PCA
## >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
pca_plot = function(object,
                    intgroup = "condition",
                    pc1.pn = "group1",
                    pc2.pn = "group2",
                    pc1 = 1,
                    pc2 = 2,
                    ntop = 500,
                    type = NULL,
                    log = NULL,
                    point.size = 1.7) {

  if (class(object) == "DESeqTransform") {
    object.exprs = assay(object)
    object.pD = colData(object)
  }
  if (class(object) == "DESeqDataSet") {
    object.exprs = assay(object)
    object.pD = colData(object)
  }
  if (class(object) == "ExpressionSet") {
    object.exprs = exprs(object)
    object.pD = pData(object)
  }
  if (class(object) == "DGEList") {
    object.exprs = object$counts
    object.pD = object$samples
  }
  if (class(object) == "SummarizedExperiment") {
    if (is.null(type)) {
      object.exprs = assay(object)
    } else {
      object.exprs = assays(object)[[type]]
    }
    if (!is.null(log)) {
      object.exprs = log(object.exprs)
    }
    object.pD = colData(object)
  }


  # calculate the variance for each gene
  rv = rowVars(object.exprs)
  # select the ntop genes by variance
  select = order(rv, decreasing = TRUE)[seq_len(min(ntop, length(rv)))]
  # perform a PCA on the data in assay(x) for the selected genes
  pca = prcomp(t(object.exprs[select, ]))
  # the contribution to the total variance for each component
  percentVar = pca$sdev ^ 2 / sum(pca$sdev ^ 2)

  if (!all(intgroup %in% names(object.pD))) {
    stop("the argument 'intgroup' should specify columns of the pheno table")
  }

 intgroup.df = as.data.frame(object.pD[, intgroup, drop = FALSE])

  # assembly the data for the plot
  d = data.frame(
    PC1 = pca$x[, 1],
    PC2 = pca$x[, 2],
    PC3 = pca$x[, 3],
    PC4 = pca$x[, 4],
    PC4 = pca$x[, 5],
    intgroup.df,
    name = colnames(object)
  )

  if (length(intgroup) == 1) {
    pca1 = ggplot(data = d, aes_string(
      x = d[, eval(pc1)],
      y = d[, eval(pc2)],
      color = intgroup[1]
    )) + geom_point(size = point.size)
    pca1 = pca1 + labs(colour = pc1.pn)
  } else {
    pca1 = ggplot(data = d,
                  aes_string(x = d[, eval(pc1)], y = d[, eval(pc2)], color = intgroup[1], shape = intgroup[2])) +
      geom_point(size = point.size)
    pca1 = pca1 + guides(color = guide_legend(order = 1),
                         shape = guide_legend(order = 2))
    pca1 = pca1 + labs(colour = pc1.pn, shape = pc2.pn)
  }
  pca1 = pca1 + xlab(paste0("PC", eval(pc1), ": ", round(percentVar[eval(pc1)] * 100), "% var"))
  pca1 = pca1 + ylab(paste0("PC", eval(pc2), ": ", round(percentVar[eval(pc2)] * 100), "% var"))
  pca1 = pca1

  if (is.numeric(d[, intgroup[1]])) {
    pca1
  } else {
    if (length(levels(d[, intgroup[1]])) > 4) {
      pca1 = pca1 + scale_colour_manual(values = col.blind.10, na.value = "black")
    } else {
      pca1 = pca1 + scale_colour_manual(values = col.blind.4, na.value = "black")
    }
    pca1
  }
}

## >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
## MDS
## >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
mds_plot = function(object,
                    intgroup = "condition",
                    dim1.pn = "group1",
                    dim2.pn = "group2",
                    dim1 = 1,
                    dim2 = 2,
                    ntop = 500) {
  # calculate the variance for each gene
  rv = rowVars(assay(object))
  # select the ntop genes by variance
  select = order(rv, decreasing = TRUE)[seq_len(min(ntop, length(rv)))]

  sampleDists = dist(t(assay(object)[select, ]))
  sampleDistMatrix = as.matrix(sampleDists)
  mdsData = data.frame(cmdscale(sampleDistMatrix, k = 4))
  mds = cbind(mdsData, as.data.frame(colData(object)))

  if (!all(intgroup %in% names(colData(object)))) {
    stop("the argument 'intgroup' should specify columns of colData(dds)")
  }

  intgroup.df = as.data.frame(colData(object)[, intgroup, drop = FALSE])

  # assembly the data for the plot
  d = data.frame(
    DIM1 = mds[, "X1"],
    DIM2 = mds[, "X2"],
    DIM3 = mds[, "X3"],
    DIM4 = mds[, "X4"],
    intgroup.df,
    name = colnames(object)
  )

  if (length(intgroup) == 1) {
    mds1 = ggplot(data = d, aes_string(
      x = d[, eval(dim1)],
      y = d[, eval(dim2)],
      color = intgroup[1]
    )) + geom_point(size = .9)
    mds1 = mds1 + labs(colour = dim1.pn)
  } else {
    mds1 = ggplot(data = d,
                  aes_string(
                    x = d[, eval(dim1)],
                    y = d[, eval(dim2)],
                    color = intgroup[1],
                    shape = intgroup[2]
                  )) +
      geom_point(size = .9)
    mds1 = mds1 + guides(color = guide_legend(order = 1),
                         shape = guide_legend(order = 2))
    mds1 = mds1 + labs(colour = dim1.pn, shape = dim2.pn)
  }
  mds1 = mds1 + xlab(paste0("DIM", eval(dim1)))
  mds1 = mds1 + ylab(paste0("DIM", eval(dim2)))
  mds1 = mds1
  if (is.numeric(d[, intgroup[1]])) {
    mds1
  } else {
    if (length(levels(d[, intgroup[1]])) > 4) {
      mds1 = mds1 + scale_colour_manual(values = col.blind.9, na.value = "black")
    } else {
      mds1 = mds1 + scale_colour_manual(values = col.blind.4, na.value = "black")
    }
    mds1
  }
}

## >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
## FOR HEATMAP COLORS
## >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
degColors <- function(ann,
                      cat_values = c("#88CCEE", "#332288"),
                      four_values = c("#44AA99", "#CC6677", "#DDCC77", "#004488"),
                      five_values = c("#77AADD","#99DDE1","#44BB99","#BBCC33","#AAAA00"),
                      six_values = c("#77AADD","#99DDE1","#44BB99","#BBCC33", "#AAAA00", "#EE8866"),
                      tol21rainbow= c("#771155", "#AA4488", "#CC99BB", "#114477", "#4477AA", "#77AADD",
                                      "#117777", "#44AAAA", "#77CCCC", "#117744", "#44AA77", "#88CCAA",
                                      "#777711", "#AAAA44", "#DDDD77", "#774411", "#AA7744", "#DDAA77",
                                      "#771122", "#AA4455", "#DD7788")) {
  col <- lapply(names(ann), function(a) {

    if (length(unique(ann[[a]])) < 3) {
      v <- cat_values[1:length(unique(ann[[a]]))]
      names(v) <- unique(ann[[a]])
      return(v)
    }
    if (length(unique(ann[[a]])) <= 4) {
      v <- four_values[1:length(unique(ann[[a]]))]
      names(v) <- unique(ann[[a]])
      return(v)
    }
    if (length(unique(ann[[a]])) == 5) {
      v <- five_values[1:length(unique(ann[[a]]))]
      names(v) <- unique(ann[[a]])
      return(v)
    }
    if (length(unique(ann[[a]])) == 6) {
      v <- six_values[1:length(unique(ann[[a]]))]
      names(v) <- unique(ann[[a]])
      return(v)
    }
    if (length(unique(ann[[a]])) > 6) {
      v <- tol21rainbow[1:length(unique(ann[[a]]))]
      names(v) <- unique(ann[[a]])
      return(v)
    }
  })
  names(col) <- names(ann)
  col
}

## >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
## MA plot
## >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
maplot = function (res,
                   thresh = alpha,
                   yinterceptUp = 1,
                   yinterceptDn = -1,
                   method = "deseq2",
                   lfc = 1) {

  if (method == "deseq2") {
    res = res[!is.na(res$pvalue), ]
    results = data.frame(res) %>% dplyr::mutate(sig = case_when(padj < thresh ~ "yes", TRUE ~ "no"),
                                              sig = factor(sig, levels = c("yes", "no")))
  }
  if (method == "limma") {
    res = res[!is.na(res$P.Value), ]
    results = data.frame(res) %>% dplyr::mutate(sig = case_when(adj.P.Val < thresh & abs(logFC) > lfc ~ "yes", TRUE ~ "no"),
                                                sig = factor(sig, levels = c("yes", "no")))
    results = results[, c("logFC", "AveExpr", "sig")]
    colnames(results) = c("log2FoldChange", "baseMean", "sig")
  }

  results %>% ggplot(aes(baseMean, log2FoldChange, col = sig)) +
    geom_point(size = rel(0.4), alpha = 0.3) +
    geom_point(
      data = results[results$sig == "yes", ],
      aes(baseMean, log2FoldChange),
      alpha = 0.3,
      size = rel(0.4)
    ) +
    xscale("log2", .format = TRUE) +
    xlab("Mean of normalized counts") + labs(colour = paste0("FDR < ", thresh))  +
    theme(aspect.ratio = 1) +
    scale_color_manual(values = c("#BB5566", "grey")) +
    geom_hline(yintercept = yinterceptUp,
               col = "#4477AA",
               lwd = .3) +
    geom_hline(yintercept = yinterceptDn,
               col = "#4477AA",
               lwd = .3) +
    geom_hline(yintercept = 0,
               col = "black",
               lwd = .3)
}

## >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
## VOLCANO plot
## >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# res = dge.results$Cig_3_vs_LPS ; thresh = alpha; yinterceptUp = 1; yinterceptDn = -1;  method = "limma"
volcano <- function (res,
                     thresh = alpha,
                     yinterceptUp = 1,
                     yinterceptDn = -1,
                     method = "deseq2") {

  if (method == "limma") {
    colnames(res)[colnames(res) == "logFC"] = "log2FoldChange"
    colnames(res)[colnames(res) == "adj.P.Val"] = "padj"
  }

  res = res[ !is.na(res$padj), ]
  results = data.frame(res) %>% dplyr::mutate(
    sig = case_when(padj < thresh ~ "yes", TRUE ~ "no"),
    sig = factor(sig, levels = c("yes", "no"))
  )

  x.axis.max = max(max(res$log2FoldChange), abs(min(res$log2FoldChange)))
  results %>% ggplot(aes(log2FoldChange, -log10(padj), col=sig)) +
    geom_point(size = rel(0.4), alpha = 0.3) +
    geom_point(data = results[results$sig=="yes", ], aes(log2FoldChange, -log10(padj)), size = rel(0.4), alpha = 0.3) +
    xlab("log2(FoldChange)") + labs(colour=paste0("FDR < ",thresh))  +
    theme(aspect.ratio = 1) + xlim(-abs(x.axis.max), abs(x.axis.max)) +
    scale_color_manual(values = c("#BB5566", "grey")) +
    geom_vline(xintercept = yinterceptUp, col = "#4477AA", lwd = .3, linetype="dashed") +
    geom_vline(xintercept = yinterceptDn, col = "#4477AA", lwd = .3, linetype="dashed") +
    geom_hline(yintercept = -log10(thresh), col = "black", lwd = .3, linetype="dashed")
}


## >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
## GET UNTION SET OF INTEREST
## >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
fromList = function (input) {
  elements = unique(unlist(input))
  data = unlist(lapply(input, function(x) {
    x <- as.vector(match(elements, x))
  }))
  data[is.na(data)] = as.integer(0)
  data[data != 0] = as.integer(1)
  data = data.frame(matrix(data, ncol = length(input), byrow = F))
  data = data[which(rowSums(data) != 0), ]
  names(data) = names(input)
  row.names(data) = elements
  return(data)
}

## >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
## COUNT PLOT
## >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
countplot = function(object, xs, n=9, genes=NULL,
                          pairs=NULL,
                          slot=1L,
                          log = FALSE,
                          de.log = FALSE,
                          z_score = FALSE,
                          type = NULL,
                          xsLab=xs,
                          color = "black",
                          ncol = 5,
                          intgroup = NULL,
                          point.size = .8,
                          font.size = 8,
                          col.p = NULL,
                          plattform = "SEQ"){

  if (is.null(col.p)) {
    col.p = c("#44AA99", "#CC6677", "#DDCC77", "#4477AA",
              "#AA4499", "#88CCEE")
  }

  if (class(object) == "DESeqTransform") {
    object.exprs = assay(object)
    object.pD = colData(object)
  }
  if (class(object) == "DESeqDataSet") {
    object.exprs = assay(object)
    object.pD = colData(object)
  }
  if (class(object) == "ExpressionSet") {
    object.exprs = exprs(object)
    object.pD = pData(object)
  }
  if (class(object) == "DGEList") {
    object.exprs = object$counts
    object.pD = object$samples
  }
  if (class(object) == "SummarizedExperiment") {
    if (is.null(type)) {
      object.exprs = assay(object)
    } else {
      object.exprs = assays(object)[[type]]
    }
    object.pD = colData(object)
  }

  metadata = data.frame(object.pD)
  if (log == TRUE) {
    if (plattform == "SEQ") {object.exprs = log2(object.exprs)}
    if (plattform == "ARRAY") {object.exprs = log2(object.exprs)}
  }
  if (de.log == TRUE) {
    if (plattform == "SEQ") {object.exprs = 2^(object.exprs)}
    if (plattform == "ARRAY") {object.exprs = 2^(object.exprs)}
  }
  if (z_score == TRUE) {
    object.exprs = normalize(object.exprs)
  }

  lvls = rownames(object.exprs[genes[1L:n], , drop = FALSE])
  dd = reshape2::melt(as.data.frame(object.exprs[genes[1L:n], , drop = FALSE]) %>%
                        mutate(gene = genes[1L:n]))
  colnames(dd) = c("gene", "sample", "count")

  dd$xs = as.factor(metadata[as.character(dd$sample), xs])
  if (is.null(intgroup) != T) {
    dd$intgroup = as.factor(metadata[as.character(dd$sample), intgroup])

  }
  dd$gene = factor(dd$gene, levels = lvls)

  if (is.null(intgroup)) {
    st = stat_summary(fun = mean, fun.min = mean, fun.max = mean, geom = "line", colour = "#555555", lwd = .5, aes(group = 1))
  } else {
    st = stat_summary(fun = mean, fun.min = mean, fun.max = mean, geom = "line", colour = "#555555", lwd = .5, aes(group = intgroup))
  }

  if (is.null(pairs) ) {
    p = ggplot(dd, aes(x = xs, y = count, color = intgroup)) +
      geom_point(size = point.size, position = position_jitter(width = 0.1)) +
      facet_wrap(~gene, scales = "free_y", ncol = ncol) +
      xlab(xsLab) + ylab(type) +
      # scale_y_continuous(trans='log10') +
      theme(aspect.ratio=1,
            axis.text.x = element_text(angle = 90, vjust=0.5, hjust=1),
            text = element_text(size = font.size))
    p = p + st
    if (!is.null(intgroup)) {
      p = p + scale_colour_manual(values = col.p, name = eval(intgroup))
    }
    # print(p)
  }
  if (!is.null(pairs) ) {
    dd[, pairs] = as.factor(metadata[as.character(dd$sample), pairs])
    p = ggplot(dd, aes_string(x = "xs", y = "count", group = pairs)) +
      geom_point(size = .5) + facet_wrap(~gene, scales = "free_y", ncol = ncol) +
      xlab(xsLab) + ylab(type) + geom_line(size = .2) +
      theme(aspect.ratio=1,
            axis.text.x = element_text(angle = 90, vjust=0.5, hjust=1, size = rel(.8)))
    # print(p)
  }
  suppressWarnings(p)
}

## >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
## Enrichment profiles
## >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
es_barplot_top = function(query, labels, title) {

  query.nbr = length(query)

  es.gene = lapply(b.a.l, function(x){
    df = x[x$GENE_NAME %in% query, ]
    df.agg = aggregate(ES ~ PROBE_NAME, df, max)
    df.max =  merge(df.agg, df)
    df.max = droplevels(df.max)
    df.max = arrange(df.max, desc(ES))
    df.max$TISSUE = factor(df.max$TISSUE, levels = unique(df.max[order(df.max$ES, decreasing = F), ]$TISSUE))
    return(df.max)
  })
  df.rbind = as.data.frame(data.table::rbindlist(es.gene))
  # df.rbind.top = df.rbind %>%
  #   group_by(GENE_NAME) %>%
  #   # filter(ES > 0)
  #   top_n(n = 20, wt = ES) %>% arrange(desc(ES))
  df.rbind.top = df.rbind %>%
    dplyr::group_by(GENE_NAME) %>%
    dplyr::arrange(desc(ES), .by_group = TRUE) %>% filter(row_number()<=20)

  df.rbind.top$TISSUE = paste0(df.rbind.top$TISSUE, " | ", df.rbind.top$RANK)

  scale_x_reordered <- function(..., sep = "___") {
    reg <- paste0(sep, ".+$")
    ggplot2::scale_x_discrete(labels = function(x) gsub(reg, "", x), ...)
  }

  reorder_within <- function(x, by, within, fun = mean, sep = "___", ...) {
    new_x <- paste(x, within, sep = sep)
    stats::reorder(new_x, by, FUN = fun)
  }

  pl = ggplot(df.rbind.top, aes(reorder_within(TISSUE, ES, GENE_NAME), ES, colour = GROUP, fill = GROUP)) +
    geom_bar(stat = 'identity') +
    scale_x_reordered() +
    coord_flip() +
    facet_wrap(.~ GENE_NAME, ncol=3, scales = "free", labeller=labeller(GENE_NAME = labels)) +
    theme(axis.text.x=element_text(angle=90,hjust=1,vjust=.5,colour='gray50')) +
    scale_fill_manual(values = body.atlas.col) + scale_colour_manual(values = body.atlas.col) +
    ylab("Enrichment Score") + xlab("") + ggtitle(title) +
    mytheme(base_size = 5) +
    theme(legend.position="bottom",
          strip.text = element_text(size = rel(1.1)),
          plot.title = element_text(hjust = .5, face = "bold", colour = "black", size = rel(2))) +
    guides(fill=guide_legend(nrow=3))
  return(list(df = df.rbind.top, pl = pl))
}

## >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Dispersions plot
# https://galaxy.pasteur.fr/datasets/264986c34bc06a86/display/?preview=True
## >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
dispersionsPlot = function(dds){

  # dispersions plot
  d <- as.data.frame(mcols(dds)[,c("baseMean", "dispGeneEst", "dispFit", "dispersion")])
  d <- d[which(d$baseMean > 0),]
  d <- data.frame(baseMean=rep(d$baseMean, 3),
                  value=c(d$dispGeneEst, d$dispersion, d$dispFit),
                  variable=factor(rep(c("dispGeneEst", "dispersion", "dispFit"), each=nrow(d)),
                                  levels=c("dispGeneEst", "dispersion", "dispFit")))
  p1 <- ggplot(d, aes(x=.data$baseMean, y=.data$value, colour=.data$variable)) +
    geom_point(size=0.01) +
    scale_x_continuous(trans = log10_trans(),
                       breaks = trans_breaks("log10", function(x) 10^x),
                       labels = trans_format("log10", math_format())) +
    scale_y_continuous(trans = log10_trans(),
                       breaks = trans_breaks("log10", function(x) 10^x),
                       labels = trans_format("log10", math_format())) +
    ylab("Dispersion") +
    xlab("Mean of normalized counts") +
    scale_colour_manual(
      values=c("Black", "#33BBEE", "#e41a1c"),
      breaks=c("dispGeneEst", "dispersion", "dispFit"),
      labels=c("Estimate", "Final", "Fit"),
      name="") +
    guides(colour = guide_legend(override.aes = list(size=2))) +
    ggtitle("Dispersions") +
    mytheme(base_size = 6) +
    theme(legend.position = c(.8, .2))

  # diagnostic of log normality
  disp <- mcols(dds)$dispGeneEst
  disp <- disp[!is.na(disp)]
  disp <- log(disp)
  mean.disp <- mean(disp,na.rm=TRUE)
  sd.disp <- sd(disp,na.rm=TRUE)
  d <- data.frame(disp)
  p2 <- ggplot(data=d, aes(x=.data$disp)) +
    geom_histogram(bins=80, fill = "grey", colour = "black", lwd = .1, aes(y=.data$..density..)) +
    scale_y_continuous(expand=expand_scale(mult=c(0.01, 0.05))) +
    xlab("Feature dispersion estimate") +
    ylab("Density") +
    ggtitle("log-normality of estimated dispersion") +
    stat_function(fun = dnorm, args = list(mean = mean.disp, sd = sd.disp)) +
    mytheme(base_size = 6)

  p3 = ggqqplot(data.frame(disp), x = "disp", size = .1) +
    mytheme_grid(base_size = 6) +
    ggtitle("qq of estimated dispersion") +
    xlab("Theoretical Quantiles") +
    ylab("Sample Quantiles")

  plot_grid(p1, p2, p3, nrow = 1, align = "h")
}

## >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
## cnetplot form Clusterprofiler (modified)
## >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# convert a list of gene IDs to igraph object.
list2graph <- function(inputList) {
  x <- list2df(inputList)
  g <- graph.data.frame(x, directed=FALSE)
  return(g)
}

list2df <- function(inputList) {
  ldf <- lapply(1:length(inputList), function(i) {
    data.frame(categoryID=rep(names(inputList[i]),
                              length(inputList[[i]])),
               Gene=inputList[[i]])
  })

  do.call('rbind', ldf)
}

pathnet = function(obj,
                   showCategory = 5,
                   foldChange   = NULL,
                   layout = "kk",
                   node_size = 4,
                   base_size = 8,
                   title_size_rel = 1.5,
                   .gc.anno = .gc.anno) {
  library(ggraph)
  library(igraph)

  geom_edge <- geom_edge_link

  geneSets = strsplit(obj$geneID, "/")
  names(geneSets) = obj$Description
  entrez.to.symbol = lapply(geneSets, function(x) {
    x = .gc.anno$GENE_SYMBOL[match(x, .gc.anno$ENTREZ)]
  })

  if (!all(is.na(unlist(entrez.to.symbol)))) {
    geneSets = entrez.to.symbol
  }

  if (length(geneSets) <= showCategory) {
    geneSets = geneSets[1:length(geneSets)]
  } else {
    geneSets = geneSets[1:showCategory]
  }

  g <- list2graph(geneSets)

  size <- sapply(geneSets, length)
  V(g)$size <- min(size)/2

  n <- length(geneSets)
  V(g)$size[1:n] <- size

  edge_layer <- geom_edge(alpha=.8, width = .2, colour='darkgrey')

  if (!is.null(foldChange)) {
    fc <- foldChange[V(g)$name[(n+1):length(V(g))]]
    V(g)$color <- NA
    V(g)$color[(n+1):length(V(g))] <- fc

    p = ggraph(g, layout=layout, circular = FALSE) +
      edge_layer +
      geom_node_point(aes_(color = ~ as.numeric(as.character(color)), size = ~ size)) +
      scale_size(name = "No. of genes", range=c(2, 6), breaks=unique(round(seq(min(size), max(size), length.out=4)))) +
      scale_colour_gradient2(name = "Standardized combined\nlog2 fold change", low = "#0077BB", mid = "white", high = "#CC3311", na.value = "grey60")
  } else {
    V(g)$color <- "#B3B3B3"
    V(g)$color[1:n] <- "#E5C494"
    p <- ggraph(g, layout=layout, circular = FALSE) +
      edge_layer +
      geom_node_point(aes_(color=~I(color), size=~size))
  }

  p = p + geom_node_text(aes_(label=~name), repel=TRUE, size = node_size)
  p = p + mytheme(base_size = base_size)
  p = p + theme(
    legend.position="top",
    axis.line=element_blank(),
    axis.text.x=element_blank(),
    axis.text.y=element_blank(),
    axis.ticks=element_blank(),
    axis.title.x=element_blank(),
    axis.title.y=element_blank(),
    panel.background=element_blank(),
    panel.border=element_blank(),
    panel.grid.major=element_blank(),
    panel.grid.minor=element_blank(),
    plot.background=element_blank(),
    plot.title = element_text(hjust = .5, face = "plain", colour = "black", size = rel(title_size_rel))
    )
  p = p + guides(size = guide_legend(order=1))

  return(p)
}

pathnet_meta = function(.cP.res, .lvls,
                        .meta.res = meta.res,
                        .gc.anno = gc.anno,
                        .categories = 3,
                        .node.size = 3,
                        .base.size = 8,
                        .title.size.rel = 1.5) {

  meta.res.tmp = .meta.res
  names(meta.res.tmp) = paste0(1:length(names(meta.res.tmp)), "_", names(meta.res.tmp))

  lfc.l = lapply(names(meta.res.tmp), function(ctrst){
    ctrst.entrez = .cP.res@geneClusters[[ctrst]]
    ctrst.symbol = .gc.anno$GENE_SYMBOL_DUPL_MARKED[match(ctrst.entrez, .gc.anno$ENTREZ)]
    m.res.ctrst = meta.res.tmp[[ctrst]]
    m.res.ctrst = m.res.ctrst[m.res.ctrst$GENE_SYMBOL %in% ctrst.symbol, ]
    setNames( m.res.ctrst$est_effect, nm = m.res.ctrst$GENE_SYMBOL)
  })
  names(lfc.l) = names(meta.res.tmp)


  cP.res.ctrsts = as.character(unique(.cP.res@compareClusterResult$Cluster))
  cP.res.ctrsts = naturalsort(cP.res.ctrsts)

  # ctrst = "5_12h"
  pathnet.l = lapply(.lvls, function(ctrst) {
    cP.res.df = .cP.res@compareClusterResult
    cP.res.df = cP.res.df[cP.res.df$Cluster %in% ctrst, ]
    pathnet(
      obj = cP.res.df,
      showCategory = .categories,
      foldChange = lfc.l[[ctrst]],
      .gc.anno = .gc.anno,
      node_size = .node.size,
      base_size = .base.size,
      title_size_rel = .title.size.rel
    ) + ggtitle(gsub("._", "", ctrst))
  })
}

