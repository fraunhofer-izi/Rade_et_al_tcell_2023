library(naturalsort)

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
exprs_med_melt = function(.exprs.mat = NULL, .ftrs = NULL, .exclude = "Th1_") {
  exprs.all = .exprs.mat[subset(.ftrs, CLUSTER != "0100")$GENE_SYMBOL, ]
  exprs.all = z_score(exprs.all)

  if (any(.ftrs$CLUSTER == "0100")) {
    exprs.0100 = .exprs.mat[subset(.ftrs, CLUSTER == "0100")$GENE_SYMBOL, ]
    exprs.0100 = exprs.0100[, !grepl(.exclude, colnames(exprs.0100))]
    exprs.0100 = z_score(exprs.0100)

    df = rbind(
      reshape2::melt(as.matrix(exprs.all)),
      reshape2::melt(as.matrix(exprs.0100))
    )
  } else {
    exprs.0100 = NULL
    df = reshape2::melt(as.matrix(exprs.all))
  }

  colnames(df) = c("GENE_SYMBOL", "GROUP", "EXPRS")
  df$HOURS = gsub(".*\\_", "", df$GROUP)
  df$CELLTYPE = gsub("\\_.+", "", df$GROUP)
  df$CLUSTER = .ftrs$CLUSTER[match(df$GENE_SYMBOL, .ftrs$GENE_SYMBOL)]
  lvls = unique(df$HOURS)
  lvls = lvls[naturalorder(as.numeric(gsub("h", "", lvls)))]
  df$HOURS = factor(df$HOURS, levels = lvls)

  return(list(df = df, exprs.all = exprs.all, exprs.0100 = exprs.0100))
}

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Information Coefficient [IC]
# https://github.com/yanwu2014/swne
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Compute Information Coefficient [IC]
# Pablo Tamayo Dec 30, 2015
#
# @param x Input vector x
# @param y Input vector y
# @param n.grid Gridsize for calculating IC
#
# @return Mutual information between x and y
#

MutualInf <-  function(x, y, n.grid = 25) {
  x.set <- !is.na(x)
  y.set <- !is.na(y)
  overlap <- x.set & y.set

  x <- x[overlap] +  0.000000001*runif(length(overlap))
  y <- y[overlap] +  0.000000001*runif(length(overlap))

  if (length(x) > 2) {
    delta = suppressWarnings(c(MASS::bcv(x), MASS::bcv(y)))
    rho <- cor(x, y)
    rho2 <- abs(rho)
    delta <- delta*(1 + (-0.75)*rho2)
    kde2d.xy <- MASS::kde2d(x, y, n = n.grid, h = delta)
    FXY <- kde2d.xy$z + .Machine$double.eps
    dx <- kde2d.xy$x[2] - kde2d.xy$x[1]
    dy <- kde2d.xy$y[2] - kde2d.xy$y[1]
    PXY <- FXY/(sum(FXY)*dx*dy)
    PX <- rowSums(PXY)*dy
    PY <- colSums(PXY)*dx
    HXY <- -sum(PXY * log(PXY))*dx*dy
    HX <- -sum(PX * log(PX))*dx
    HY <- -sum(PY * log(PY))*dy
    PX <- matrix(PX, nrow=n.grid, ncol=n.grid)
    PY <- matrix(PY, byrow = TRUE, nrow=n.grid, ncol=n.grid)
    MI <- sum(PXY * log(PXY/(PX*PY)))*dx*dy
    IC = MI
    IC <- sign(rho) * sqrt(1 - exp(- 2 * MI))
    if (is.na(IC)) IC <- 0
  } else {
    IC <- 0
  }
  return(IC)
}

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Ammotate metagenes by median "expression" peak and reorder H and W
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
anno_metagenes = function(obj, rank, reorder = TRUE, group = NULL, .order = NULL) {
  nmf.obj = obj$estim.r
  se.obj = obj$se

  if(!is.null(obj$estim.r.random)) {
    estim.r.random = obj$estim.r.random
  } else {
    estim.r.random = NULL
  }

  h = NMF::coef(nmf.obj$fit[[rank]])
  colnames(h) = colData(se.obj)[[group]]
  groups = unique(colnames(h))
  rownames(h) = paste0("M_",rep(1:nrow(h),1))

  h.grp.med = reshape2::melt(h) %>%
    dplyr::group_by(Var1, Var2) %>%
    dplyr::summarize(median = median(value)) %>%
    data.frame() %>% `colnames<-`(c("metagene", "group", "median"))

  metagenes = h.grp.med %>%
    dplyr::group_by(metagene) %>%
    dplyr::mutate(rank = rank(dplyr::desc(median))) %>%
    dplyr::filter(rank == 1) %>%
    dplyr::summarize(label = paste0("M_", paste0(group, collapse = '_'))) %>%
    data.frame()
  metagenes = metagenes[naturalsort::naturalorder(metagenes$label), ]

  if (!is.null(.order)) {
    x = gsub("\\_.*", "", gsub("M_", "", metagenes$label))
    metagenes = metagenes[order(match(x,.order)), ]

  }

  # metagenes = metagenes[naturalsort::naturalorder(metagenes$label), ]
  metagenes$prim_order = as.numeric(as.character(rownames(metagenes)))
  metagenes$metagene = paste0("M_",rep(1:nrow(metagenes),1))

  suffix = gsub("M_", "", metagenes$metagene)
  for (i in 1:nrow(metagenes)) {
    binary = ifelse(suffix %in% gsub("M_", "", metagenes[i, "metagene"]), 1, 0)
    metagenes[i, "binary"] = paste(binary,collapse=" ")
  }

  if (reorder == TRUE) {
    NMF::coef(nmf.obj$fit[[rank]]) = NMF::coef(nmf.obj$fit[[rank]])[metagenes$prim_order, ]
    NMF::basis(nmf.obj$fit[[rank]]) = NMF::basis(nmf.obj$fit[[rank]])[, metagenes$prim_order]
    metagenes$prim_order = NULL
    return(list(estim.r = nmf.obj, estim.r.random = estim.r.random, se = se.obj, anno = metagenes))
  } else {
    metagenes$prim_order = NULL
    return(list(estim.r = nmf.obj, estim.r.random = estim.r.random, se = se.obj, anno = metagenes))
  }
}

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Plot metagene profile (H) over time
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
kinetic_mg_profile = function(obj, rank, metag.col, font.size = 8) {

  nmf.obj = obj$estim.r
  se.obj = obj$se
  h = coef(nmf.obj$fit[[rank]])
  colnames(h) = colData(se.obj)$HOURS
  rownames(h) = paste0("M",rep(1:nrow(h),1))

  for (i in 1:ncol(h)) {
    h[ ,i] = h[ ,i] / sum(h[ ,i])
  }

  h.melt = reshape2::melt(as.matrix(h))
  colnames(h.melt) = c("METAGENE", "HOURS", "WEIGHT")
  lvls = as.character(naturalsort(unique(h.melt$HOURS)))
  h.melt$HOURS = factor(h.melt$HOURS, levels = lvls)

  ggplot(h.melt, aes(x = HOURS, y = WEIGHT, group = METAGENE, colour = METAGENE)) +
    stat_summary(fun = median, fun.min = median, fun.max = median,
                 geom = "line", size = .5) +
    geom_pointrange(mapping = aes(x = HOURS, y = WEIGHT),
                    stat = "summary",
                    fun.min = function(z) {quantile(z,0.25)},
                    fun.max = function(z) {quantile(z,0.75)},
                    fun = median, size = .2) +
    mytheme(base_size = font.size) +
    scale_color_manual(values = metag.col) +
    theme(
      axis.title.y=element_blank(),
      axis.title.x=element_blank())
}

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Coef (H) heatmap
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
coef_heatmap = function(.obj = NULL, .rank = NULL, .anno_col = NULL, .cluster_cols = T ) {
  library(pheatmap)

  mat.h = coef(.obj$estim.r$fit[[.rank]])
  mat.h = sweep(mat.h, 2L, colSums(mat.h, na.rm =T), '/', check.margin = FALSE)
  rownames(mat.h) = paste0("M",rep(1:nrow(mat.h),1))

  hm.features = c(.anno_col)
  anno.col = data.frame(colData(.obj$se)[, hm.features, drop = F])
  anno.col.colors = degColors(anno.col)

  pheatmap(
    mat               = mat.h,
    color             = colorRampPalette(c("skyblue",  "tomato3"))(50),
    annotation_col    = anno.col,
    border_color      = NA,
    show_rownames     = T,
    show_colnames     = F,
    cluster_cols      = .cluster_cols,
    cluster_rows      = F,
    annotation_colors = anno.col.colors,
    drop_levels       = TRUE,
    fontsize          = 8,
    scale             = "none",
    treeheight_col    = 10
  )
}


# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Consensus maps (heatmaps) and diff nmf stats
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
FrobErrors = function(nmf.fit, se){
  # get all W-matrices
  W.list <- lapply(nmf.fit, basis)
  # get all H-matrices
  H.list = lapply(nmf.fit, coef)
  fes <- sapply(seq(length(W.list)), function(oi){
    fe <- norm(exprs(se) - (W.list[[oi]] %*% H.list[[oi]]), type = "F") / norm(exprs(se), type = "F")
  })
  return(fes)
}

nmf_stats = function(obj, anno, hm.fonsize = 5.5, main.title = "Title") {
  library("ComplexHeatmap")
  estim.r = obj$estim.r
  se = obj$se

  if(!is.null(obj$estim.r.random)) {
    estim.r.random = obj$estim.r.random
  }

  h.grid.l = list()

  for (i in names(estim.r$fit)) {
    consensus.m = estim.r$fit[[as.character(i)]]@consensus
    # distance and cluster
    dend = hclust(as.dist(1-consensus.m), method = "average")

    # top annotation bar
    anno.top = colData(se)[, anno, drop = F]
    anno.top.col = degColors(anno.top, tol21rainbow = tol10qualitative)
    top.anno = HeatmapAnnotation(
      df = anno.top,
      col = anno.top.col,
      show_legend = F,
      simple_anno_size = unit(.2, "cm"),
      annotation_name_gp = gpar(fontsize = hm.fonsize - 2))

    lvls = naturalsort(unique(anno.top$HOURS))
    top.anno@anno_list$HOURS@color_mapping@levels = lvls

    # cophenetic
    coph = round(estim.r$measures$cophenetic, 3)[as.integer(as.character(i)) - 1]

    # heatmap
    h = Heatmap(consensus.m,
                cluster_rows = dend,
                cluster_columns = dend,
                show_row_names = FALSE,
                show_column_names = FALSE,
                column_dend_reorder = FALSE,
                row_dend_reorder= FALSE,
                border = TRUE,
                width = unit(4.5, "cm"), height = unit(4.5, "cm"),
                column_title_gp = gpar(fontsize = hm.fonsize),
                row_title_gp = gpar(fontsize = hm.fonsize),
                col = rev(magma(100)[30:100]),
                top_annotation = top.anno,
                row_title = "Samples",
                column_title = paste0("Samples\nCophenetic coef = ", coph),
                column_title_side = "bottom",
                show_heatmap_legend = F,
                row_dend_width = unit(.2, "cm"),
                column_dend_height = unit(.2, "cm"),
                heatmap_legend_param = list(title = "",
                                            legend_height = unit(3, "cm"),
                                            grid_width = unit(.2, "cm"),
                                            labels_gp = gpar(fontsize = hm.fonsize))
    )
    # grid for heatmaps
    h.gTree = grid.grabExpr(draw(h))
    title = ggdraw() + draw_label(paste0("k=", i), size = hm.fonsize)
    h.grid = plot_grid(title, h.gTree, rel_heights = c(.1, .9), ncol = 1)
    h.grid.l[[i]] = h.grid

    if (i == names(estim.r$fit)[length(names(estim.r$fit))]) {
      anno.leg = lapply(h@top_annotation@anno_list, function(x) {
        Legend(labels = x@color_mapping@levels,
               title = x@color_mapping@name,
               legend_gp = gpar(fill = x@color_mapping@colors),
               labels_gp = gpar(fontsize = hm.fonsize),
               title_gp = gpar(fontsize = hm.fonsize))})

      anno.leg$HOURS@grob$children[[2]]$children[[2]]$gp$fill = anno.top.col$HOURS

      dist.leg = color_mapping_legend(
        h@matrix_color_mapping,
        plot = F, title = "",
        labels_gp = gpar(fontsize = hm.fonsize))
      if (length(anno) == 1) {
        anno.leg.all = packLegend(dist.leg, anno.leg[[1]])
      }
      if (length(anno) == 2) {
        anno.leg.all = packLegend(dist.leg, anno.leg[[1]], anno.leg[[2]])
      }
      if (length(anno) == 3) {
        anno.leg.all = packLegend(dist.leg, anno.leg[[1]], anno.leg[[2]], anno.leg[[3]])
      }
      leg.gTree = grid.grabExpr(draw(anno.leg.all))

      consensus.maps = plot_grid(plotlist = h.grid.l, nrow = 3)
      consensus.maps.all = plot_grid(consensus.maps, leg.gTree, rel_widths = c(.8, .2))
    }

  }

  nmf.stats = estim.r$measures %>%
    dplyr::select(rank, cophenetic, dispersion, silhouette.consensus) %>%
    dplyr::mutate(run = "actual")

  coph.p = ggplot(nmf.stats, aes(x = rank, y = cophenetic)) +
    geom_line() +
    geom_point() +
    mytheme_grid(base_size = hm.fonsize) +
    theme(aspect.ratio = 1, axis.text.y = element_text(angle = 90, hjust = .5)) +
    ylab("Cophenetic coefficient")

  silh = ggplot(nmf.stats, aes(x = rank, y = silhouette.consensus)) +
    geom_line() +
    geom_point() +
    mytheme_grid(base_size = hm.fonsize) +
    theme(aspect.ratio = 1,axis.text.y = element_text(angle = 90, hjust = .5)) +
    ylab("Avg silhouette width")

  stats.grid = plot_grid(plot_grid(coph.p, NULL, silh, rel_widths = c(1, .03, 1), nrow = 1), NULL, rel_widths = c(.8, .2))

  if(!is.null(obj$estim.r.random)) {
    nmf.stats.r = estim.r.random$measures %>%
      dplyr::select(rank, cophenetic, dispersion, silhouette.consensus) %>%
      dplyr::mutate(run = "randomized")

    coph.p = ggplot(rbind(nmf.stats, nmf.stats.r), aes(x = rank, y = cophenetic)) +
      geom_line() +
      geom_point() +
      mytheme_grid(base_size = hm.fonsize) +
      theme(panel.spacing = unit(1, "lines")) +
      # theme(aspect.ratio = 1, axis.text.y = element_text(angle = 90, hjust = .5)) +
      ylab("Cophenetic coefficient") +
      ylim(c(NA, 1)) +
      facet_wrap( ~ run, scales = "free_y")

    silh = ggplot(rbind(nmf.stats, nmf.stats.r), aes(x = rank, y = silhouette.consensus)) +
      geom_line() +
      geom_point() +
      mytheme_grid(base_size = hm.fonsize) +
      theme(panel.spacing = unit(1, "lines")) +
      ylab("Avg silhouette width") +
      ylim(c(NA, 1)) +
      facet_wrap( ~ run, scales = "free_y")

    stats.grid = plot_grid(plot_grid(coph.p, NULL, silh, rel_widths = c(1, .1, 1), nrow = 1), NULL, rel_widths = c(.9, .05))
  }

  # now add the title
  title = ggdraw() +
    draw_label(
      main.title,
      fontface = 'bold',
      x = .32,
      hjust = 0
    )

  res = plot_grid(title, consensus.maps.all, NULL, stats.grid, nrow = 4,rel_heights = c(.04, .8, .03, .24), align = "vh")
  return(res)
}

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Function to order binary data matrix
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
orderBinary <- function(matrix) {
  col.sum <- apply(matrix, 2, sum)
  unlist(sapply(unique(col.sum), function(s) which(col.sum == s)))
}

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Group average
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
aveByGrp = function(.obj, .ave = "mean") {
  grp.ave.l = list()
  grp.sd.l = list()
  for (i in paste0("^", unique(colnames(.obj)), "$")) {

    grp = .obj[, grepl(i , colnames(.obj))]

    if (.ave == "mean") {
      grp.ave = data.frame(rowMeans(grp), row.names = rownames(grp))
    } else {
      grp.ave = data.frame(rowMedians(grp), row.names = rownames(grp))
    }
    colnames(grp.ave) = unique(colnames(grp))
    grp.ave.l[[unique(colnames(grp))]] = grp.ave

  }
  grp.ave.df = do.call("cbind", grp.ave.l)
  return(grp.ave.df)
}

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Average by group -> Arrays scaled by row based to quantile from RNA-Seq
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
aveByGrp_qscaled = function(.obj, .tp = "_0h|_12h|_24h|_48h|_72h", .qu = NULL, .type = "mean", counts = NULL, anno = NULL) {
  exprs.cpm = assays(.obj)[[counts]]
  colnames(exprs.cpm) = gsub(".+\\.", "", colnames(exprs.cpm))
  aveByGroup = aveByGrp(.obj = exprs.cpm, .ave = .type)

  if (!is.null(anno)){
    rownames(aveByGroup) = gc.anno$GENE_SYMBOL_DUPL_MARKED[match(rownames(aveByGroup), gc.anno$ENSEMBL_ID_ABBR)]
  }
  colnames(aveByGroup) = gsub("05h", "0.5h", colnames(aveByGroup))

  aveByGroup.th2 = aveByGroup[, grepl("Th2_", colnames(aveByGroup))]
  aveByGroup.th1 = aveByGroup[, grepl("Th1_", colnames(aveByGroup))]
  aveByGroup.seq = aveByGroup[, !grepl("Th1_|Th2_", colnames(aveByGroup))]

  scaling.fctr.th2 = rowQuantiles(as.matrix(aveByGroup.th2), probs = .qu) /
    rowQuantiles(as.matrix(aveByGroup.seq), probs = .qu)
  aveByGroup.th2.sc = aveByGroup.th2/ scaling.fctr.th2

  aveByGroup.seq.tmp = aveByGroup.seq[, grepl(.tp, colnames(aveByGroup.seq))]
  scaling.fctr.th1 = rowQuantiles(as.matrix(aveByGroup.th1), probs = .qu) /
    rowQuantiles(as.matrix(aveByGroup.seq.tmp), probs = .qu)
  aveByGroup.th1.sc = aveByGroup.th1/ scaling.fctr.th1

  aveByGroup = cbind(
    aveByGroup.seq,
    aveByGroup.th2.sc,
    aveByGroup.th1.sc
  )
  return(aveByGroup)
}


# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Validates gene embeddings by plotting cluster logFC vs gene factor loading logFC
# https://github.com/yanwu2014/swne
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
EmbeddedFeatures <- function(W, norm.counts, genes.embed, sample.groups,
                             min.cluster.logfc = 1.5, min.factor.logfc = 1.5,
                             font.size = 12, label.size = 3, title = NULL,
                             label = T) {

  if(!all(rownames(W) == rownames(norm.counts))) {
    stop("Rownames of W must match rownames of norm.counts")
  }

  colnames(W) = paste0("F",rep(1:ncol(W),1))

  basis.norm = W/colSums(W); eps = 1e-4
  gene.factor.logfc = log2(apply(basis.norm, 1, function(x) {
    max.i <- which.max(x)
    (x[[max.i]] + eps)/(mean(x[-1*max.i]) + eps)
  }))

  gene.cell.logfc = log2(apply(norm.counts, 1, function(x) {
    cl.avg = tapply(x, sample.groups, mean)
    cl.max = names(cl.avg[which.max(cl.avg)])
    (mean(x[sample.groups == cl.max]) + eps)/(mean(x[sample.groups != cl.max]) + eps)
  }))

  gene.cell.logfc.cl.max = apply(norm.counts, 1, function(x) {
    cl.avg = tapply(x, sample.groups, mean)
    cl.max = names(cl.avg[which.max(cl.avg)])
    cl.max
  })

  gene.cell.logfc.cl.ave = apply(norm.counts, 1, function(x) {
    log10(mean(x) + 1)
  })

  gene.logfc.df = data.frame(CLUSTER = gene.cell.logfc, FACTOR = gene.factor.logfc,
                             CL_MAX = gene.cell.logfc.cl.max, MEAN = gene.cell.logfc.cl.ave)
  gene.logfc.df$EMBEDDED = factor(rownames(gene.logfc.df) %in% genes.embed)
  gene.logfc.df = gene.logfc.df[order(gene.logfc.df$EMBEDDED),]
  # rownames(gene.logfc.df) = paste0(rownames(gene.logfc.df),"(", gene.cell.logfc.grp.names, ")")

  gene.logfc.label.point = subset(gene.logfc.df, EMBEDDED == "TRUE")
  # gene.logfc.label.name = subset(gene.logfc.label.point, FACTOR > min.factor.logfc)
  # gene.logfc.label.name = subset(gene.logfc.label.name, CLUSTER > min.cluster.logfc)
  gene.logfc.label.name = subset(gene.logfc.label.point, FACTOR > 3)
  gene.logfc.label.name = subset(gene.logfc.label.name, CLUSTER > 1)
  gene.logfc.label.name$NAME = paste0(rownames(gene.logfc.label.name),"_(", gene.logfc.label.name$CL_MAX, ")")

  gene.logfc.label.point = gene.logfc.label.point %>%
    dplyr::mutate(
      LITERATURE = case_when(
        rownames(gene.logfc.label.point) %in% unlist(act.mrk) ~ "Activation",
        rownames(gene.logfc.label.point) %in% unlist(anergy.mrk) ~ "Anergy",
        TRUE ~ "No"
      )
    )
  gene.logfc.label.point = gene.logfc.label.point %>%
    dplyr::mutate(
      COLOR = case_when(
        LITERATURE == "Activation" ~ "#DDAA33",
        LITERATURE == "Anergy" ~ "#EE6677",
        LITERATURE == "No" ~ "darkred"
      )
    )
  # keywords = c("Activation", "Anergy", "No")
  # jColors = data.frame(KEYWORD = keywords,
  #                      COLOR = c("#DDAA33","black","darkred"),
  #                      stringsAsFactors = F)

  gg.obj = gene.logfc.df %>% ggplot(aes(FACTOR, CLUSTER)) +
    geom_hex(bins = 50, alpha = 1) +
    xlab("Max metagene logFC") + ylab("Max sampleGroup logFC") +
    scale_fill_gradient(low = "skyblue", high = "tomato3") +
    mytheme + theme(text = element_text(size = font.size), aspect.ratio = 1,
                    plot.title = element_text(vjust=-.5, face = "bold"))

  gg.obj = gg.obj + geom_point(data = gene.logfc.label.point, size = 2,
                               aes(color = LITERATURE)) +
    scale_color_manual(values = c(Activation = "#DDAA33", Anergy = "#EE6677",
                                  No = "darkred"))
  if (label == T) {
    gg.obj = gg.obj + ggrepel::geom_text_repel(data = gene.logfc.label.name,
                                               mapping = aes(FACTOR, CLUSTER, label = NAME),
                                               size = label.size, box.padding = 0.15)
  }
  gg.obj = gg.obj + geom_vline(xintercept = min.factor.logfc,
                               linetype = "dotted", size = 1.25, color = "black") +
    geom_hline(yintercept = min.cluster.logfc, linetype = "dotted", size = 1.25, color = "black") +
    ggtitle(title)

  gg.obj
}

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
  return(apply(mat, 2, function(y) (rank(y))/ length(y)))
}

# rank_norm = function(mat) {
#   return(apply(mat, 2, function(y) 10000*(rank(y) -1)/ length(y)))
# }

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

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Temporal Correlation Distance
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#check if series have equal length, a requisite of some functions
.check.equal.length.ts <- function(x,y) {
  if (length(x) != length(y)) {
    stop("Time series must have the same length")
  }
}

corrtemporder1 <- function (x, y) {
  p <- length(x)
  sum((x[2:p] - x[1:(p-1)]) * (y[2:p] - y[1:(p-1)])) / ( sqrt( sum((x[2:p] - x[1:(p-1)])^2) ) * sqrt( sum((y[2:p] - y[1:(p-1)])^2) ))
}
##CHOUAKRIA-DOUZAL
cort = function( x, y, k=2) {
  .check.equal.length.ts(x,y)
  corrt = corrtemporder1(x,y)
  typedist = 1
  typedist <- as.numeric( dist(rbind(x,y)) )
  # typedist = dtw(x,y, window.type = "sakoechiba", window.size = 2, distance.only=T)$distance
  (2/( 1+ exp(k*corrt)))*typedist
}

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Autocorrelaton
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#weighted distance of acf and pacf coefficients
.internal.autocorr.dist <- function(rhox, rhoy, p=NULL, omega=NULL) {
  if ( length(rhox) != length(rhoy) ) {#check compatible coefficient vectors
    stop("The amount of autocorrelation coefficients must be the same, maybe lag.max greater than the length of one of the series")
  }
  if (is.null(omega)) { #if there is no weighting matrix
    if (!is.null(p)) { #check if there is gemoetrical decay parameter
      omega <- diag(p*(1-p)**(1:length(rhox)))
    }
    else { #no weightinh matrix and no geomtrical decay parameter, use identity matrix
      omega <- diag(length(rhox))
    }
  }
  sqrt(t(rhox - rhoy) %*% omega %*% (rhox - rhoy)) #weighted euclidean distance
}


diss.ACF <- function(x, y ,  p=NULL,  omega=NULL, lag.max=50) {
  rhox <- acf(x, lag.max=lag.max, plot=FALSE)$acf[-1]
  rhoy <- acf(y, lag.max=lag.max, plot=FALSE)$acf[-1]
  .internal.autocorr.dist( rhox, rhoy, p, omega)
}

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# short time series (sts) distance measure.
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
STSDistance <- function(x, y, tx=NULL, ty=NULL) {

  # If no index is specified then evenly samples series are assumed.
  if (is.null(tx) & is.null(ty)) {
    tx <- c(1:length(x))
    ty <- tx
  }
  if (is.null(tx)) {
    tx <- ty
  }
  if (is.null(ty)) {
    ty <- tx
  }

  if (class(try(STSInitialCheck(x, y, tx, ty))) == "try-error") {
    return(NA)
  } else {

    # The STS distance is calculated.
    d <- sqrt(sum((diff(x) / diff(tx) - diff(y) / diff(ty)) ^ 2))
    return(d)
  }

}

#  This function checks for possible initial errors:
STSInitialCheck <- function(x, y, tx, ty) {

  if (! is.numeric(x) | ! is.numeric(y)) {
    stop('The series must be numeric', call.=FALSE)
  }
  if (! is.vector(x) | ! is.vector(y)) {
    stop('The series must be univariate vectors', call.=FALSE)
  }
  if (length(x) <= 1 | length(y) <= 1) {
    stop('The series must have a more than one point', call.=FALSE)
  }
  if (length(x) != length(y)) {
    stop('Both series must have the same length', call.=FALSE)
  }
  if (any(is.na(x)) | any(is.na(y))) {
    stop('There are missing values in the series', call.=FALSE)
  }
  if (! missing(tx) & ! missing(ty)) {
    if (any(tx<=0) | any(ty<=0)) {
      stop('The temporal indice must always be positive', call.=FALSE)
    }
    if (any(diff(tx) != diff(ty))) {
      stop('The sampling rate must be equal in both series', call.=FALSE)
    }
    if (any(diff(tx) <= 0) | any(diff(ty) <= 0)) {
      stop('The temporal index must be ascending.', call.=FALSE)
    }

    if ((length(tx) != length(x)) | (length(ty) != length(y))) {
      stop('The length of the time indice must be equal to the length of the series', call.=FALSE)
    }
  }
}

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
logit_scale = function(matrix, k = 8) {
  sc = apply(matrix, 1, function(x) {
    # x = (x- min(x)) /(max(x)-min(x))
    # norm = 1 / (1 + exp(-k * (x-0.5)))
    q <- as.numeric(quantile(x, .95))
    # x / q
    norm <- 2 / (1 + exp((-2) * x / q)) - 1
  })
  return(t(sc))
}

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# 3 Metagenes to 4 Metagenes
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
pseudo_3_to_4_mg_th1 = function(obj) {
  cl = do.call(rbind, strsplit(as.character(obj$cluster), split = ""))
  cl = data.frame(cl, row.names = rownames(obj))

  cl.s = cl[cl$X2 == 0, ]
  cl = cl[cl$X2 != 0, ]

  cmb.0.1 = data.frame(M1 = cl.s[, 1], M2 = 0, M3 = 0, M4 = cl.s[, ncol(cl.s)], row.names = rownames(cl.s))
  cmb.0.1 = apply(cmb.0.1, 1, function(x) {paste(x, collapse="")})

  cmb.0.2 = data.frame(M1 = cl.s[, 1], M2 = 1, M3 = 0, M4 = cl.s[, ncol(cl.s)], row.names = rownames(cl.s))
  cmb.0.2 = apply(cmb.0.2, 1, function(x) {paste(x, collapse="")})

  cmb.1 = data.frame(M1 = cl[, 1], M2 = 1, M3 = 1, M4 = cl[, ncol(cl)], row.names = rownames(cl))
  cmb.1 = apply(cmb.1, 1, function(x) {paste(x, collapse="")})

  cmb.2 = data.frame(M1 = cl[, 1], M2 = 1, M3 = 0, M4 = cl[, ncol(cl)], row.names = rownames(cl))
  cmb.2 = apply(cmb.2, 1, function(x) {paste(x, collapse="")})

  cmb.3 = data.frame(M1 = cl[, 1], M2 = 0, M3 = 1, M4 = cl[, ncol(cl)], row.names = rownames(cl))
  cmb.3 = apply(cmb.3, 1, function(x) {paste(x, collapse="")})

  obj.0.1 = obj; obj.0.1$cluster = cmb.0.1[match(rownames(obj.0.1), names(cmb.0.1))]
  obj.0.2 = obj; obj.0.2$cluster = cmb.0.2[match(rownames(obj.0.2), names(cmb.0.2))]
  obj.1 = obj; obj.1$cluster = cmb.1[match(rownames(obj.1), names(cmb.1))]
  obj.2 = obj; obj.2$cluster = cmb.2[match(rownames(obj.2), names(cmb.2))]
  obj.3 = obj; obj.3$cluster = cmb.3[match(rownames(obj.3), names(cmb.3))]

  obj.m = rbind(obj.0.1, obj.0.2, obj.1, obj.2, obj.3)
  obj.m = obj.m[!is.na(obj.m$cluster), ]

  return(obj.m)
}

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Find optimal number of factors
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#' KL divergence. Pseudocounts added to avoid NAs
kl_div <- function(x, y, pseudocount = 1e-12) {
  x <- x + pseudocount; y <- y + pseudocount;
  stopifnot(length(x) == length(y))
  return(x*log(x/y) - x + y)
}

PlotFactorSelection <- function(k.err.diff, font.size = 12) {
  if(is.list(k.err.diff)) {
    k.err.diff <- k.err.diff[["err"]]
  }
  err.del.df <- data.frame(y = k.err.diff, x = factor(names(k.err.diff), levels = names(k.err.diff)))

  ggplot(data = err.del.df, aes(x, y)) +
    geom_line(aes(group = 1)) + geom_point(size = 2.0) +
    theme_classic() + xlab("Number of factors") + ylab("Error reduction above noise") +
    theme(axis.text.x = element_text(hjust = 1, size = font.size, angle = 90, color = "black"),
          axis.text = element_text(size = font.size, color = "black"),
          axis.title = element_text(size = font.size, color = "black"),
          legend.position = c(0.8,0.8),
          legend.text = element_text(size = font.size),
          legend.title = element_text(size = font.size)) +
    geom_hline(yintercept = 0.0, linetype = "dashed", color = "darkred")
}

FindNumFactors = function(nmf.fit, k.range = seq(2,10,1), loss = "mkl", do.plot = TRUE){

  obj = nmf.fit$estim.r
  obj.random = nmf.fit$estim.r.random

  X = assays(nmf.fit$se)$nmf
  X.random = assays(nmf.fit$se)$random

  # get all W-matrices
  W.list <- lapply(obj$fit, basis)
  W.list.random <- lapply(obj.random$fit, basis)

  # get all H-matrices
  H.list = lapply(obj$fit, coef)
  H.list.random = lapply(obj.random$fit, coef)

  X.hat <- lapply(seq(length(W.list)), function(oi){
    W.list[[oi]] %*% H.list[[oi]]
  })
  X.hat.random <- lapply(seq(length(W.list)), function(oi){
    W.list.random[[oi]] %*% H.list.random[[oi]]
  })

  names(X.hat) = as.character(k.range)
  names(X.hat.random) = as.character(k.range)

  k.err = sapply(k.range, function(k) {

    if (loss == "mse") {
      err  <- mean((X.hat[[as.character(k)]] - X)^2)
      err.rand <- mean((X.hat.random[[as.character(k)]] - X.random)^2)
    } else if (loss == "mkl") {
      err = mean(kl_div(X.hat[[as.character(k)]], X))
      err.rand = mean(kl_div(X.hat.random[[as.character(k)]], X.random))
      # library(philentropy)
      # X.1 = c(X.hat[[as.character(k)]]) / sum(c(X.hat[[as.character(k)]]))
      # X.2 = c(X) / sum(c(X))
      # X.1.r = c(X.hat.random[[as.character(k)]]) / sum(c(X.hat.random[[as.character(k)]]))
      # X.2.r = c(X.random) / sum(c(X.random))
      # err = KL(rbind(X.1, X.2))
      # err.rand = KL(rbind(X.1.r, X.2.r))
    }

    return(rbind(err, err.rand))
  })

  rownames(k.err) <- c("err", "err.rand")
  colnames(k.err) <- k.range

  err.del <- sapply(1:(ncol(k.err) - 1), function(i) k.err[,i] - k.err[,i + 1])
  colnames(err.del) <- colnames(k.err)[2:length(colnames(k.err))]
  err.del.diff <- err.del[1,] - err.del[2,]

  if (sum(err.del.diff < 0) > 0) {
    min.idx <- min(which(err.del.diff < 0))
  } else {
    min.idx <- which.min(err.del.diff) + 1
  }

  res <- list()
  res$err <- err.del.diff
  res$k <- k.range[[min.idx]]

  if (do.plot) {
    print(PlotFactorSelection(err.del.diff, font.size = 10))
  }

  return(res)
}

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Feature selection: CLUSTERING OF THE BASIS MATRIX (Rowwise, HCLUST)
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
perform_row_hclust = function(matrix, logitscale) {

  k.row = apply(matrix, 1, function(x) {
    if (logitscale == TRUE) {
      x = (x- min(x)) /(max(x)-min(x))
      norm = 1 / (1 + exp(-6 * (x-0.5)))
      # q <- as.numeric(stats::quantile(x, .95))
      # norm <- 2 / (1 + exp((-2) * x / q)) - 1
    }
    else {
      norm = (x- min(x)) /(max(x)-min(x))
      # norm = x / sum(x)
    }
    k.cl = cutree(hclust(dist(norm)), k = 2)

    centroid.1 = mean(norm[k.cl == 1])
    centroid.2 = mean(norm[k.cl == 2])
    factor.mean = mean(norm)
    bss = (2*(factor.mean - centroid.1)^2) + (2*(factor.mean - centroid.2)^2)

    sse.1 = sum((norm[k.cl == 1] - mean(norm[k.cl == 1]))^2)
    sse.2 = sum((norm[k.cl == 2] - mean(norm[k.cl == 2]))^2)

    # centroids = c(centroid.1, centroid.2)
    # max.c = which.max(centroids)
    # if(mean(centroids[-1*max.c]) == 0) {print("UNDERFLOW :-(")}
    # logfc = log2(centroids[max.c]/centroids[-1*max.c])

    d = dist(as.data.frame(norm))
    sil = silhouette(k.cl, d)[, 3]
    count.zero = sil == 0
    count.zero = length(count.zero[count.zero == TRUE])
    if (count.zero == 1) {
      sil.min = min(sil[sil != 0])
      sil.mean = mean(sil[sil != 0])
      ftr = "yes"
    } else {
      sil.min = min(sil)
      sil.mean = mean(sil)
      ftr = "no"
    }

    cl.deltaCenter = mean(x[k.cl == 1]) - mean(x[k.cl == 2])
    cl.deltaCenter.norm = mean(norm[k.cl == 1]) - mean(norm[k.cl == 2])

    return(list(attribution = k.cl,
                sil_mean = sil.mean,
                sil_min = sil.min,
                cl_deltaCenter_norm = cl.deltaCenter.norm,
                cl_deltaCenter = cl.deltaCenter,
                centroid_1 = centroid.1,
                centroid_2 = centroid.2,
                sse_1 = sse.1,
                sse_2 = sse.2,
                ftr = ftr,
                bss = bss))
  })
  return(k.row)
}

compute_feature_stats_hclust = function(obj, rank, ftrs = "all", rank_filter = c(), logitscale = FALSE) {
  nmf.obj = obj$estim.r
  W = basis(nmf.obj$fit[[rank]])
  if (length(rank_filter) != 0) {
    W = W[, rank_filter]
  }
  if (ftrs != "all") {
    W = W[ftrs, ]
  }

  k.row = perform_row_hclust(matrix = W, logitscale)

  k.silmean = unlist(lapply(k.row, function(r) r$sil_mean))
  k.silmin = unlist(lapply(k.row, function(r) r$sil_min))

  k.ftr = unlist(lapply(k.row, function(r) r$ftr))

  k.centroid.1 = unlist(lapply(k.row, function(r) r$centroid_1))
  k.centroid.2 = unlist(lapply(k.row, function(r) r$centroid_2))

  k.sse.1 = unlist(lapply(k.row, function(r) r$sse_1))
  k.sse.2 = unlist(lapply(k.row, function(r) r$sse_2))

  k.bss = unlist(lapply(k.row, function(r) r$bss))
  # Delta between cluster means
  k.deltaMean = unlist(lapply(k.row, function(r) r$cl_deltaCenter))
  # Delta between cluster means (transformed data!)
  k.deltaMean.norm = unlist(lapply(k.row, function(r) r$cl_deltaCenter_norm))
  # Extract Signature combinations generated by hclust
  k.attribution = lapply(k.row, function(r) abs(r$attribution))
  # if (length(rank_filter) != 0) {
  #   k.attribution = lapply(k.attribution, function(r) r[rank_filter])
  # }
  k.attribution = do.call(rbind, k.attribution)
  k.ids <- apply(k.attribution, 1, function(r) paste(r, collapse = ""))
  # Set FeatureStats
  feature.stats = DataFrame("cluster" = k.ids,
                            "delta_mean" = k.deltaMean,
                            "delta_mean_norm" = k.deltaMean.norm,
                            "mean_sil" = k.silmean,
                            "min_sil" = k.silmin,
                            "centroid_1" = k.centroid.1,
                            "centroid_2" = k.centroid.2,
                            "sse_1" = k.sse.1,
                            "sse_2" = k.sse.2,
                            "ftr" = k.ftr,
                            "bss" = k.bss)

  # Re-write cluster ids in a more useful binary code
  i = which(feature.stats$delta_mean_norm > 0 & grepl(feature.stats$cluster, pattern = "2"))
  feature.stats$cluster[i] = gsub("2", "0", feature.stats$cluster[i])
  i = which(feature.stats$delta_mean_norm < 0 & grepl(feature.stats$cluster, pattern = "2"))
  feature.stats$cluster[i] <- gsub("2", "1", gsub("1", "0", feature.stats$cluster[i]))
  feature.stats$GENE_SYMBOL = rownames(feature.stats)

  return(feature.stats)
}

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
anno_w = function(.obj, .rank, .rank_filter = NULL) {

  mat.w = basis(.obj$estim.r$fit[[.rank]])
  if(!is.null(.rank_filter)) {
    mat.w =  mat.w[, .rank_filter]
  }
  for (i in 1:nrow(mat.w)) {
    mat.w[i, ] = (mat.w[i, ]- min(mat.w[i, ])) /(max(mat.w[i, ])-min(mat.w[i, ]))
  }
  # mat.w = sweep(mat.w, 1L, rowSums(mat.w, na.rm = T), '/', check.margin = FALSE)
  mat.w = mat.w > 0.5

  mat.w = as.data.frame(mat.w * 1)
  m = apply(mat.w, 1, function(x){
    paste0(x, collapse = "")
  })
  m = data.frame(cluster = m)
  m$GENE_SYMBOL = rownames(m)
  # print(m["CD4", ])
  # print(m["IL2", ])
  # print(m["IL2RA", ])
  # print(m["TNFRSF9", ])
  m
}

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Feature selection: features from nmf (method based on nmf package)
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
basis_marker = function(
  object,
  colnorm = FALSE,
  qu = FALSE,
  lfc = FALSE,
  eps = 1e-20) {

  object = object + eps
  s = apply(object, 1, function(g){
    g = abs(g)
    p_i = g/sum(g)
    crossprod(p_i, log2(p_i))
    # g = object["IL2", ]
    #   g <- abs(c(0.99, .01, .1, .01))
    #   g
    #   p_i <- g/sum(g)
    #   # p_i = c(.97,.01,.01,.01)
    #   print(p_i)
    #   print(log2(p_i + eps))
    #   print(crossprod(p_i, log2(p_i)))
    #   sum(p_i * log2(p_i))
    #   print(1 + crossprod(p_i, log2(p_i)) / log2(ncol(object)))
    #   print("-----------")
  })
  s = 1 + s / log2(ncol(object))

  if (qu == TRUE) {th = quantile(s)[4]} else{th = median(s)}
  print(summary(s))
  print(paste0("NAs: ", s[is.na(s)]))
  sel1 = s >= th
  print(table(sel1))

  if (colnorm == TRUE) {
    basis.norm = colnorm(object)
  } else {
    basis.norm = object
  }

  if (lfc == TRUE) {
    gene.factor.logfc = log2(apply(object, 1, function(x) {
      max.i <- which.max(x)
      if(mean(x[-1*max.i]) == 0) {print("NOOOOOOO")}
      (x[[max.i]])/(mean(x[-1*max.i]))
    }))
    sel2 = gene.factor.logfc >= 2
    print(table(sel2))

    cut = sel1 & sel2
    print(table(cut))
    basis.filterd = object[cut, ]
  } else {
    cut = sel1
    print(table(cut))
    basis.filterd = object[cut, ]
  }

  pIndx = unlist(apply(basis.filterd,1,which.max))
  marker.list <- list()
  for(i in sort(unique(pIndx))) {
    marker.list[[i]] = names(pIndx[pIndx==i])
  }

  marker.list

}

## >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
## # Feature selection: Cogaps - PatterMarkers
## https://academic.oup.com/bioinformatics/article/33/12/1892/2975325
## https://doi.org/10.1093/bioinformatics/btx058
## Orginal function: patternMarkers()
## https://cran.r-project.org/web/packages/NMF/NMF.pdf
## >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

## !!! AENDERUNG DER LISTE "gBYp". Die Authoren machen hier folgendes:
## Hat ein Gen (z.B. RP11-143J12.2) ein Punkit mit Suffix hinten dran, wird
## nur der Prefix behalten. Da die finale Tabelle "ssgenes.th" nach der Liste
## "gBYp" gefiltert wird fallen alle Gene mit erwaehnter Suffix raus (da "ssgenes.th"
## diese Suffix noch hat). DAS HABE ICH GEAENDERT !!!
## Vermutlicher Beweggrund der Authoren: Nur Proteine sollen behalten werden.


cogaps_marker = function(object, threshold="all", lp=NA) {
  nGenes <- nrow(basis(object))
  nPatterns <- ncol(basis(object))

  # find the A with the highest magnitude
  Arowmax <- t(apply(basis(object), 1, function(x) x/max(x)))

  if (!is.na(lp))
  {
    if (length(lp) != nPatterns)
    {
      warning("lp length must equal the number of columns of the Amatrix")
    }
    sstat <- apply(Arowmax, 1, function(x) sqrt(t(x-lp)%*%(x-lp)))
    ssranks <-rank(sstat)
    ssgenes.th <- names(sort(sstat,decreasing=FALSE,na.last=TRUE))
  }
  else
  {
    # determine which genes are most associated with each pattern
    sstat <- matrix(NA, nrow=nGenes, ncol=nPatterns,
                    dimnames=dimnames(basis(object)))
    ssranks <- matrix(NA, nrow=nGenes, ncol=nPatterns,
                      dimnames=dimnames(basis(object)))
    ssgenes <- matrix(NA, nrow=nGenes, ncol=nPatterns, dimnames=NULL)

    i = 1
    for (i in 1:nPatterns)
    {
      lp <- rep(0,nPatterns)
      lp[i] <- 1


      lp
      head(sstat)
      head(ssranks)
      head(Arowmax)


      a = Arowmax[1,] - lp
      a = c(0,.5,.5,.5)
      sum(a)
      sqrt(crossprod(a,a))
      sqrt(t(a)%*%(a))

      sstat[,i] <- unlist(apply(Arowmax, 1, function(x) sqrt(t(x-lp)%*%(x-lp))))
      ssranks[,i]<-rank(sstat[,i])
      ## ssgenes (column = factor with genes by increasing ranks)
      ssgenes[,i]<-names(sort(sstat[,i],decreasing=FALSE,na.last=TRUE))
    }
    if (threshold=="cut")
    {
      # which(ssranks[ssgenes[,x],x] = Genes ordered (increasing) by rank for current factor
      geneThresh <- sapply(1:nPatterns,function(x) min(which(ssranks[ssgenes[,x],x] > apply(ssranks[ssgenes[,x],],1,min))))
      ssgenes.th <- sapply(1:nPatterns,function(x) ssgenes[1:geneThresh[x],x])
    }
    else if (threshold=="all")
    {
      pIndx<-unlist(apply(sstat,1,which.min))
      gBYp <- list()
      for(i in sort(unique(pIndx))) {
        # gBYp[[i]]<-sapply(strsplit(names(pIndx[pIndx==i]),"[.]"),function(x) x[[1]][1])
        gBYp[[i]] = names(pIndx[pIndx==i])
      }
      ssgenes.th <- lapply(1:max(sort(unique(pIndx))), function(x)
        ssgenes[which(ssgenes[,x] %in% gBYp[[x]]),x])
    }
    else
    {
      stop("Threshold arguement not viable option")
    }
  }
  return(list("PatternMarkers"=ssgenes.th, "PatternRanks"=ssranks,
              "PatternMarkerScores"=sstat))
}

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
add_pseudo_mg_4to5 = function(cl) {
  l = list()
  for (i in 1:length(cl)) {
    spl = strsplit(cl[i], "")[[1]]
    r = NULL
    if (spl[3] == 1 & spl[4] == 1) {
      r = paste0(substring(cl[i], 1, 3), 1, substring(cl[i], 4, 4))
    }
    if (spl[3] == 0 & spl[4] == 0) {
      r = paste0(substring(cl[i], 1, 3), 0, substring(cl[i], 4, 4))
    }
    if (spl[3] == 1 & spl[4] == 0) {
      r = paste0(substring(cl[i], 1, 3), 0, substring(cl[i], 4, 4))
    }
    if (spl[3] == 0 & spl[4] == 1) {
      r = paste0(substring(cl[i], 1, 3), 0, substring(cl[i], 4, 4))
    }
    if (is.null(r)) {
      r = cl[i]
    }

    stopifnot(length(r) != 3)

    l[[i]] = r
  }
  l
}

ch_mg_5 = function(cl) {
  l = list()
  for (i in 1:length(cl)) {
    spl = strsplit(cl[i], "")[[1]]
    r = NULL
    if (spl[3] == 1 & spl[5] == 0) {
      r = paste0(substring(cl[i], 1, 3), 0, substring(cl[i], 5, 5))
    }
    if (spl[3] == 0 & spl[5] == 1) {
      r = paste0(substring(cl[i], 1, 3), 0, substring(cl[i], 5, 5))
    }
    if (is.null(r)) {
      r = cl[i]
    }

    l[[i]] = r
  }
  l
}