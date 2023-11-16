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

nmf_stats = function(obj, anno, hm.fonsize = 5.5, main.title = "Title", .r = 99) {
  library("ComplexHeatmap")
  estim.r = obj$estim.r
  se = obj$se

  if(!is.null(obj$estim.r.random)) {
    estim.r.random = obj$estim.r.random
  }

  h.grid.l = list()

  # i = "2"
  for (i in names(estim.r$fit)) {
    consensus.m = estim.r$fit[[as.character(i)]]@consensus
    # distance and cluster
    dend = hclust(as.dist(1-consensus.m), method = "average")

    # top annotation bar
    anno.top = colData(se)[, anno, drop = F]
    anno.top.col = degColors(anno.top, tol21rainbow = tol10qualitative)

    tmp = anno.top.col$Time_point
    tmp = time.col[match(names(tmp), names(time.col))]
    anno.top.col$Time_point = tmp

    top.anno = HeatmapAnnotation(
      df = anno.top,
      col = anno.top.col,
      show_legend = F,
      simple_anno_size = unit(.2, "cm"),
      annotation_name_gp = gpar(fontsize = hm.fonsize - 2))

    lvls = naturalsort(unique(anno.top$Time_point))
    top.anno@anno_list$Time_point@color_mapping@levels = lvls

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
                # col = c("#f2f2f7", scico(30, palette = 'acton', direction = -1)),
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

      anno.leg$Time_point@grob$children[[2]]$children[[2]]$gp$fill = anno.top.col$Time_point

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

  nmf.stats$using = nmf.stats$rank == .r

  coph.p =
    ggplot(nmf.stats, aes(x = rank, y = cophenetic)) +
    geom_line() +
    geom_point(aes(colour = using)) +
    scale_colour_manual(values = c("black", "red")) +
    mytheme_grid(base_size = hm.fonsize) +
    theme(aspect.ratio = 1, axis.text.y = element_text(angle = 90, hjust = .5), legend.position = "none") +
    ylab("Cophenetic coefficient")

  silh = ggplot(nmf.stats, aes(x = rank, y = silhouette.consensus)) +
    geom_line() +
    geom_point(aes(colour = using)) +
    scale_colour_manual(values = c("black", "red")) +
    mytheme_grid(base_size = hm.fonsize) +
    theme(aspect.ratio = 1,axis.text.y = element_text(angle = 90, hjust = .5), legend.position = "none") +
    ylab("Avg silhouette width")

  stats.grid = plot_grid(plot_grid(coph.p, NULL, silh, rel_widths = c(1, .03, 1), nrow = 1), NULL, rel_widths = c(.8, .2))

  if(!is.null(obj$estim.r.random)) {
    nmf.stats.r = estim.r.random$measures %>%
      dplyr::select(rank, cophenetic, dispersion, silhouette.consensus) %>%
      dplyr::mutate(run = "randomized")
    nmf.stats.r$using = FALSE

    coph.p = ggplot(rbind(nmf.stats, nmf.stats.r), aes(x = rank, y = cophenetic)) +
      geom_line() +
      geom_point(aes(colour = using)) +
      scale_colour_manual(values = c("black", "red")) +
      mytheme_grid(base_size = hm.fonsize) +
      theme(panel.spacing = unit(1, "lines"), legend.position = "none") +
      # theme(aspect.ratio = 1, axis.text.y = element_text(angle = 90, hjust = .5)) +
      ylab("Cophenetic coefficient") +
      ylim(c(NA, 1)) +
      facet_wrap( ~ run, scales = "free_y")

    silh = ggplot(rbind(nmf.stats, nmf.stats.r), aes(x = rank, y = silhouette.consensus)) +
      geom_line() +
      geom_point(aes(colour = using)) +
      scale_colour_manual(values = c("black", "red")) +
      mytheme_grid(base_size = hm.fonsize) +
      theme(panel.spacing = unit(1, "lines"), legend.position = "none") +
      ylab("Avg silhouette width") +
      ylim(c(NA, 1)) +
      facet_wrap( ~ run, scales = "free_y")

    stats.grid = plot_grid(plot_grid(
      coph.p, NULL, silh, rel_widths = c(1, .1, 1), nrow = 1, labels = c("B", "", "C"), label_size = 10, vjust = 0
    ), NULL, rel_widths = c(.9, .05))
  }

  # now add the title
  title = ggdraw() +
    draw_label(
      main.title,
      fontface = 'bold',
      x = .32,
      hjust = 0
    )

  res = plot_grid(
    title, consensus.maps.all, NULL, stats.grid, nrow = 4,
    rel_heights = c(.04, .8, .03, .15), align = "vh",
    labels = c("", "A", "", ""), label_size = 10
  )
  return(res)
}
