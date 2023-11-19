mg.labels = c(
  "M1", "M2", "M3", "M5",
  "M1/M2", "M1/M3", "M1/M5",
  "M2/M3", "M2/M5", "M3/M5",
  "M1/M2/M3", "M2/M3/M5",
  "M4"
)

cl.code = c(
  "1000", "0100", "0010", "0001",
  "1100", "1010", "1001",
  "0110", "0101", "0011",
  "1110", "0111",
  "00010"
)

names(mg.labels) = cl.code

mg.col = c(
  "#000000", "#004488", "#DDAA33", "#BB5566",
  "white", "white", "white",
  "white", "white", "white",
  "white", "white",
  "white")

names(mg.col) = cl.code

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# LIBRARIES
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
.cran_packages = c("ggplot2","devtools", "yaml", "naturalsort", "reshape2",
                   "patchwork", "tidyr", "ggrepel", "cowplot", "ggpubr",
                   "broom", "dplyr", "rcartocolor", "scico")
.bioc_packages = c("SummarizedExperiment")

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

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# THEME: GGPLOT
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
mytheme = function(base_size = 8, base_family = "") {
  half_line <- base_size/2
  theme_light(base_size = base_size, base_family = base_family) +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_rect(fill = "transparent",colour = NA),
      plot.background = element_rect(fill = "transparent",colour = NA),
      axis.ticks.length = unit(half_line / 2.2, "pt"),
      axis.ticks = element_line(size = .2, color="black"),
      strip.background = element_rect(fill = NA, colour = NA),
      strip.text.x = element_text(size = rel(1), colour = "black"),
      strip.text.y = element_text(size = rel(1), colour = "black"),
      strip.text = element_text(size = rel(1), colour = "black"),
      axis.text = element_text(size = rel(1), colour = "black"),
      axis.title = element_text(size = rel(1), colour = "black"),
      legend.title = element_text(colour = "black", size = rel(1)),
      panel.border = element_rect(fill = NA, colour = "black", size = .25),
      legend.key.size = unit(1, "lines"),
      legend.text = element_text(size = rel(1), colour = "black"),
      legend.key = element_rect(colour = NA, fill = NA),
      legend.background = element_rect(colour = NA, fill = NA),
      plot.title = element_text(hjust = 0, face = "plain", colour = "black", size = rel(1)),
      plot.subtitle = element_text(colour = "black", size = rel(.85))
    )
}

mytheme_grid = function(base_size = 8, base_family = "") {
  half_line <- base_size/2
  theme_light(base_size = base_size, base_family = base_family) +
    theme(
      panel.background = element_rect(fill = "transparent",colour = NA),
      plot.background = element_rect(fill = "transparent",colour = NA),
      axis.ticks.length = unit(half_line / 2.2, "pt"),
      axis.ticks = element_line(size = .2, color="black"),
      strip.background = element_rect(fill = NA, colour = NA),
      strip.text.x = element_text(size = rel(1), colour = "black"),
      strip.text.y = element_text(size = rel(1), colour = "black"),
      strip.text = element_text(size = rel(1), colour = "black"),
      axis.text = element_text(size = rel(1), colour = "black"),
      axis.title = element_text(size = rel(1), colour = "black"),
      legend.title = element_text(colour = "black", size = rel(1)),
      panel.border = element_rect(fill = NA, colour = "black", size = .25),
      legend.key.size = unit(1, "lines"),
      legend.text = element_text(size = rel(1), colour = "black"),
      legend.key = element_rect(colour = NA, fill = NA),
      legend.background = element_rect(colour = NA, fill = NA),
      plot.title = element_text(hjust = 0, face = "plain", colour = "black", size = rel(1)),
      plot.subtitle = element_text(colour = "black", size = rel(.85))
    )
}

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# https://tradeblotter.wordpress.com/2013/02/28/the-paul-tol-21-color-salute/
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
tol2qualitative=c("#4477AA", "#CC6677")
tol3qualitative=c("#4477AA", "#DDCC77", "#CC6677")
tol4qualitative=c("#4477AA", "#117733", "#DDCC77", "#CC6677")
tol5qualitative=c("#332288", "#88CCEE", "#117733", "#DDCC77", "#CC6677")
tol6qualitative=c("#332288", "#88CCEE", "#117733", "#DDCC77", "#CC6677","#AA4499")
tol7qualitative=c("#332288", "#88CCEE", "#44AA99", "#117733", "#DDCC77", "#CC6677","#AA4499")
tol8qualitative=c("#332288", "#88CCEE", "#44AA99", "#117733", "#999933", "#DDCC77", "#CC6677","#AA4499")
tol9qualitative=c("#332288", "#88CCEE", "#44AA99", "#117733", "#999933", "#DDCC77", "#CC6677", "#882255", "#AA4499")
tol10qualitative=c("#332288", "#88CCEE", "#44AA99", "#117733", "#999933", "#DDCC77", "#661100", "#CC6677", "#882255", "#AA4499")
tol11qualitative=c("#332288", "#6699CC", "#88CCEE", "#44AA99", "#117733", "#999933", "#DDCC77", "#661100", "#CC6677", "#882255", "#AA4499")
tol12qualitative=c("#332288", "#6699CC", "#88CCEE", "#44AA99", "#117733", "#999933", "#DDCC77", "#661100", "#CC6677", "#AA4466", "#882255", "#AA4499")
tol14rainbow=c("#882E72", "#B178A6", "#D6C1DE", "#1965B0", "#5289C7", "#7BAFDE", "#4EB265", "#90C987", "#CAE0AB", "#F7EE55", "#F6C141", "#F1932D", "#E8601C", "#DC050C")
tol15rainbow=c("#114477", "#4477AA", "#77AADD", "#117755", "#44AA88", "#99CCBB", "#777711", "#AAAA44", "#DDDD77", "#771111", "#AA4444", "#DD7777", "#771144", "#AA4477", "#DD77AA")
tol18rainbow=c("#771155", "#AA4488", "#CC99BB", "#114477", "#4477AA", "#77AADD", "#117777", "#44AAAA", "#77CCCC", "#777711", "#AAAA44", "#DDDD77", "#774411", "#AA7744", "#DDAA77", "#771122", "#AA4455", "#DD7788")
# ...and finally, the Paul Tol 21-color salute
tol21rainbow= c("#771155", "#AA4488", "#CC99BB", "#114477", "#4477AA", "#77AADD", "#117777", "#44AAAA", "#77CCCC", "#117744", "#44AA77", "#88CCAA", "#777711", "#AAAA44", "#DDDD77", "#774411", "#AA7744", "#DDAA77", "#771122", "#AA4455", "#DD7788")

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# ggplot stuff
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
align_legend <- function(p, hjust = 0.5)
{
  # extract legend
  g <- cowplot::plot_to_gtable(p)
  grobs <- g$grobs
  legend_index <- which(sapply(grobs, function(x) x$name) == "guide-box")
  legend <- grobs[[legend_index]]

  # extract guides table
  guides_index <- which(sapply(legend$grobs, function(x) x$name) == "layout")

  # there can be multiple guides within one legend box
  for (gi in guides_index) {
    guides <- legend$grobs[[gi]]

    # add extra column for spacing
    # guides$width[5] is the extra spacing from the end of the legend text
    # to the end of the legend title. If we instead distribute it by `hjust:(1-hjust)` on
    # both sides, we get an aligned legend
    spacing <- guides$width[5]
    guides <- gtable::gtable_add_cols(guides, hjust*spacing, 1)
    guides$widths[6] <- (1-hjust)*spacing
    title_index <- guides$layout$name == "title"
    guides$layout$l[title_index] <- 2

    # reconstruct guides and write back
    legend$grobs[[gi]] <- guides
  }

  # reconstruct legend and write back
  g$grobs[[legend_index]] <- legend
  g
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
#
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
z_score = function(z) {
  rowmean = apply(z, 1, mean, na.rm=TRUE)
  rowsd = apply(z, 1, sd, na.rm=TRUE)
  rv = sweep(z, 1, rowmean,"-")
  rv = sweep(rv, 1, rowsd, "/")
  return(rv)
}

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Vote-counting stats
# Modified version from https://github.com/csbl-usp/MetaVolcanoR
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
votecount_mv <- function(diffexp = list(),
                         foldchangecol = "logFC",
                         genenamecol = "GENE_SYMBOL",
                         geneidcol = NULL,
                         foldchange = 0) {

  nstud = length(diffexp)

  # Defining DEGs
  diffexp = lapply(diffexp, function(x) {
    dplyr::mutate(x, deg = ifelse(as.numeric(!!as.name(foldchangecol)) > 0, 1, -1))
  })

  if(is.null(geneidcol)) {
    geneidcol = genenamecol
  }

  # Testing if geneIDs are unique
  gid = vapply(diffexp, function(g) {
    length(unique(g[[geneidcol]])) == nrow(g)
  }, logical(1))

  if(all(gid)) {

    # Subsetting the diffexp inputs (columns)
    diffexp = lapply(diffexp, function(x)
      dplyr::select(x,
                    dplyr::matches(paste(
                      c(foldchangecol,
                        geneidcol, '^deg$'), collapse = '|'
                    ))))


    # merging DEG results
    diffexp = rename_col(diffexp, geneidcol)
    meta_diffexp = Reduce(function(x, y)
      merge(x = x, y = y, by = geneidcol, all = TRUE), diffexp)

    genecol = geneidcol
  } else {
    stop("the geneidcol contains duplicated values")
  }

  # --- Defining new vars for visualization
  meta_diffexp %>%
    dplyr::select(dplyr::matches("deg_")) %>%
    data.matrix -> n_deg

  meta_diffexp[['ndeg']] <- rowSums(n_deg^2, na.rm = TRUE)
  meta_diffexp[['ddeg']] <- rowSums(n_deg, na.rm = TRUE)

  return(meta_diffexp)
}

rename_col = function(diffexp, genecol) {
  ns <- names(diffexp)
  des <- lapply(seq(diffexp), function(nstudy) {
    dex <- diffexp[[nstudy]]
    colnames(dex) <- paste(colnames(dex), nstudy, sep = "_")
    colnames(dex)[grep(genecol, colnames(dex))] <- genecol
    return(dex)
  })
  names(des) <- ns
  return(des)
}

cum_freq_data <- function(meta_diffexp, nstud) {
  data.frame("DEGs" = vapply(0:nstud, function(idx) {
    length(which(meta_diffexp[['ndeg']] >= idx))
  }, numeric(1)),
  "ndatasets" = 0:nstud
  )
}

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Topconfects plot
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
topconfect_plot = function(confects, limits=NULL) {

  tab = confects
  limits = c(NA,NA)

  min_effect <- min(0, tab$est_effect, na.rm=TRUE)
  max_effect <- max(0, tab$est_effect, na.rm=TRUE)

  if (min_effect == max_effect) {
    min_effect <- -1
    max_effect <- 1
  }

  if (is.na(limits[1]) & is.na(limits[2])) {
    max_abs_effect <- max(-min_effect,max_effect)
    limits <- c(-max_abs_effect*1.05, max_abs_effect*1.05)
  } else if (is.na(limits[1])) {
    limits[1] <- min_effect * 1.05
  } else if (is.na(limits[2])) {
    limits[2] <- max_effect * 1.05
  }

  assert_that(is.numeric(limits), length(limits) == 2)

  tab$confect_from <- limits[1]
  tab$confect_to <- limits[2]
  positive <- !is.na(tab$confect) & tab$est_effect > 0
  tab$confect_from[positive] <- tab$confect[positive]
  negative <- !is.na(tab$confect) & tab$est_effect < 0
  tab$confect_to[negative] <- tab$confect[negative]

  tab$ORDER = paste0((seq(nrow(tab))), "_", tab$GROUP)
  tab$ORDER = factor(tab$ORDER, levels = rev(tab$ORDER))

  head(data.frame(tab))

  lock_order  <- function(var, rev = FALSE) {
    var <- forcats::fct_inorder(var)
    if (rev) var <- forcats::fct_rev(var)
    var
  }

  p = ggplot(tab, aes(x = est_effect, y = ORDER)) +
    geom_vline(xintercept = 0, lwd = .1) +
    # annotate("rect", xmin = -lfc, xmax = lfc, ymin = -Inf, ymax = Inf, alpha=1, fill="#DDDDDD") +
    geom_segment(aes(yend = ORDER, x = confect_from, xend = confect_to), size = .2) +
    geom_point(shape=23, fill="black") +
    scale_x_continuous(expand=c(0,0), limits = limits, oob = function(a,b) a) +
    scale_y_discrete(labels = setNames( tab$GENE_SYMBOL,tab$ORDER)) +
    mytheme_grid(base_size = 5) +
    theme(axis.ticks = element_blank(),
          panel.spacing = unit(2, "lines")) +
    # facet_wrap( ~ lock_order(GROUP, rev = FALSE), ncol=7, scales = "free") +
    labs(x = "Combined effect size", y = NULL)
  return(p)
}

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Forest plots for top genes
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
lock_order  <- function(var, rev = TRUE) {
  var <- forcats::fct_inorder(var)
  if (rev) var <- forcats::fct_rev(var)
  var
}

guide_axis_label_trans <- function(label_trans = identity, ...) {
  axis_guide <- guide_axis(...)
  axis_guide$label_trans <- rlang::as_function(label_trans)
  class(axis_guide) <- c("guide_axis_trans", class(axis_guide))
  axis_guide
}

guide_train.guide_axis_trans <- function(x, ...) {
  trained <- NextMethod()
  trained$key$.label <- x$label_trans(trained$key$.label)
  trained
}

forest_plots = function(.meta.data = NULL ,
                        .x.max = NULL,
                        .title = NULL,
                        .i.2 = NULL,
                        .tau.2 = NULL,
                        .re.col = NULL,
                        .range = c(1, 3),
                        .err.lwd = .5,
                        .base.size = 7) {
  .meta.data %>%
    ggplot(aes(x = estimate, y = lock_order(study), shape = type, col = type)) +
    geom_point(aes(size = weight), alpha = 1) +
    scale_shape_manual(values = c(15, 18)) +
    scale_size_continuous(range = .range) +
    scale_colour_manual(values = c("#BBBBBB", .re.col)) +
    geom_errorbarh(aes(xmin = conf.low, xmax = conf.high), height = 0, alpha = 1, lwd = .err.lwd) +
    guides(y.sec = guide_axis_label_trans(~ paste(rev(.meta.data$sec_axis)))) +
    xlim(-.x.max, .x.max) +
    geom_vline(xintercept = 0, linetype = "solid", lwd = .05) +
    geom_vline(xintercept = .meta.data[.meta.data$type == "summary", ]$estimate, linetype = "dotted", lwd = .2) +
    mytheme_grid(.base.size) +
    theme(aspect.ratio = 1,
          axis.text.x = element_text(size = rel(.8), colour = "black"),
          legend.position = "none",
          plot.caption = element_text(hjust = 0, size = rel(1)),
          axis.ticks = element_blank()) +
    xlab("Hedge's standardized\nlog2 fold change") +
    ylab(NULL) +
    # scale_x_continuous(expand = c(0, 0), limits = c(-.x.max, .x.max)) +
    ggtitle(.title)
}

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Top genes for each Metagene by decreased logFC (from Limma)
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
mg_tops = function(.obj, .meta.study, .u.cl, .ftrs.ave, .metafor.res, .ntop, .col, .legend.title, .font.size) {

  top.l = list()
  # cl = "0100"
  for (cl in .u.cl) {

    obj.u.cl = subset(.obj, cluster == cl)
    st = do.call("rbind", .meta.study)
    st$HOURS = gsub("_vs.+", "", gsub("HOURS_", "", rownames(st)))
    st = st[st$GENE_SYMBOL %in% rownames(obj.u.cl), ]

    st = st[!is.na(st$confect), ]

    if (cl == "100" | cl == "1000" | cl == "10000") {
      st = subset(st, confect < 0)
    } else {
      st = subset(st, confect > 0)
    }
    st$confect_abs = abs(st$confect)
    st.tops = st %>%
      dplyr::group_by(GENE_SYMBOL) %>%
      dplyr::filter(rank == min(rank)) %>%
      dplyr::filter(confect_abs == max(confect_abs)) %>%
      dplyr::arrange(desc(confect_abs), .by_group = F) %>%
      data.frame()

    st.tops = head(st.tops, n = .ntop)
    st.tops$CLUSTER = cl
    top.l[[cl]] = st.tops
  }

  df = do.call("rbind", top.l)
  df$CLUSTER = factor(df$CLUSTER, levels = .u.cl)

  # median expresseion
  colnames(.ftrs.ave) = gsub(".+_", "", colnames(.ftrs.ave))
  ftrs.ave.norm = apply(.ftrs.ave[rownames(.obj), ], 2, normalize_vector, method = "rank")

  ftrs.ave.melt = reshape2::melt(as.matrix(ftrs.ave.norm))
  colnames(ftrs.ave.melt) = c("GENE_SYMBOL", "HOURS", "EXPRS")
  ftrs.ave.melt$ID = paste0(ftrs.ave.melt$GENE_SYMBOL, "_", ftrs.ave.melt$HOURS)

  df = df %>% mutate(
    TMP = case_when(
      CLUSTER == .u.cl[1] ~ "0h",
      TRUE ~ HOURS
    ))
  df$ID = paste0(df$GENE_SYMBOL, "_", df$TMP)
  df$ID_META = paste0(df$GENE_SYMBOL, "_", df$HOURS)
  df$AVE_EXPRS = ftrs.ave.melt$EXPRS[match(df$ID, ftrs.ave.melt$ID)]

  df$LABELS = paste0(df$GENE_SYMBOL, " | ", df$HOURS)

  lock_order  <- function(var, rev = TRUE) {
    var <- forcats::fct_inorder(var)
    if (rev) var <- forcats::fct_rev(var)
    var
  }
  my_colors = carto_pal(7, "BurgYl")
  df$COL = mg.col[match(df$CLUSTER, names(mg.col))]
  df$CLUSTER_SYMBOL = names(mg.lbls)[match(df$CLUSTER, mg.lbls)]

  df = df %>% mutate(x_max_x = case_when(confect > 0 ~ confect, TRUE ~ -Inf))
  df = df %>% mutate(x_min_x = case_when(confect < 0 ~ confect, TRUE ~ Inf))

  p = df %>%
    ggplot(aes(x = logFC, y = lock_order(LABELS))) +
    scale_color_gradientn(colours = my_colors, limits = c(0,1), breaks = c(0,1)) +
    geom_errorbarh(aes(xmin = x_min_x, xmax = x_max_x), position=position_nudge(y = 0), height = 0, lwd = .2) +
    geom_point(aes(logFC, col = AVE_EXPRS), position=position_nudge(y = 0)) +
    mytheme(base_size = .font.size) +
    theme(legend.position = "bottom",
          axis.ticks.y =  element_blank()) +
    guides(color = guide_colorbar(title.position = 'left', title.hjust = .5, title.vjust = 1, barwidth = unit(6.5, 'lines'), barheight = unit(.5, 'lines'))) +
    xlab(NULL) + ylab(NULL) + labs(color = .legend.title) +
    facet_grid(CLUSTER_SYMBOL ~ ., scales = "free") +
    geom_vline(aes(xintercept = Inf), color = factor(df$COL), size=3) +
    geom_vline(aes(xintercept = 0), size=.1)
  return(p)
}

## >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
## cnetplot function from Clusterprofiler (modified)
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
      scale_colour_gradient2(name = "Combined effect size", low = "#0077BB", mid = "white", high = "#CC3311", na.value = "grey60")
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

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Sample and Gene embeddngs based on W and H from NMF
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
embedding_pattern = function(H, alpha.exp = 1, proj.method = "sammon", dist.use = "cosine", rescale.coords = "bounded") {

  ix = Matrix::colSums(H, na.rm = TRUE) > 0
  if (!all(ix)) {
    print("At least on sample has zero values across all components. These samples were removed")
    H = H[,ix]
  }

  H.coords = get_factor_coords(H, method = proj.method, distance = dist.use, rescale = rescale.coords)
  H.coords = data.frame(H.coords)
  H.coords$name = rownames(H.coords)
  rownames(H.coords) = NULL;

  sample.coords = get_sample_coords(H, H.coords, alpha = alpha.exp)
  return(list(H.coords = H.coords, sample.coords = data.frame(sample.coords)))
}

cos_sim = function(ix)
{
  A = H[ix[1],]
  B = H[ix[2],]
  A.na = is.na(A)
  B.na = is.na(B)
  A = A[A.na == F & B.na == F]
  B = B[A.na == F & B.na == F]

  return( sum(A*B, na.rm = T)/sqrt(sum(A^2, na.rm = T)*sum(B^2, na.rm = T)) )
}

get_factor_coords = function(H, method, distance, rescale) {
  H.t = t(H)

  if (method == "sammon") {
    if (distance == "IC") {
      H.dist = sqrt(2*(1 - as.matrix(usedist::dist_make(H, distance_fcn = MutualInf))))
    } else if (distance == "cosine") {
      n = nrow(H)
      cmb = expand.grid(i=1:n, j=1:n)
      # ix = cmb[4, ]
      # ix = as.integer(ix)
      # H.dist = 1 - as.dist(matrix(apply(cmb, 1, cos_sim),n,n))
      H.dist = 1 - proxy::simil(H.t, method = "cosine", by_rows = F)
    } else if (distance == "euclidean") {
      H.dist = proxy::dist(H.t, method = "Euclidean", p = 1, by_rows = F)
    } else if (distance == "pearson") {
      H.dist <- sqrt(2*(1 - cor(H.t)))
    }
    H.coords = MASS::sammon(H.dist, k = 2, niter = 500)$points
    rownames(H.coords) = colnames(H.t)

  }
  if (method == "pca") {
    H.coords = prcomp(H)$x[, c(1,2)]
  }

  H.coords = apply(H.coords, 2, normalize_vector, method = rescale)
  colnames(H.coords) <- c("x","y")

  return(H.coords)
}

get_sample_coords <- function(H, H.coords, alpha) {
  sample.coords <- t(sapply(1:ncol(H), function(i) {
    H.non.na = !is.na(H[, i])
    H.i = H[, i][H.non.na]
    pull.sum = sum(H.i^alpha)
    x <- sum(H.coords[, 1][H.non.na] * ( H.i^alpha)) / pull.sum
    y <- sum(H.coords[, 2][H.non.na] * ( H.i^alpha)) / pull.sum
    return(c(x,y))
  }))
  colnames(sample.coords) <- c("x","y")
  rownames(sample.coords) <- colnames(H)
  return(sample.coords)
}

ftr.embedding = function(swne.embedding, feature.assoc, features.embed, alpha.exp = 1, scale.cols = T, overwrite = T) {
  feature.assoc = t(feature.assoc[features.embed,])
  stopifnot(nrow(swne.embedding$H.coords) == nrow(feature.assoc))

  if (scale.cols) {
    feature.assoc <- apply(feature.assoc, 2, normalize_vector, method = "bounded")
  }

  feature.coords <- get_sample_coords(feature.assoc, swne.embedding$H.coords, alpha = alpha.exp)
  feature.coords <- data.frame(feature.coords)
  feature.coords$name <- rownames(feature.coords); rownames(feature.coords) <- NULL;

  if (overwrite || is.null(swne.embedding$feature.coords)) {
    swne.embedding$feature.coords <- feature.coords
  } else {
    swne.embedding$feature.coords <- rbind(swne.embedding$feature.coords, feature.coords)
  }
  return(swne.embedding)
}

normalize_vector <- function(x, method = "scale", n_ranks = 10000) {
  stopifnot(method %in% c("rank", "scale", "bounded"))
  if (method == "scale") {
    x.scale <- (x - mean(x, na.rm = T))/sd(x, na.rm = T)
  } else if (method == "rank") {
    x.scale <- rank(x) / length(x)
  } else if (method == "bounded") {
    x.max <- max(x, na.rm = T)
    x.min <- min(x, na.rm = T)
    x.scale <- (x - x.min)/(x.max - x.min)
  }
  return(x.scale)
}

