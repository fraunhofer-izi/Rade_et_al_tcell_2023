# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Libraries & Functions
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
.cran_packages = c("devtools","reshape2", "ggrepel", "stringr", "dplyr", "scales", "data.table",
                   "yaml", "cowplot", "scico", "ggtext", "circlize", "magick")
.bioc_packages = c("SummarizedExperiment", "ComplexHeatmap", "edgeR")

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

library(clusterProfiler)

source("code/R/eda-plots.R")
source("paper-helper.R")
theme_set(mytheme())

colors_dutch <- c(
  '#FFC312','#C4E538','#12CBC4','#FDA7DF','#ED4C67',
  '#F79F1F','#A3CB38','#1289A7','#D980FA','#B53471',
  '#EE5A24','#009432','#0652DD','#9980FA','#833471',
  '#EA2027','#006266','#1B1464','#5758BB','#6F1E51'
)

colors_spanish <- c(
  '#40407a','#706fd3','#f7f1e3','#34ace0','#33d9b2',
  '#2c2c54','#474787','#aaa69d','#227093','#218c74',
  '#ff5252','#ff793f','#d1ccc0','#ffb142','#ffda79',
  '#b33939','#cd6133','#84817a','#cc8e35','#ccae62'
)

tol12qualitative=c(
  "#332288", "#6699CC", "#88CCEE", "#44AA99", "#117733",
  "#999933", "#DDCC77", "#661100", "#CC6677", "#EE99AA",
  "#997700", "#AA4499"
)

custom_colors = c(colors_dutch, colors_spanish)


# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
manifest = yaml.load_file("assets/manifest.yaml")
work = manifest$workdata

se.izi = readRDS(paste0(work, manifest$limma_izi_thp_th0_act))$se

dge.act = readRDS(paste0(work, manifest$limma_izi_thp_th0_act))
dge.ctr = readRDS(paste0(work, manifest$limma_izi_ctrl))

dge.act.res = dge.act$dge.res
dge.ctr.res = dge.ctr$dge.res
rem = names(dge.ctr.res)[!grepl("^ACT_vs_CTR", names(dge.ctr.res))]
dge.ctr.res = dge.ctr.res[rem]

gc.anno = readRDS(paste0(work, manifest$gencode_v29_features))

font.size = 7
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Heatmap
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
plot_ht = function(.obj, gs = NULL) {

  mat = z_score(assays(.obj)$cpm[gs, ])

  hm.features = c("GROUP_BASE", "DONOR", "HOURS")
  anno.col = as.data.frame(colData(.obj)[,hm.features, drop = F])
  anno.col.colors = degColors(anno.col, cat_values = c("#555555", "#DDDDDD"))

  anno.col = anno.col %>% mutate(
    GROUP_BASE = case_when(
      GROUP_BASE == "naive" ~ "Unactivated (0h)",
      GROUP_BASE == "control" ~ "Negative control",
      GROUP_BASE == "activated" ~ "Activated"
      # TRUE ~ as.character(GROUP_BASE)
    )
  )
  anno.col$GROUP_BASE = factor(anno.col$GROUP_BASE, levels = c("Unactivated (0h)", "Negative control", "Activated"))
  colnames(anno.col) = c("Condition", "Donor", "Time point")
  anno.col = droplevels(anno.col)

  group = c("#000000", "#959595", "#009988")
  names(group) = c("Unactivated (0h)", "Negative control", "Activated")

  time.point = c("#BBBBBB", "#999933", "#DDCC77", "#CC6677", "#882255", "#AA4499")
  names(time.point) = c("0h", "6h", "12h", "24h", "48h", "72h")

  Donor = c("#4477AA", "#66CCEE", "#228833", "#CCBB44")
  names(Donor) = c("D3", "D4", "D5", "D6")


  anno.col.colors = list(
    "Condition" = group,
    "Donor" = Donor,
    "Time point" = time.point
  )

  anno.col.colors$Condition = anno.col.colors$Condition[names(anno.col.colors$Condition) %in% anno.col$Condition]

  ht.fonsize = 8
  ht.anno.col = HeatmapAnnotation(
    df = anno.col,
    col = anno.col.colors,
    gap = unit(rep(1,length(anno.col)), "mm"),
    annotation_legend_param = list(title_gp = gpar(fontsize = ht.fonsize),
                                   # title=c("a", "v"),
                                   labels_gp = gpar(fontsize = ht.fonsize),
                                   grid_height = unit(3, "mm"), grid_width = unit(3, "mm"),
                                   title = NULL,
                                   ncol = 1),
    annotation_height = unit(rep(2,length(anno.col)), "mm"),
    simple_anno_size_adjust = TRUE,
    annotation_name_side = "right",
    show_annotation_name = T,
    annotation_name_gp = gpar(fontsize = ht.fonsize),
    na_col = "black"
  )

  paletteLength = 50
  # use floor and ceiling to deal with even/odd length pallettelengths
  myBreaks = c(
    seq(min(mat), 0, length.out = ceiling(paletteLength/2)),
    seq(max(mat)/paletteLength, max(mat), length.out=floor(paletteLength/2))
  )

  range.mid = subset(myBreaks, myBreaks <= 2 & myBreaks >= -2)
  colpal.mid = scico(length(range.mid), palette = 'vik', begin = .2, end = .8)
  bord = (paletteLength - length(range.mid)) / 2
  colpal = c(rep(colpal.mid[1], bord), colpal.mid, rep(colpal.mid[length(colpal.mid)], bord))

  # for complexHeatmap use always colorRamp2
  # https://jokergoo.github.io/ComplexHeatmap-reference/book/a-single-heatmap.html
  col_fun = colorRamp2(myBreaks[1:length(colpal)], colpal)

  ht1 =
    Heatmap(
      mat, name = "mat",
      use_raster = TRUE, raster_by_magick = T,
      show_row_names = FALSE,
      show_column_names = FALSE,
      column_names_gp = gpar(fontsize = ht.fonsize),
      col = col_fun,
      column_dend_gp = gpar(lwd = .5),
      column_dend_height = unit(5, "mm"),
      row_dend_gp = gpar(lwd = .5),
      row_dend_width = unit(5, "mm"),
      clustering_method_columns = "ward.D",
      clustering_method_rows = "ward.D",
      heatmap_legend_param = list(color_bar = "continuous", ncol = 1, direction = "vertical",
                                  grid_width  = unit(.2, "cm"),
                                  legend_height = unit(1, "cm"),
                                  title = "Z-scores",
                                  title_position = "topleft",
                                  labels_gp = gpar(fontsize = ht.fonsize),
                                  title_gp = gpar(fontsize = ht.fonsize)),
      row_dend_reorder = T, column_dend_reorder = F,
      top_annotation = ht.anno.col, border = TRUE, border_gp = gpar(lwd = 1)
    )
  ht1
  # hm = draw(ht1, merge_legend = TRUE)
}

dge.res.act = data.frame(rbindlist(dge.act$dge.res, idcol = "CONTRAST"))
dge.res.act = dge.res.act %>% .[!is.na(dge.res.act$confect), ]
dge.res.act = dge.res.act[abs(dge.res.act$confect) > log2(1.2) , ]
dge.res.act.ftrs = unique(dge.res.act$ENSEMBL_ID_ABBR)

f = names(dge.ctr$dge.res)[!grepl("^ACT_vs_CTR", names(dge.ctr$dge.res))]
dge.res.ctrl = data.frame(rbindlist(dge.ctr$dge.res[f], idcol = "CONTRAST"))
dge.res.ctrl = dge.res.ctrl %>% .[!is.na(dge.res.ctrl$confect), ]
dge.res.ctrl = dge.res.ctrl[abs(dge.res.ctrl$confect) > log2(1.2) , ]
dge.res.ctrl.ftrs = unique(dge.res.ctrl$ENSEMBL_ID_ABBR)

.obj = se.izi
.obj = .obj[, .obj$GROUP_BASE != "activated"]
colData(.obj) = droplevels(colData(.obj))
# ht.ctrl.sign = plot_ht(.obj, gs = dge.res.ctrl.ftrs[1:100])
ht.ctrl.sign = plot_ht(.obj, dge.res.ctrl.ftrs)

.obj = se.izi
.obj = .obj[, .obj$GROUP_BASE != "control"]
colData(.obj) = droplevels(colData(.obj))
# ht.act.sign = plot_ht(.obj, dge.res.act.ftrs[1:100])
ht.act.sign = plot_ht(.obj, dge.res.act.ftrs)

.obj = se.izi
# ht.int.sign = plot_ht(.obj, intersect(dge.res.act.ftrs, dge.res.ctrl.ftrs)[1:100])
ht.int.sign = plot_ht(.obj, intersect(dge.res.act.ftrs, dge.res.ctrl.ftrs))

# If an error message appears, it may be due to imagemagick
# https://github.com/ImageMagick/ImageMagick6/issues/56
# User: yoursunny
ggsave2(filename="../figures/additional/fig_5_heatmap.png",
        plot = plot_grid(
          grid.grabExpr(draw(ht.act.sign, merge_legend = T)),
          grid.grabExpr(draw(ht.ctrl.sign, merge_legend = T)),
          grid.grabExpr(draw(ht.int.sign, merge_legend = T)),
          ncol = 2,
          labels = c("A", "B", "C"), label_fontface = "bold", label_size = 10
        ),
        width = 210,
        height = 290,
        dpi = 100,
        bg = "white",
        units = "mm")

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Define Reactome levels
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
reactomePathways = read.delim("enrichment_analysis/reactome_hierarchy/ReactomePathways.txt", header = F, sep = "\t", check.names = F)
colnames(reactomePathways) = c("ID", "DESCRIPTION", "SPECIES")
reactomePathways = reactomePathways[reactomePathways$SPECIES == "Homo sapiens", ]

reactome.rel = read.delim("enrichment_analysis/reactome_hierarchy/ReactomePathwaysRelation.txt", header = F, sep = "\t", check.names = F)
colnames(reactome.rel) = c("R1", "R2")
reactome.rel = reactome.rel[grepl("\\-HSA\\-", reactome.rel$R2), ]

reactome.paths = reactomePathways %>% dplyr::select(ID)

reactome.paths = base::merge(reactome.paths, reactome.rel, by.x = "ID", by.y = "R2", all.x = T)
colnames(reactome.paths)[colnames(reactome.paths) == "R1"] = "H1"

for (i in 1:11) {
  x1 = paste0("H",i)
  x2 = paste0("H",i+1)
  reactome.paths = base::merge(reactome.paths, reactome.rel, by.x = x1, by.y = "R2", all.x = T)
  colnames(reactome.paths)[colnames(reactome.paths) == "R1"] = x2

}

reactome.paths$H12 = NULL
reactome.paths = rev(reactome.paths)
reactome.paths = reactome.paths[!duplicated(reactome.paths$ID), ]

r.l = list()
for (i in 1:nrow(reactome.paths)) {
  p = as.character(reactome.paths[i,])
  p = p[!is.na(p)]

  p.top = p[length(p)]
  if (length(p) == 1) {
    r.l[[p[1]]] = c(p[1], p[1])
  }
  else if (length(p) == 2) {
    r.l[[p[1]]] = c(p[1], p[2])
  }
  else if (length(p) > 2) {
    r.l[[p[1]]] = c(p[length(p) -1], p[length(p)])
  }
  else {
    r.l[[p[1]]] = "bra"
  }
}
r = as.data.frame(do.call("rbind", r.l))
r$TOP_1 = reactomePathways$DESCRIPTION[match(r$V1, reactomePathways$ID)]
r$TOP_2 = reactomePathways$DESCRIPTION[match(r$V2, reactomePathways$ID)]

reactomePathways$TOP_1 = r$TOP_1[match(reactomePathways$ID, rownames(r))]
reactomePathways$TOP_2 = r$TOP_2[match(reactomePathways$ID, rownames(r))]

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# ORA
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
ora_plot = function(ora.obj, .run = "go") {

  res.cP.cc.go.simp = ora.obj
  cp.res = res.cP.cc.go.simp@compareClusterResult

  cp.res = cp.res %>%
    tidyr::separate(GeneRatio, c('num_query_hits','num_query_all'), sep='/',remove = F, convert = T) %>%
    tidyr::separate(BgRatio, c('num_background_hits','num_background_all'), sep='/',remove = F, convert = T) %>%
    dplyr::mutate(
      enrichment_factor =
        round(as.numeric((num_query_hits/num_query_all) /
                           (num_background_hits/num_background_all)),
              digits = 1))

  cp.res = mutate(cp.res, richFactor = Count /num_background_hits)
  cp.res = cp.res[cp.res$Count >= 5, ]
  cp.res = cp.res[cp.res$p.adjust < 0.01, ]
  cp.res$geneID = NULL

  cp.res.l = split(cp.res, cp.res$Cluster)

  tops.l = lapply(cp.res.l, function(x){
    x = x[order(x$richFactor, decreasing = T), ]
    # x = x[order(x$p.adjust, decreasing = T), ]
    x = head(x, 20)
    x
  })
  tops.l = base::Filter(function(x) dim(x)[1] > 0, tops.l)

  tops = as.data.frame(data.table::rbindlist(tops.l))
  res.df = subset(cp.res, ID %in% tops$ID)

  lvls = naturalsort(names(tops.l))
  l = list()
  for (i in lvls) { l[[i]] = subset(res.df, Cluster == i) }
  res.df = do.call("rbind", l)
  res.df$Description = factor(res.df$Description, levels = rev(unique(res.df$Description)))

  t.lvls = c("6h", "12h", "24h", "48h", "72h")
  res.df$CTRST = factor(res.df$CTRST, t.lvls)

  facet.lbls = c(
    "Act.\n(up)", "Act.\n(down)",
    "Neg.ctrl.\n(up)", "Neg.ctrl\n(down)"
    # "Act. and Neg.ctrl.\n(diff)"
  )
  names(facet.lbls) = c(
    "ACT_UP", "ACT_DWN",
    "CTRL_ALL_UP", "CTRL_ALL_DWN"
    # "INT_DIFF"
  )

  res.df$GROUP = factor(res.df$GROUP, levels = names(facet.lbls))

  l.tmp = list()
  for (grp in names(facet.lbls)[!names(facet.lbls) %in% unique(res.df$GROUP)]) {
    tmp.1 = data.frame(
      Cluster = paste0(t.lvls, ".", grp),
      CTRST = t.lvls,
      GROUP = grp
    )
    tmp.2 = res.df[1:5, 4:ncol(res.df)]
    tmp.2[, 1:ncol(tmp.2)] = NA

    l.tmp[[grp]] = cbind(tmp.1, tmp.2)
  }
  l.tmp = do.call("rbind", l.tmp)

  res.df = rbind(res.df, l.tmp)

  if (.run == "go") {

    pl =
      ggplot(res.df, aes(x = CTRST, y = Description)) +
      geom_point(aes(fill = p.adjust, size = richFactor), pch=21, stroke=.3) +
      scale_size(range = c(1,4)) +
      mytheme_grid(base_size = font.size) +
      theme(
        panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank(),
        axis.ticks = element_blank()
      ) +
      facet_wrap( ~ GROUP, nrow = 1, labeller = as_labeller(facet.lbls)) +
      scale_fill_scico(palette = 'bilbao', begin = 0.1, end = .8, direction = -1) +
      # scale_colour_gradient(high = "#332288", low = "#DDCC77", name = "Adjusted\np-value") +
      guides(
        fill = guide_colorbar(
          title.position = 'top',
          title = "Adjusted\np-value",
          title.vjust = 1, barwidth = unit(.4, 'lines'), barheight = unit(5, 'lines')
        )
      ) +
      guides(size = guide_legend(ncol = 1, override.aes =, title.hjust = .5, title = "Rich factor", order = 1)) +
      xlab("Contrast") +
      ylab(NULL)
  }

  if (.run == "reactome") {

    res.df$TOP  = reactomePathways$TOP_2[match(res.df$ID, reactomePathways$ID)]
    top.level = unique(res.df$TOP)
    top.level.col = c(tol12qualitative, "grey", "#555555", "black", "orange", "green", "red")[1:length(top.level)]
    # top.level.col = c(colors_spanish)[1:length(top.level)]
    names(top.level.col) = top.level

   pl =
      ggplot(res.df, aes(x = CTRST, y = Description, size = richFactor, colour = TOP)) +
      geom_point() +
      scale_size(range = c(1,4)) +
     scale_colour_discrete(na.translate = F, type = top.level.col) +
      mytheme_grid(base_size = font.size) +
      theme(
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.ticks = element_blank()
      ) +
      facet_wrap( ~ GROUP, nrow = 1, labeller = as_labeller(facet.lbls)) +
      guides(colour = guide_legend(ncol = 1, override.aes = list(size=4), title.hjust = 0, title =  "Top level of the\nReactome pathway hierarchy")) +
      guides(size = guide_legend(ncol = 1, override.aes =, title.hjust = .5, title = "Rich factor", order = 1)) +
      xlab("Contrast") +
      ylab(NULL)
  }
  pl
}

ora.obj.go = readRDS("enrichment_analysis/figure_5_supps_ora_go_act_neg_ctrl_comb.rds")
ora.obj.reactome = readRDS("enrichment_analysis/figure_5_supps_ora_reactome_act_neg_ctrl_comb.rds")

ora_plot(ora.obj = ora.obj.go, .run = "go")
ora_plot(ora.obj = ora.obj.reactome, .run = "reactome")


ggsave2(filename="../figures/additional/fig_5_act_neg_ctrl_ora_go.pdf",
        plot = ora_plot(ora.obj = ora.obj.go, .run = "go"),
        width = 310,
        height = 390,
        dpi = 100,
        bg = "white",
        units = "mm")

ggsave2(filename="../figures/additional/fig_5_act_neg_ctrl_ora_reactome.pdf",
        plot = ora_plot(ora.obj = ora.obj.reactome, .run = "reactome"),
        width = 360,
        height = 390,
        dpi = 100,
        bg = "white",
        units = "mm")



