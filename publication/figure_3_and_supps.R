# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# LIBRARIES
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
.cran_packages = c("ggplot2","devtools", "yaml", "naturalsort", "reshape2",
                   "patchwork", "tidyr", "ggrepel", "cowplot", "ggpubr",
                   "broom", "dplyr", "rcartocolor", "NMF")
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
# Helper functions
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
source("code/eda/meta-effect//helper.R")
source("paper-helper.R")

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Read data
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
manifest = yaml.load_file("assets/manifest.yaml")
work = manifest$workdata

meta.prep = readRDS(paste0(work, manifest$meta_se))
obj.meta = readRDS(paste0(work, manifest$meta_all))

meta.res = obj.meta$meta.res # metafor results

obj.studies.all = readRDS(paste0(work, manifest$meta_effect_sizes_all))
meta.studies = obj.studies.all$studies
# meta.studies = obj.meta$studies # all DGE genes

gc.anno = readRDS(paste0(work, manifest$gencode_v29_features))

# Median expression for each phenotype and time point (RNA-Seq and Array)
# the cpm assay in SummarizedExperiment is in fact rma intensities for arrays (confusing, use MultiAssayExperiment for future projects)
exprs.cpm = assays(meta.prep)$cpm

colnames(exprs.cpm) = gsub(".+\\.", "", colnames(exprs.cpm))
ftrs.ave = grp_ave_sd(.obj = exprs.cpm, .ave = "mean")$ave
rownames(ftrs.ave) = gc.anno$GENE_SYMBOL_DUPL_MARKED[match(rownames(ftrs.ave), gc.anno$ENSEMBL_ID_ABBR)]

# Metafor results
m.res = lapply(meta.res, function(ctrst) {
  ctrst %>% dplyr::select(GENE_SYMBOL, signcon, ntimes, est_effect, est_pvalue, rank, confect, GROUP)
})
meta.res.df = do.call("rbind", m.res)

# After assignment of genes to Metagenes
ftrs.work = readRDS(paste0(work, manifest$pipeline_discovery))
ftrs.th0.wrk = ftrs.work$ftrs.th0.wrk
ftrs.th17.wrk = ftrs.work$ftrs.th17.wrk
ftrs.itreg.wrk = ftrs.work$ftrs.itreg.wrk
ftrs.th2.wrk = ftrs.work$ftrs.th2.wrk
ftrs.th2.wrk = ftrs.work$ftrs.th2.wrk
ftrs.th1.wrk = ftrs.work$ftrs.th1.wrk
nmf.th0 = ftrs.work$nmf.th0
nmf.th17 = ftrs.work$nmf.th17
nmf.itreg = ftrs.work$nmf.itreg
nmf.th2 = ftrs.work$nmf.th2
nmf.th1 = ftrs.work$nmf.th1

nmf.v1.obj = readRDS(paste0(work, manifest$pipeline_verification_v1))
nmf.v2.obj = readRDS(paste0(work, manifest$pipeline_verification_v2))
ftrs.work.v1 = nmf.v1.obj$ftrs.v1.wrk.all
ftrs.work.v2 = nmf.v2.obj$ftrs.v2.wrk.all

nmf.v1 = nmf.v1.obj$nmf.v1
nmf.v2 = nmf.v2.obj$nmf.v2

rank.th0 = ftrs.work$rank.th0
rank.th17 = ftrs.work$rank.th17
rank.itreg = ftrs.work$rank.itreg
rank.th2 = ftrs.work$rank.th2
rank.th1 = ftrs.work$rank.th1

rank.v1 = nmf.v1.obj$rank.v1
rank.v2 = nmf.v2.obj$rank.v2

metag.col.3 = c("#555555", "#DDAA33", "#BB5566")
names(metag.col.3) = c("M1", "M2", "M3")

metag.col.7 = c("#555555", "#004488", "#DDAA33", "#BB5566", "#117733", "#BEBEBE", "green")
names(metag.col.7) = c("M1", "M2", "M3", "M4", "M5", "M6")

mg.col = c(
  "#000000", "#DDAA33", "#BB5566",
  "#000000", "#004488", "#DDAA33", "#BB5566",
  "#000000", "#004488", "#DDAA33", "#117733", "#BB5566")

names(mg.col) = c(
  "100", "010", "001",
  "1000", "0100", "0010", "0001",
  "10000", "01000", "00100", "00010", "00001")

mg.lbls = names(mg.col)
names(mg.lbls) = c(
  "M1", "M2", "M5",
  "M1", "M2", "M3", "M5",
  "M1", "M2", "M3", "M4", "M5")

title.lab = c("M1", "M2", "M3", "M4", "M1/M2", "M2/M3", "M3/M4", "M1/M2/M3")
names(title.lab) = c("1000", "0100", "0010", "0001", "1100", "0110", "0011", "1110")

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Top genes for each Metagene by decreased logFC (from Limma)
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
ntop = 10
font.size = 8
top.title.s = 1

top.th0 = mg_tops(
  .obj = ftrs.th0.wrk,
  .meta.study = meta.studies$Th0,
  .u.cl = c("1000", "0100", "0010", "0001"),
  .ftrs.ave = ftrs.ave[, grepl("Th0_", colnames(ftrs.ave))],
  .metafor.res = meta.res.df,
  .ntop = ntop,
  .col = metag.col.7,
  .legend.title = "Mean expression (rank normalized)",
  .font.size = font.size)

top.itreg = mg_tops(
  .obj = ftrs.itreg.wrk,
  .meta.study = meta.studies$iTreg,
  .u.cl = c("1000", "0100", "0010", "0001"),
  .ftrs.ave = ftrs.ave[, grepl("iTreg_", colnames(ftrs.ave))],
  .metafor.res = meta.res.df,
  .ntop = ntop,
  .col = metag.col.7,
  .legend.title = "Mean expression (rank normalized)",
  .font.size = font.size)

top.th17 = mg_tops(
  .obj = ftrs.th17.wrk,
  .meta.study = meta.studies$Th17,
  .u.cl = c("1000", "0100", "0010", "0001"),
  .ftrs.ave = ftrs.ave[, grepl("Th17_", colnames(ftrs.ave))],
  .metafor.res = meta.res.df,
  .ntop = ntop,
  .col = metag.col.7,
  .legend.title = "Mean expression (rank normalized)",
  .font.size = font.size)

top.th2 = mg_tops(
  .obj = ftrs.th2.wrk,
  .meta.study = meta.studies$Th2,
  .u.cl = c("10000", "01000", "00100", "00010", "00001"),
  .ftrs.ave = ftrs.ave[, grepl("Th2_", colnames(ftrs.ave))],
  .metafor.res = meta.res.df,
  .ntop = ntop,
  .col = metag.col.7,
  .legend.title = "Mean expression (rank normalized)",
  .font.size = font.size)

top.th1 = mg_tops(
  .obj = ftrs.th1.wrk,
  .meta.study = meta.studies$Th1,
  .u.cl = c("100", "010", "001"),
  .ftrs.ave = ftrs.ave[, grepl("Th1_", colnames(ftrs.ave))],
  .metafor.res = meta.res.df,
  .ntop = ntop,
  .col = metag.col.7,
  .legend.title = "Mean expression (rank normalized)",
  .font.size = font.size)

x.limits =  max(c(
  max(abs(top.th0$data$logFC)),
  max(abs(top.itreg$data$logFC)),
  max(abs(top.th17$data$logFC)),
  max(abs(top.th2$data$logFC)),
  max(abs(top.th1$data$logFC))
))

s = .5

topX.legend = cowplot::get_legend(top.th0)
top.th0 = top.th0 + xlim(c(-x.limits - s, x.limits + s)) + theme(legend.position = "none")
top.th0 = top.th0 + ggtitle("Th0") + theme(plot.title = element_text(hjust = 0.5, size = rel(top.title.s)))

top.itreg = top.itreg + xlim(c(-x.limits - s, x.limits +  s)) + theme(legend.position = "none")
top.itreg = top.itreg + ggtitle("iTreg") + theme(plot.title = element_text(hjust = 0.5, size = rel(top.title.s)))

top.th17 = top.th17 + xlim(c(-x.limits - s, x.limits + s)) + theme(legend.position = "none")
top.th17 = top.th17 + ggtitle("Th17") + theme(plot.title = element_text(hjust = 0.5, size = rel(top.title.s)))

top.th2 = top.th2 + xlim(c(-x.limits - s, x.limits + s)) + theme(legend.position = "none")
top.th2 = top.th2 + ggtitle("Th2") + theme(plot.title = element_text(hjust = 0.5, size = rel(top.title.s)))

top.th1 = top.th1 + xlim(c(-x.limits - s, x.limits + s)) + theme(legend.position = "none")
top.th1 = top.th1 + ggtitle("Th1") + theme(plot.title = element_text(hjust = 0.5, size = rel(top.title.s)))

# (top.th0 | top.th17 | top.itreg  | top.th2 | top.th1)

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Plot metagene profile (H) over time
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
kinetic_mg_profile = function(obj, rank, metag.col, font.size = 8) {

  nmf.obj = obj$estim.r
  se.obj = obj$se
  h = coef(nmf.obj$fit[[rank]])
  colnames(h) = colData(se.obj)$HOURS
  rownames(h) = paste0("M_", gsub(" ", "", obj$anno$binary))

  for (i in 1:ncol(h)) {
    h[ ,i] = h[ ,i] / sum(h[ ,i])
  }

  h.melt = reshape2::melt(as.matrix(h))
  colnames(h.melt) = c("METAGENE", "HOURS", "WEIGHT")
  # h.melt$HOURS = paste0(h.melt$HOURS, "h")
  h.melt$HOURS = factor(h.melt$HOURS, levels = c("0","0.5","1","2", "4", "6", "8", "12","24","48","72"))
  h.melt$METAGENE = as.character(h.melt$METAGENE)
  h.melt$METAGENE = gsub("M_", "", h.melt$METAGENE)
  h.melt$METAGENE = factor(h.melt$METAGENE, levels = unique(h.melt$METAGENE))

  ggplot(h.melt, aes(x = HOURS, y = WEIGHT, group = METAGENE, colour = METAGENE)) +
    stat_summary(fun = median, fun.min = median, fun.max = median, geom = "line", size = .5) +
    geom_pointrange(mapping = aes(x = HOURS, y = WEIGHT),
                    stat = "summary",
                    fun.min = function(z) {quantile(z,0.25)},
                    fun.max = function(z) {quantile(z,0.75)},
                    fun = median, size = .1) +
    xlab("Hours after activation") +
    ylab("Metagene weight") +
    mytheme(base_size = font.size) +
    scale_y_continuous(limits=c(0, 1), breaks = c(0, 1)) +
    scale_color_manual(values = metag.col, labels = unique(names(mg.lbls)[match(h.melt$METAGENE, mg.lbls)])) +
    theme(
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = rel(1)),
      legend.position = "none",
      aspect.ratio = .618,
      axis.title.x = element_text(margin = margin(3,0,0,0), size = rel(1)),
      plot.margin = margin(0, 0, .5, 0, "cm"),
      # axis.title.x = element_text(size = rel(1)),
      plot.title = element_text(hjust = 0.5),
      axis.ticks = element_blank())
}

exp.pr.th0 = kinetic_mg_profile(obj = nmf.th0, rank = rank.th0, metag.col = mg.col) + ggtitle("Th0") +  xlab("")
exp.pr.th17 = kinetic_mg_profile(obj = nmf.th17, rank = rank.th17, metag.col = mg.col) + ggtitle("Th17")  +  xlab("") + ylab("")
exp.pr.itreg = kinetic_mg_profile(obj = nmf.itreg, rank = rank.itreg, metag.col = mg.col) + ggtitle("iTreg") + ylab("")
exp.pr.th2 = kinetic_mg_profile(obj = nmf.th2, rank = rank.th2, metag.col = mg.col) + ggtitle("Th2")  +  xlab("") + ylab("")
exp.pr.th1 = kinetic_mg_profile(obj = nmf.th1, rank = rank.th1, metag.col = mg.col) + ggtitle("Th1")  +  xlab("") + ylab("")

exp.pr.v1 = kinetic_mg_profile(obj = nmf.v1, rank = rank.v1, metag.col = mg.col) + ggtitle("Memory T-cell Verificaton Set")
exp.pr.v2 = kinetic_mg_profile(obj = nmf.v2, rank = rank.v2, metag.col = mg.col) + ggtitle("Pan T-cell Verificaton Set")  + ylab("")

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Fig Main
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
leg.df = data.frame(
  x = c(1,2,3,4,5),
  label = c(
    "M1: sustained repressed\n0h/0.5h",
    "M2: intermediate\n2h",
    "M3: intermediate\n4h/6h/12h",
    "M4: intermediate\n24h",
    "M5: late\n72h"
  )
)

legend =
  cowplot::get_legend(
  ggplot(leg.df, aes(x, fill = label)) +
    geom_bar() +
    scale_fill_manual(values = c("#000000", "#004488", "#DDAA33", "#117733", "#BB5566")) +
    labs(fill = "Metagene pattern peaks") +
    mytheme(base_size = 8) +
    theme(legend.direction="horizontal", legend.key.height = unit(.2,"cm"), legend.box="horizontal")
)

r1 = .9
r2 = 3.5

gr = plot_grid(
  plot_grid(exp.pr.th0, top.th0, ncol = 1, rel_heights= c(r1, r2), labels = c("A", "B"), label_fontface = "bold", label_size = 10, vjust  = -1),
  plot_grid(exp.pr.th17, top.th17, ncol = 1, rel_heights= c(r1, r2)),
  plot_grid(exp.pr.itreg, top.itreg, ncol = 1, rel_heights= c(r1, r2)),
  plot_grid(exp.pr.th2, top.th2, ncol = 1, rel_heights= c(r1, r2)),
  plot_grid(exp.pr.th1, top.th1, ncol = 1, rel_heights= c(r1, r2)),
  ncol = 5, scale = .99
)
# gr = exp.pr.th0 / top.th0 + plot_layout(heights  = c(r1, r2)) |
#   exp.pr.th17 / top.th17 + plot_layout(heights  = c(r1, r2)) |
#   exp.pr.itreg / top.itreg + plot_layout(heights  = c(r1, r2)) |
#   exp.pr.th2 / top.th2 + plot_layout(heights  = c(r1, r2)) |
#   exp.pr.th1 / top.th1  + plot_layout(heights  = c(r1, r2))

gr =
  ggdraw(add_sub(gr, "Log2 fold change", size = 8, y = 1, vpadding = grid::unit(.25, "lines")))

if (ntop == 10) {
  ggsave2(
    filename="../figures/main/fig_3.pdf",
    plot = ggdraw(legend) / gr + topX.legend + plot_layout(heights  = c(.05, .9, .03)),
    width = 297,
    height = 240,
    dpi = 300,
    units = "mm"
  )
  ggsave2(
    filename="../figures/main/fig_3.png",
    plot = ggdraw(legend) / gr + topX.legend + plot_layout(heights  = c(.05, .9, .03)),
    width = 297,
    height = 240,
    dpi = 300,
    units = "mm"
  )
}

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Supps (V1 and V2)
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
leg.df = data.frame(
  x = c(1,2,3,4,5),
  label = c(
    "M1: sustained repressed\n0h",
    "M2: intermediate\n2h/4h",
    "M3: intermediate\n/6h/8h/12h",
    "M4: intermediate\n24h",
    "M5: late\n72h"
  )
)

meta.prep.all = readRDS(paste0(work, manifest$meta_se_w_izi_arcelus))
exprs.cpm.all = assays(meta.prep.all)$cpm

colnames(exprs.cpm.all) = gsub(".+\\.", "", colnames(exprs.cpm.all))
ftrs.ave.all = grp_ave_sd(.obj = exprs.cpm.all, .ave = "mean")$ave
rownames(ftrs.ave.all) = gc.anno$GENE_SYMBOL_DUPL_MARKED[match(rownames(ftrs.ave.all), gc.anno$ENSEMBL_ID_ABBR)]

legend = cowplot::get_legend(
  ggplot(leg.df, aes(x, fill = label)) +
    geom_bar() +
    scale_fill_manual(values = c("#000000", "#004488", "#DDAA33", "#117733", "#BB5566")) +
    labs(fill = "Metagene pattern peaks") +
    mytheme(base_size = 8) +
    theme(legend.direction="horizontal", legend.key.height = unit(.2,"cm"), legend.box="horizontal")
)

ntop = 25
font.size = 9
top.title.s = 1

top.v1 = mg_tops(
  .obj = ftrs.work.v1,
  .meta.study = meta.studies$Th0Cd4Mem,
  .u.cl = c("10000", "01000", "00100", "00010", "00001"),
  .ftrs.ave = ftrs.ave.all[, grepl("Th0Cd4Mem_", colnames(ftrs.ave.all))],
  .metafor.res = meta.res.df,
  .ntop = 20,
  .col = metag.col.7,
  .legend.title = "Mean expression (rank normalized)",
  .font.size = font.size)

top.v2 = mg_tops(
  .obj = ftrs.work.v2,
  .meta.study = meta.studies$Th0Cd2,
  .u.cl = c("100", "010", "001"),
  .ftrs.ave = ftrs.ave.all[, grepl("Th0Cd2", colnames(ftrs.ave.all))],
  .metafor.res = meta.res.df,
  .ntop = 20,
  .col = metag.col.7,
  .legend.title = "Mean expression (rank normalized)",
  .font.size = font.size)

x.limits =  max(c(
  max(abs(top.v1$data$logFC)),
  max(abs(top.v2$data$logFC))
))

s = .5

topX.legend = cowplot::get_legend(top.v1)
top.v1 = top.v1 + xlim(c(-x.limits - s, x.limits + s)) + theme(legend.position = "none")
top.v1 = top.v1 + ggtitle("Memory T-cell Verificaton Set") + theme(plot.title = element_text(hjust = 0.5, size = rel(top.title.s)))
top.v2 = top.v2 + xlim(c(-x.limits - s, x.limits +  s)) + theme(legend.position = "none")
top.v2 = top.v2 + ggtitle("Pan T-cell Verificaton Set") + theme(plot.title = element_text(hjust = 0.5, size = rel(top.title.s)))

r1 = .7
r2 = 4
gr = plot_grid(
  plot_grid(exp.pr.v1, top.v1, ncol = 1, rel_heights= c(r1, r2), labels = c("A", "B"), label_fontface = "bold", label_size = 10, vjust  = -1),
  plot_grid(exp.pr.v2, top.v2, ncol = 1, rel_heights= c(r1, r2)),
  ncol = 2, scale = .95
)

ggsave2(filename="../figures/additional/fig_3_top20_v1v2.pdf",
        # plot = ggdraw(legend) / plot_grid(exp.pr.v1, NULL, exp.pr.v2, rel_widths = c(1, .1, 1), nrow = 1) + plot_layout(heights  = c(.15, .85)),
        plot = ggdraw(legend) / gr + topX.legend + plot_layout(heights  = c(.02, .9, .02)),
        width = 160,
        height =400,
        dpi = 100,
        units = "mm")

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
if (ntop != 10) {
  r1 = .5
  r2 = 4

  gr = plot_grid(
    plot_grid(exp.pr.th0, top.th0, ncol = 1, rel_heights= c(r1, r2), labels = c("A", "B"), label_fontface = "bold", label_size = 10, vjust  = -1),
    plot_grid(exp.pr.th17, top.th17, ncol = 1, rel_heights= c(r1, r2)),
    plot_grid(exp.pr.itreg, top.itreg, ncol = 1, rel_heights= c(r1, r2)),
    plot_grid(exp.pr.th2, top.th2, ncol = 1, rel_heights= c(r1, r2)),
    plot_grid(exp.pr.th1, top.th1, ncol = 1, rel_heights= c(r1, r2)),
    ncol = 5, scale = .99
  )

  ggsave2(filename="../figures/additional/fig_3_top20.pdf",
          plot = ggdraw(legend) / gr + topX.legend + plot_layout(heights  = c(.02, .9, .02)),
          width = 320,
          height = 450,
          dpi = 100,
          units = "mm")

}
