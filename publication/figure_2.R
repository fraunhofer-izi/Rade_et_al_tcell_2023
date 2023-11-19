# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# LIBRARIES
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
.cran_packages = c("ggplot2","devtools", "yaml", "naturalsort", "dplyr", "data.table",
                   "metafor", "reshape2", "tidyr", "ggrepel", "cowplot", "broom", "ggrastr")
.bioc_packages = c("SummarizedExperiment", "limma", "ComplexHeatmap", "ReactomePA")

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

library("tidymeta")
library("mbmisc")

library(clusterProfiler)
library(ReactomePA)

source("paper-helper.R")
source("code/eda/meta-effect//helper.R")

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Read data
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
manifest = yaml.load_file("assets/manifest.yaml")
work = manifest$workdata
obj = readRDS(paste0(work, manifest$meta_all))
meta.in = obj$meta.in
meta.res = obj$meta.res
meta.res.obj = readRDS(paste0(work, manifest$meta_objects))
se = readRDS(paste0(work, manifest$meta_se))

m.res = lapply(meta.res, function(ctrst) {
  ctrst %>%
    dplyr::select(GENE_SYMBOL, est_pvalue, est_pvalue_adj, est_ci_l, est_ci_r,
                  est_effect, est_se, het_I2, rank, confect, confect_fdr_zero,
                  GROUP, signcon, ntimes)
})
m.res.df = do.call("rbind", m.res)
m.res.df$ID = paste0(m.res.df$GENE_SYMBOL, "_", m.res.df$GROUP)

lfc = 0
font.size = 8

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Top genes from meta-analysis (based on Topconfects)
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
tops = subset(m.res.df, abs(confect) > lfc)
tops$ID = paste0(tops$GENE_SYMBOL, "_", tops$GROUP)
tops.1 = tops[tops$GROUP == "05h" |tops$GROUP == "1h" |tops$GROUP == "2h" | tops$GROUP == "4h" |tops$GROUP == "6h", ]
tops.1 = subset(tops.1, ntimes >= 4)
tops.2 = tops[tops$GROUP == "12h" |tops$GROUP == "24h" |tops$GROUP == "48h" |tops$GROUP == "72h", ]
tops.2 = subset(tops.2, ntimes >= 5)
tops = rbind(tops.1, tops.2)

tops$GROUP = gsub("05h", "0.5h", tops$GROUP)

tops = tops %>% group_by(GROUP) %>% arrange(rank) %>% slice_head(n=10)  %>% data.frame()
tops = tops[naturalorder(tops$GROUP), ]
# tops = tops[seq(dim(tops)[1],1),] # Reverse

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Rev. cumulative distrobution
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# ctrst = "12h"
rev.cum.sum = lapply(names(meta.in), function(ctrst){
  if (ctrst %in% c("12h", "24h", "48h", "72h")) {
    nbr = 5
  } else {
    nbr = 4
  }

  df = votecount_mv(meta.in[[ctrst]])
  df = cum_freq_data(df, nbr)
  df$TIMEPOINT = ctrst
  df
})

rev.cum.sum.m = do.call("rbind", rev.cum.sum)
rev.cum.sum.m$TIMEPOINT = gsub("05h", "0.5h", rev.cum.sum.m$TIMEPOINT)
rev.cum.sum.m$TIMEPOINT = gsub("h", "", rev.cum.sum.m$TIMEPOINT)
lvls = unique(naturalsort(rev.cum.sum.m$TIMEPOINT))
rev.cum.sum.m$TIMEPOINT = factor(rev.cum.sum.m$TIMEPOINT, levels = lvls)
rev.cum.sum.m = rev.cum.sum.m[rev.cum.sum.m$ndatasets != 0, ]
rev.cum.sum.m$DEGs = rev.cum.sum.m$DEGs / 1000

cust.leg.labels = c(
  expression("n">="1"),
  expression("n">="2"),
  expression("n">="3"),
  expression("n">="4"),
  expression("n=5")
)
names(cust.leg.labels) = c(1, 2, 3, 4, 5)

rev.c.pl =
  ggplot(rev.cum.sum.m, aes(x = TIMEPOINT, y = DEGs, colour = as.factor(ndatasets), group = as.factor(ndatasets), label = DEGs)) +
  geom_line(size = .5) +
  geom_point(size = .7) +
  mytheme_grid(base_size = font.size) +
  theme(
    aspect.ratio = 1,
    axis.ticks = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "bottom",
    legend.key.size=unit(.1,"cm")
  ) +
  ggtitle("\n") +
  scale_color_manual(
    values = c("#555555", "#6699CC", "#997700", "#EE99AA", "#BB5566"),
    labels = cust.leg.labels
  ) +
  guides(color = guide_legend(
    title.position = "top",
    title = expression(paste("No. of ", CD4^{"+"},' populations')),
    nrow = 1,
    override.aes = aes(label = ""))
  ) +
  labs(
    x = "Hours after activation",
    y = "No. of DE genes (thousand)"
  ) +
  ylim(0, 18.1) +
  geom_text_repel(
    data = (subset(rev.cum.sum.m, ndatasets != 5 & TIMEPOINT != .5)),
    size = 2.5,
    nudge_y = 1,
    segment.size = .2,
    box.padding = 0.15,
    min.segment.length = 0
  ) +
    geom_text_repel(
      data = (subset(rev.cum.sum.m, ndatasets == 5 & TIMEPOINT != .5)),
      size = 2.5,
      nudge_y = -1.5,
      # nudge_x = -.2,
      segment.size = .2,
      box.padding = 0.15,
      min.segment.length = 0
    ) +
    geom_text_repel(
      data = (subset(rev.cum.sum.m, TIMEPOINT == .5)),
      size = 2.5,
      nudge_y = 6,
      # nudge_x = -.2,
      segment.size = .2,
      box.padding = 0.15,
      min.segment.length = 0
    )

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Vote-counting
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
vote.df = lapply(meta.in, function(x){
  df = votecount_mv(x)
  df %>% dplyr::select(ddeg, ndeg)
})
vote.df = do.call("rbind", vote.df)

label.df = vote.df %>%
  dplyr::count(ndeg, ddeg) %>%
  data.frame()

vote.pl =
  ggplot(vote.df, aes(x = ddeg, col = as.factor(ndeg), y = as.factor(ndeg))) +
  geom_jitter(shape = ".", width = 0.45, height = 0.45, alpha = .3) +
  scale_color_manual(values = c("black", "#959595", "#004488", "#117733", "#BB5566")) +
  scale_x_continuous(
    limits = c(-5.5, 5.5),
    breaks = c(-5, -4, -3, -2, -1, 0, 1, 2, 3, 4, 5),
    labels = c(-5, -4, -3, -2, -1, 0, 1, 2, 3, 4, 5)
  ) +
  labs(
    x = "Fold change sign consistency",
    y = "No. of times as DE"
  ) +
  mytheme(base_size = font.size) +
  theme(
    aspect.ratio = 1,
    axis.ticks = element_blank(),
    legend.position = "none"
  ) +
  geom_label(data = label.df, aes(label = n), size = 2.3) +
  ggtitle("\n")

vote.pl = rasterize(vote.pl, layers='Point', dpi=200)

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# logFC vs confect;  I^2 and results
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
y.axis.max = max(max(m.res.df$confect, na.rm = T), abs(min(m.res.df$confect, na.rm = T)))

tmp = subset(m.res.df, abs(confect) > lfc)
tmp$LFC_DIRECTION = tmp$confect > 0
tmp$ntimes.pl = paste0("k=", tmp$ntimes)
tmp = as.data.frame.matrix(table(tmp$ntimes.pl, tmp$LFC_DIRECTION))
colnames(tmp) = c("down", "up")

pl.table = ggtexttable(tmp, theme = ttheme("blank", base_size = 6.5))

tops.pl = tops[order(abs(tops$confect), decreasing = T), ]
# tops.pl.1 = subset(tops.pl, GROUP == "6h")[1:3, ]
# tops.pl.2 = subset(tops.pl, GROUP == "72h")[1:3, ]
tops.pl.1 = subset(tops.pl, ntimes == 4)[1:5, ]
tops.pl.2 = subset(tops.pl, ntimes == 5)[1:5, ]
tops.pl = rbind(tops.pl.1, tops.pl.2)

set.seed(42)
confect.logFC =
  m.res.df %>% ggplot(aes(x = est_effect, y = ifelse(is.na(confect), 0, confect), col = as.factor(ntimes))) +
  geom_point(shape = ".", alpha = .8) +
  labs(x = "Combined effect size", y = "Confect") +
  mytheme_grid(base_size = font.size) +
  theme(
    aspect.ratio = 1,
    axis.ticks = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "none"
  ) +
  scale_color_manual(values = c("#959595", "#004488", "#117733", "#BB5566")) +
  scale_y_continuous(limits=c(min(m.res.df$confect, na.rm = T), abs(y.axis.max)), breaks = c(-4, -2, 0, 2, 4), labels = c(-4, -2, 0, 2, 4)) +
  geom_label_repel(
    data = m.res.df[m.res.df$ID %in% tops.pl.2$ID, ], label = paste0(tops.pl.2$GENE_SYMBOL, " (", tops.pl.2$GROUP, ")"),
    size = 1.9,
    direction = "y",
    hjust = .5,
    nudge_x = -30,
    segment.size = .2,
    box.padding = 0.15,
    label.padding = .12,
    fill = NA,
    min.segment.length = 0
  ) +
  geom_label_repel(
      data = m.res.df[m.res.df$ID %in% tops.pl.1$ID, ], label = paste0(tops.pl.1$GENE_SYMBOL, " (", tops.pl.1$GROUP, ")"),
      size = 1.9,
      direction = "y",
      hjust = 0.5,
      nudge_x = 10,
      nudge_y = -2,
      segment.size = .2,
      box.padding = 0.15,
      label.padding = .12,
      fill = NA,
      min.segment.length = 0
    ) +
  xlim(c(min(m.res.df$est_effect) - 0, max(m.res.df$est_effect) + 0))

confect.logFC = rasterize(confect.logFC, layers='Point', dpi=200)

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Plot results
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Legend
leg.df = vote.df
leg.df$ndeg = paste0("n=", leg.df$ndeg)
leg.df$ddeg = paste0("d=", leg.df$ddeg)

legend = cowplot::get_legend(
  ggplot(leg.df, aes(ddeg, ddeg, col = ndeg)) +
  geom_point() +
  scale_colour_manual(values = c("black", "#959595", "#004488", "#117733", "#BB5566")) +
  labs(col = expression(paste("No. of ", CD4^{"+"},' populations'))) +
  mytheme(base_size = font.size) +
    theme(legend.direction="horizontal", legend.key.height = unit(.2,"cm"), legend.box="horizontal")+
  guides(colour = guide_legend(nrow = 1, title.position = "top", override.aes = list(size=3)))
)

a.b.c = plot_grid(
  rev.c.pl + theme(legend.position = "none"),
  NULL,
  vote.pl,
  NULL,
  confect.logFC,
  nrow = 1,
  rel_widths = c(1, 0.1, 1, 0.1, 1),
  align = "vh",
  labels = c("A", "", "B", "", "C"),
  label_fontface = "bold",
  label_size = 10,
  vjust  = 2
)

a.b.c = plot_grid(
  a.b.c,
  NULL,
  plot_grid(NULL, get_legend(rev.c.pl), NULL, legend, NULL, nrow = 1, rel_widths = c(.035, 1, .55, 1, .5)),
  nrow = 3, rel_heights = c(.92, .02, .12)
)

# a.b.c = plot_grid(a.b.c, NULL, legend, nrow = 3, rel_heights = c(.92, .02, .12))

set.seed(1234)
ggsave(
  filename="../figures/main/fig_2.png",
  a.b.c,
  width = 210,
  height = 82,
  dpi = 300,
  bg = "white",
  units = "mm"
)

set.seed(12346)
ggsave(
  filename="../figures/main/fig_2.pdf",
  a.b.c,
  width = 210,
  height = 82,
  dpi = 300,
  bg = "white",
  units = "mm"
)
