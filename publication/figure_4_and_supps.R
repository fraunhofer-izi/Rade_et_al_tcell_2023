.cran_packages = c("miscTools", "yaml", "NMF", "ggrepel", "ggplot2", "dplyr", "data.table",
                   "viridis", "cowplot", "naturalsort", "scico", "scales", "inlmisc")
.bioc_packages = c("SummarizedExperiment", "pheatmap", "ComplexHeatmap")

# backports, vctrs

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
library(ReactomePA)

# source("code/R/eda-plots.R")
source("paper-helper.R")
theme_set(mytheme())

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
manifest = yaml.load_file("assets/manifest.yaml")
work = manifest$workdata

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# NMF results
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# After assignment of genes to Metagenes
ftrs.work = readRDS(paste0(work, manifest$pipeline_discovery))
ftrs.th0.wrk = ftrs.work$ftrs.th0.wrk
ftrs.th17.wrk = ftrs.work$ftrs.th17.wrk
ftrs.itreg.wrk = ftrs.work$ftrs.itreg.wrk
ftrs.th2.wrk = ftrs.work$ftrs.th2.wrk
ftrs.th1.wrk = ftrs.work$ftrs.th1.wrk

nmf.th0 = ftrs.work$nmf.th0
nmf.th17 = ftrs.work$nmf.th17
nmf.itreg = ftrs.work$nmf.itreg
nmf.th2 = ftrs.work$nmf.th2
nmf.th1 = ftrs.work$nmf.th1
rank.th0 = ftrs.work$rank.th0
rank.th17 = ftrs.work$rank.th17
rank.itreg = ftrs.work$rank.itreg
rank.th2 = ftrs.work$rank.th2
rank.th1 = ftrs.work$rank.th1

ftrs.work.v = readRDS(paste0(work, manifest$pipeline_verification_v1))
ftrs.v1.wrk = ftrs.work.v$ftrs.v1.wrk
nmf.v1 = ftrs.work.v$nmf.v1
rank.v1 = ftrs.work.v$rank.v1

# Results form metafor
obj.meta = readRDS(paste0(work, manifest$meta_all))
# meta.res = obj.meta$meta.res

gc.anno = readRDS(paste0(work, manifest$gencode_v29_features))

signtrs = readRDS(paste0(work, manifest$signatures_all))

run = "d"

if (run == "d") {
  q.signtrs = subset(signtrs,  Housekeeping == TRUE)
  table(q.signtrs$CLUSTER)
  dim(q.signtrs)
  meta.prep = readRDS(paste0(work, manifest$meta_se))
  rownames(meta.prep) = rowData(meta.prep)$GENE_SYMBOL_DUPL_MARKED

  ora.obj.go = readRDS("enrichment_analysis/figure_4_ora_go_d.rds")
  ora.obj.reactome = readRDS("enrichment_analysis/figure_4_ora_reactome_d.rds")
}
if (run == "v1") {
  q.signtrs = subset(signtrs, v1_pass == TRUE & Housekeeping == TRUE)
  meta.prep = readRDS(paste0(work, manifest$meta_se_w_izi_arcelus))
  rownames(meta.prep) = rowData(meta.prep)$GENE_SYMBOL_DUPL_MARKED
  meta.prep = meta.prep[, !grepl("Th0Cd2" , colnames(meta.prep))]

  ora.obj.go = readRDS("enrichment_analysis/figure_4_ora_go_v1.rds")
  ora.obj.reactome = readRDS("enrichment_analysis/figure_4_ora_reactome_v1.rds")

  supp.fig = "fig_4_v1"
}
if (run == "v2") {
  q.signtrs = subset(signtrs, v1_pass == TRUE & v2_pass == TRUE & v2_ctr_pass == TRUE & Housekeeping == TRUE)
  q.signtrs.v1 = subset(signtrs, v1_pass == TRUE & Housekeeping == TRUE)
  meta.prep = readRDS(paste0(work, manifest$meta_se_w_izi_arcelus))
  rownames(meta.prep) = rowData(meta.prep)$GENE_SYMBOL_DUPL_MARKED

  ora.obj.go = readRDS("enrichment_analysis/figure_4_ora_go_v1v2.rds")
  ora.obj.reactome = readRDS("enrichment_analysis/figure_4_ora_reactome_v1v2.rds")

  supp.fig = "fig_4_v1v2"
}

q.signtrs$LABEL = mg.labels[match(q.signtrs$CLUSTER, names(mg.labels))]

rank_vector <- function(x) {
  r = x
  r[r!=0] <- rank(r[r!=0])
  # r = dense_rank(x) - 1
  r / max(r)
}
meta.exprs =  apply(2^assays(meta.prep)$cpm, 2, rank_vector)

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Coef & Basis
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
fit.th0 = nmf.th0$estim.r$fit[[rank.th0]]
h.th0 = coef(fit.th0);
# h.th0 = insertRow(h.th0, 4, NA)

fit.th17 = nmf.th17$estim.r$fit[[rank.th17]]
h.th17 = coef(fit.th17);
# h.th17 = insertRow(h.th17, 4, NA)

fit.itreg = nmf.itreg$estim.r$fit[[rank.itreg]]
h.itreg = coef(fit.itreg)
# h.itreg = insertRow(h.itreg, 4, NA)

fit.th2 = nmf.th2$estim.r$fit[[rank.th2]]
h.th2 = coef(fit.th2)
h.th2 = h.th2[c(1,2,3,5), ]

fit.th1 = nmf.th1$estim.r$fit[[rank.th1]]
h.th1 = coef(fit.th1)
h.th1 = insertRow(h.th1, 2, NA);
# h.th1 = insertRow(h.th1, 4, NA)

fit.v1 = nmf.v1$estim.r$fit[[rank.v1]]
h.v1 = coef(fit.v1)
# h.v1 = insertRow(h.v1, 4, NA)

dim(h.th0);dim(h.th17);dim(h.itreg);dim(h.th2);dim(h.th1);dim(h.v1)
h.merge = cbind(h.th0, h.th17, h.itreg, h.th2, h.th1)
# rownames(h.merge) = paste0("M", seq(1, nrow(h.merge)))
rownames(h.merge) = c("M1", "M2", "M3", "M5")
###

b.th0 = basis(fit.th0); # rownames(b.th0) = paste0("Th0_", rownames(b.th0))
b.th17 = basis(fit.th17); # rownames(b.th17) = paste0("Th17_", rownames(b.th17))
b.itreg = basis(fit.itreg); # rownames(b.itreg) = paste0("iTreg_", rownames(b.itreg))
b.th1 = basis(fit.th1); b.th1 = insertCol(b.th1, 2, NA); # rownames(b.th1) = paste0("Th1_", rownames(b.th1))
b.th2 = basis(fit.th2)[, c(1,2,3,5)]; # rownames(b.th2) = paste0("Th2_", rownames(b.th2))

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Metagene landscape
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# alpha.exp: Increase this > 1.0 to move the cells closer to the factors. Values > 2 start to distort the data.
swne.embedding = embedding_pattern(H = h.merge, alpha.exp = 1, proj.method = "sammon", dist.use = "euclidean", rescale.coords = "scale")
swne.embedding = ftr.embedding(swne.embedding, b.th0, rownames(b.th0), scale.cols = T, overwrite = T, alpha.exp = 1)
swne.embedding = ftr.embedding(swne.embedding, b.th17, rownames(b.th17), scale.cols = T, overwrite = F, alpha.exp = 1)
swne.embedding = ftr.embedding(swne.embedding, b.itreg, rownames(b.itreg), scale.cols = T, overwrite = F, alpha.exp = 1)
swne.embedding = ftr.embedding(swne.embedding, b.th1, rownames(b.th1), scale.cols = T, overwrite = F, alpha.exp = 1)
swne.embedding = ftr.embedding(swne.embedding, b.th2, rownames(b.th2), scale.cols = T, overwrite = F, alpha.exp = 1)

sample.coords = swne.embedding$sample.coords
sample.coords = data.frame(merge(sample.coords, colData(meta.prep), by = "row.names"))

H.coords = swne.embedding$H.coords
F.coords = swne.embedding$feature.coords
# F.coords = F.coords[F.coords$name %in% q.signtrs$GENE_SYMBOL, ]

F.coords.ave = F.coords %>%
  dplyr::group_by(name) %>%
  dplyr::summarise(across(everything(), list(mean)))
colnames(F.coords.ave) = c("name", "x", "y")

samp.coords.size = 1.5
samp.coords.sh.size = 2
h.coords.size = 6
h.coords.sh.size = 6.5
emb.pl.size = 8

sample.coords$HOURS = gsub("0h", "0", sample.coords$HOURS)
sample.coords$HOURS = gsub("h", "", sample.coords$HOURS)
sample.coords$HOURS = gsub("05", "0.5", sample.coords$HOURS)
lvls = naturalsort(unique(sample.coords$HOURS))
sample.coords$HOURS = factor(sample.coords$HOURS, levels = lvls)

sample.emb.pl =
  ggplot() +
  stat_density_2d(data = F.coords, aes(x, y, fill = ..density..), alpha = .2, geom = "raster", contour = F, show.legend = F) +
  scale_fill_scico(palette = "roma", direction = -1, begin = 0, end = .8) +
  geom_point(data = sample.coords, aes(x, y, colour = HOURS), alpha = 1, size = samp.coords.size) +
  geom_point(data = sample.coords, aes(x, y), shape = 1, size = samp.coords.sh.size, colour = "#555555", stroke = .1, alpha = .7) +
  geom_point(data = H.coords, aes(x, y), size = h.coords.size, color = c("#000000", "#004488", "#DDAA33", "#BB5566")) + # "#BB5566"
  geom_point(data = H.coords, aes(x, y), shape = 1, size = h.coords.sh.size, colour = "#555555", stroke = .1) +
  geom_text(data = H.coords, aes(x, y, label=name), size = 3, colour = "white" ) +
  geom_point(data = subset(F.coords.ave, name == "CD4" | name == "IL2" | name == "IL2RA" | name == "NFKBID"), aes(x, y), alpha = 1, size = 1) +
  ggrepel::geom_label_repel(data = subset(F.coords.ave, name == "CD4" | name == "IL2" | name == "IL2RA" | name == "NFKBID"), mapping = aes(x, y, label = name), fill = "white", colour = "#555555", size = 3, box.padding = .5, label.padding = 0.15, min.segment.length = unit(0, 'lines'), nudge_y = .05) +
  # ggrepel::geom_text_repel(data = H.coords, mapping = aes(x, y, label = name), colour = "#555555", size = 2.5, box.padding = 0.5) +
  mytheme(base_size = emb.pl.size) +
  theme(
    legend.position = "right",
    plot.title = element_text(hjust = 0.5, face = "plain", colour = "black", size = rel(1)),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.ticks = element_blank(),
    axis.line = element_blank(),
    axis.text = element_blank()
  )+
  scale_color_manual(values = c("#BBBBBB", tol9qualitative)) +
  guides(colour = guide_legend(ncol = 1, byrow=TRUE,  override.aes = list(size=3))) +
  labs(col = "Hours after activation") +
  ggtitle(NULL) +
  scale_x_reverse()

ftrs_exprs_emb = function(.sample_coords = NULL, .meta_exprs = NULL, .signtrs = NULL, .ftr = NULL, .col.pal = "romaO") {

  ftrs.exprs = .meta_exprs[.ftr, ]
  pl.title = subset(.signtrs, GENE_SYMBOL == .ftr)
  pl.title = paste0("\n\n", pl.title$GENE_SYMBOL, " | ", pl.title$LABEL)
  .sample_coords$FTR_EXPRS = ftrs.exprs[match(.sample_coords$Row.names,  names(ftrs.exprs))]

  ggplot() +
    geom_point(data = .sample_coords, aes(x, y, colour = FTR_EXPRS), alpha = 1, size = .8) +
    geom_point(data = H.coords, aes(x, y), size = h.coords.size - 2, color = c("#000000", "#004488", "#DDAA33", "#BB5566")) +
    geom_point(data = H.coords, aes(x, y), shape = 1, size = h.coords.sh.size - 2, colour = "#555555", stroke = .1) +
    geom_text(data = H.coords, aes(x, y, label=name), size = 2, colour = "white" ) +
    # ggrepel::geom_text_repel(data = H.coords, mapping = aes(x, y, label = name), colour = "#555555", size = 2.5, box.padding = 0.3) +
    mytheme(base_size = 8) +
    theme(
      legend.position = "right",
      legend.title = element_text(size = rel(.9)),
      plot.title = element_text(hjust = 0.5, face = "plain", colour = "black", size = rel(1)),
      axis.title.x = element_blank(), axis.title.y = element_blank(),
      axis.ticks = element_blank(), axis.line = element_blank(),
      axis.text = element_blank())+
    guides(color = guide_colorbar(title.position = 'top', title = "Gene expression\n(rank normalized)", title.hjust = 0, barwidth = unit(.3, 'lines'), barheight = unit(5, 'lines'))) +
    scale_colour_scico(palette = .col.pal, direction = -1, begin = 0.1, end = 1, breaks= pretty_breaks()) +
    # scale_colour_gradientn(colours = inlmisc::GetColors(256, bias = 1.5)) +
    ggtitle(pl.title)+
    scale_x_reverse()
}

exprs.1 = ftrs_exprs_emb(.sample_coords = sample.coords, .meta_exprs = meta.exprs, .signtrs = q.signtrs, .ftr = "IL2")
exprs.2 = ftrs_exprs_emb(.sample_coords = sample.coords, .meta_exprs = meta.exprs, .signtrs = q.signtrs, .ftr = "NFKBID")
exprs.3 = ftrs_exprs_emb(.sample_coords = sample.coords, .meta_exprs = meta.exprs, .signtrs = q.signtrs, .ftr = "CD4")
exprs.4 = ftrs_exprs_emb(.sample_coords = sample.coords, .meta_exprs = meta.exprs, .signtrs = q.signtrs, .ftr = "IL2RA")

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Kinetic boxplots
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
ftrs = q.signtrs

exprs.cpm = assays(meta.prep)$fsqn

colnames(exprs.cpm) = gsub(".+\\.", "", colnames(exprs.cpm))
exprs.med = aveByGrp(exprs.cpm, .ave = "median")

obj = exprs_med_melt(.exprs.mat = exprs.med, .ftrs = ftrs, .exclude = "Th1_") # also standardize
df = obj$df

if (run == "d") {
  factor.col = c("#88CCEE", "#44AA99","#117733", "#999933", "#DDCC77")
  names(factor.col) = c("Th0", "iTreg", "Th17", "Th2", "Th1")
  legend.title = expression(paste(CD4^{"+"},' population'))
}
if (run == "v1") {
  df[df$CELLTYPE == "Th0Cd4Mem", ]$CELLTYPE = "Th0 (Memory)"
  factor.col = c("#88CCEE", "#4477AA", "#44AA99","#117733", "#999933", "#DDCC77")
  names(factor.col) = c("Th0", "Th0 (Memory)", "iTreg", "Th17", "Th2", "Th1")
  legend.title = expression(paste(CD4^{"+"},' population'))
}
if (run == "v2") {
  df[df$CELLTYPE == "Th0Cd4Mem", ]$CELLTYPE = "Th0 (Memory)"
  df[df$CELLTYPE == "Th0Cd2", ]$CELLTYPE = "Th0 (Pan)"
  factor.col = c("#88CCEE", "#4477AA", "#332288", "#44AA99","#117733", "#999933", "#DDCC77")
  names(factor.col) = c("Th0", "Th0 (Memory)", "Th0 (Pan)", "iTreg", "Th17", "Th2", "Th1")
  legend.title = expression(paste(CD4^{"+"},"/Pan", ' population'))
}

df$CELLTYPE = factor(df$CELLTYPE , levels = names(factor.col))

df$HOURS = gsub("h", "", df$HOURS)
df$HOURS = gsub("05", "0.5", df$HOURS)
lvls = naturalsort(unique(df$HOURS))
df$HOURS = factor(df$HOURS, levels = lvls)


plot.order = names(mg.labels[names(mg.labels) %in% unique(ftrs$CLUSTER)])
title.lab = mg.labels[names(mg.labels) %in% unique(ftrs$CLUSTER)]

ftrs.nbr = table(df[!duplicated(df$GENE_SYMBOL), ]$CLUSTER)
title.lab = title.lab[names(title.lab) %in% names(ftrs.nbr)]
ftrs.nbr = ftrs.nbr[names(title.lab)]
title.lab[names(ftrs.nbr)] = paste0(title.lab, " (", ftrs.nbr, " genes)" )

kinetic_boxplot = function(obj, title = NULL) {
  ggplot(obj, aes(HOURS, EXPRS, group = interaction(CELLTYPE, HOURS))) +
    geom_boxplot(
      aes(color = CELLTYPE, fill = CELLTYPE),
      position = position_dodge(preserve = "single", width=.86),
      outlier.colour = NULL,
      outlier.size = .01,
      lwd = 0.1
     ) +
    stat_summary(
      geom = "crossbar",
      width = 0.6,
      lwd = .1,
      fatten=0,
      color="black",
      position = position_dodge(preserve = "single", width=.86),
      fun.data = function(x){c(y=median(x), ymin=median(x), ymax=median(x))}
  ) +
  guides(colour = "none") +
  scale_fill_manual(values = factor.col) +
  scale_colour_manual(values = factor.col) +
  mytheme(base_size = 8) +
  theme(
    axis.text.x = element_text(angle = 90, vjust = .5, hjust = 1, size = rel(1)),
    plot.title = element_text(hjust = 0.5, face = "bold", colour = "black", size = rel(1.0)),
    axis.ticks.x = element_blank(),
    panel.spacing = unit(1, "lines")
  ) +
  xlab("Hours after activation") + ylab("Standardized expression") + labs(fill = legend.title) +
  ggtitle(title) +
  facet_wrap( ~ CLUSTER, nrow = 2, scales = "free", labeller = labeller(CLUSTER = title.lab))
}

df.identity = df[df$CLUSTER %in% names(title.lab)[1:4], ]
df.identity$CLUSTER = factor(df.identity$CLUSTER, levels = plot.order)
kinetic.pl.identity = kinetic_boxplot(obj = df.identity, title = "Identity Metagenes")

df.identity$COL = mg.col[match(df.identity$CLUSTER, names(mg.col))]
kinetic.pl.identity = kinetic.pl.identity +
  geom_hline(aes(yintercept = Inf), color = factor(df.identity$COL), size=2)

df.shared = df[!df$CLUSTER %in% names(title.lab)[1:4], ]
df.shared$CLUSTER = factor(df.shared$CLUSTER, levels = plot.order)
kinetic.pl.shared = kinetic_boxplot(df.shared, "Shared Metagenes")

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
ora_plot = function(ora.obj, .go = TRUE) {
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
  cp.res = cp.res[cp.res$Count >= 4, ]
  cp.res = cp.res[cp.res$p.adjust < 0.05, ]
  cp.res$geneID = NULL

  cp.res.l = split(cp.res, cp.res$Cluster)

  tops.l = lapply(cp.res.l, function(x){
    x = x[order(x$richFactor, decreasing = T), ]
    # x = x[order(x$p.adjust, decreasing = F), ]
    x = head(x, 10)
    x
  })
  tops.l = base::Filter(function(x) dim(x)[1] > 0, tops.l)

  tops = as.data.frame(data.table::rbindlist(tops.l))
  res.df = subset(cp.res, ID %in% tops$ID)

  res.df$LABEL = mg.labels[match(res.df$Cluster, names(mg.labels))]
  lvls = unname(mg.labels[mg.labels %in% res.df$LABEL])
  res.df$LABEL = factor(res.df$LABEL, levels = lvls)
  l = list()
  for (i in lvls) { l[[i]] = subset(res.df, LABEL == i) }
  res.df = do.call("rbind", l)
  res.df$Description = factor(res.df$Description, levels = rev(unique(res.df$Description)))


  if (.go == TRUE) {
    ggplot(res.df, aes(x = LABEL, y = Description)) +
      geom_point(aes(fill = p.adjust, size = richFactor), pch=21, stroke=.3) +
      scale_size(range = c(2, 5)) +
      mytheme_grid(base_size = 8) +
      theme(
        legend.position="right",
        plot.margin = margin(10, 15, 4, 4, "pt"),
        panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank(),
        axis.ticks = element_blank()
      ) +
      scale_fill_scico(palette = 'bilbao', begin = 0.1, end = .8, direction = -1) +
      # scale_colour_gradient(high = "#332288", low = "#DDCC77", name = "Adjusted\np-value") +
      guides(fill = guide_colorbar(title.position = 'top', title = "Adjusted\np-value", title.vjust = 1, barwidth = unit(.4, 'lines'), barheight = unit(5, 'lines'))) +
      guides(size = guide_legend(ncol = 1, override.aes =, title.hjust = .5, title = "Rich factor", order = 1)) +
      xlab("Metagene") +
      ylab(NULL)
  } else {
    res.df$TOP  = reactomePathways$TOP_2[match(res.df$ID, reactomePathways$ID)]
    top.level = unique(res.df$TOP)
    top.level.col = c(tol12qualitative, "grey", "555555", "black")[1:length(top.level)]
    names(top.level.col) = top.level

    ggplot(res.df, aes(x = LABEL, y = Description, size = richFactor, colour = TOP)) +
      geom_point() +
      scale_size(range = c(2, 5)) +
      mytheme_grid(base_size = 8) +
      theme(
        legend.position="right",
        plot.margin = margin(10, 15, 4, 4, "pt"),
        panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank(),
        axis.ticks = element_blank()
      ) +
      scale_colour_manual(values = top.level.col) +
      guides(colour = guide_legend(
        ncol = 1, override.aes = list(size=2), title.hjust = 0, title =  "Top level of the\nReactome pathway hierarchy"
      )) +
      guides(size = guide_legend(ncol = 1, override.aes =, title.hjust = .5, title = "Rich factor", order = 1)) +
      xlab("Metagene") +
      ylab(NULL)
  }

}
# ora_plot(ora.obj.go)

ora = ora_plot(ora.obj.go)
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Top genes (ranked by confect)
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
ftrs = q.signtrs

# abbr = "M1"
# .cl = "1000"
# .ntop = 10
consensus_tops = function(.cl, .ntop = t, .sort = "confect") {

  abbr = unname(.cl)
  .cl = names(.cl)
  print(.cl)

  st.tops = subset(ftrs, CLUSTER == .cl)
  st.tops = st.tops %>% dplyr::arrange(desc(abs(confect)), .by_group = F)

  df = head(st.tops, n = .ntop)
  df$COL = mg.col[match(df$CLUSTER, names(mg.col))]
  df$ID = paste0(df$GENE_SYMBOL, "_", df$max_time_point)

  df = df %>% mutate(x_max_x = case_when(confect > 0 ~ confect, TRUE ~ -Inf))
  df = df %>% mutate(x_min_x = case_when(confect < 0 ~ confect, TRUE ~ Inf))

  df$GENE_SYMBOL = gsub("\\_dupl", "", df$GENE_SYMBOL)
  df$GENE_LABEL = paste0(df$GENE_SYMBOL, " | ", df$max_time_point)
  df$CLUSTER_NAME = abbr
  df
}

gep_tops = function(.obj, .title = NULL) {
  ggplot(.obj, aes(x = est_effect, y = lock_order(GENE_LABEL))) +
    facet_wrap( ~ CLUSTER_NAME, scales = "free", nrow = 2) +
    geom_errorbarh(aes(xmin = x_min_x, xmax = x_max_x), position=position_nudge(y = 0), height = 0, lwd = .2, col = "#555555") +
    {if(run == "v1")geom_point(aes(v1_g_value), position=position_nudge(y = 0), size = 1.2, col = "#44AA99")} +
    {if(run == "v2")geom_point(aes(v1_g_value), position=position_nudge(y = 0), size = 1.2, col = "#44AA99")} +
    {if(run == "v2")geom_point(aes(v2_g_value), position=position_nudge(y = 0), size = 1.2, col = "#117733")}+
    geom_point(aes(est_effect), position=position_nudge(y = 0), shape = 23, fill = NA, size = 1.2, col = "#555555") +
    mytheme(base_size = 8) +
    theme(
      legend.position = "none",
      axis.ticks.y  = element_blank(),
      plot.margin = margin(4, 4, 4, 4, "pt"),
      panel.spacing = unit(1, "lines"),
      plot.title = element_text(hjust = 0.5, face = "bold", colour = "black", size = rel(1))
    ) +
    xlab("Combined effect size") + ylab(NULL) +
    ggtitle(.title) +
    geom_vline(aes(xintercept = 0),color = "white", size=.1) +
    geom_hline(aes(yintercept = Inf), color = .obj$COL, size=2) +
    scale_x_continuous(breaks= pretty_breaks(4))
}

title.lab = mg.labels[names(mg.labels) %in% unique(ftrs$CLUSTER)]

i.gep = title.lab[1:4]
identity.gep = lapply(1:length(i.gep), function(x){ consensus_tops(i.gep[x], 15) })
identity.gep = do.call("rbind", identity.gep)
identity.gep$CLUSTER_NAME = factor(identity.gep$CLUSTER_NAME, levels = i.gep)

s.gep = title.lab[5:length(title.lab)]
shared.gep = lapply(1:length(s.gep), function(x){ consensus_tops(s.gep[x], 15) })
shared.gep = do.call("rbind", shared.gep)
shared.gep$CLUSTER_NAME = factor(shared.gep$CLUSTER_NAME, levels = s.gep)

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# FIGURE 3 and Supps (Discovery)
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
if(run == "d") {

  # Supps: Top ranked identify and shared GEPs
  i.s = plot_grid(gep_tops(identity.gep, "Identity Metagenes"), NULL, gep_tops(shared.gep, "Shared Metagenes"), nrow = 1, rel_widths = c(.3, .05, .5))
  ggsave2(filename="../figures/additional/fig_4_tops_shared_identity.pdf",
          plot = i.s,
          width = 410,
          height = 180,
          dpi = 100,
          bg = "white",
          units = "mm")

  b = plot_grid(
    exprs.1,
    exprs.2,
    exprs.3,
    exprs.4,
    nrow = 2
  )

  fig.a.b =
    plot_grid(
      NULL,
      sample.emb.pl + ggtitle("\n"),
      NULL,
      b,
      NULL,
      nrow = 1,
      rel_widths = c(.028, .52, .03, .5, .02),
      labels = c("  A", "", "  B"),
      label_fontface = "bold",
      label_size = 10
    )

  fig.d = plot_grid(
    kinetic.pl.identity + theme(legend.position = "none"),
    NULL,
    kinetic.pl.shared + theme(legend.position = "none"),
    NULL,
    nrow = 1,
    align = "h",
    rel_widths = c(.44, .05, .65, .02),
    labels = c("  C", "", ""), label_fontface = "bold", label_size = 10
  )

  fig.d.legend = get_legend(kinetic.pl.shared + theme(legend.position = "bottom") + guides(fill = guide_legend(nrow = 1)))
  fig.d.all = plot_grid(fig.d, fig.d.legend,  nrow = 2, rel_heights = c(.9, .1))

  d.e = plot_grid(
    gep_tops(identity.gep),
    NULL,
    ora,
    ncol = 3,
    # align = "h",
    rel_widths = c(.4, .05, .6),
    labels = c("  D", "", "    E"), label_fontface = "bold", label_size = 10, vjust  = 0.5
  )

  ggsave2(filename="../figures/main/fig_4.png",
          # plot = plot_grid(fig.a.b, NULL, fig.d.all, NULL, d.e, nrow = 5, rel_heights = c(.5, 0.03, .53, 0.02, .65)),
          plot = plot_grid(fig.a.b, NULL, fig.d.all, NULL, d.e, nrow = 5, rel_heights = c(0.29, 0.02, 0.31, 0.02, .38)),
          width = 297,
          height = 370,
          dpi = 300,
          bg = "white",
          units = "mm")

  ggsave2(filename="../figures/main/fig_4.pdf",
          # plot = plot_grid(fig.a.b, NULL, fig.d.all, NULL, d.e, nrow = 5, rel_heights = c(.5, 0.03, .53, 0.02, .65)),
          plot = plot_grid(fig.a.b, NULL, fig.d.all, NULL, d.e, nrow = 5, rel_heights = c(0.29, 0.02, 0.31, 0.02, .38)),
          width = 297,
          height = 370,
          dpi = 300,
          bg = "white",
          units = "mm")

  ggsave2(filename="../figures/additional/additional_file_2.png",
          plot = fig.d.all,
          width = 297,
          height = 110,
          dpi = 90,
          bg = "white",
          units = "mm")

  ggsave2(filename="../figures/additional/fig_4_reactome.pdf",
          plot = ora_plot(ora.obj.reactome, .go = F),
          width = 297,
          height = 120,
          dpi = 90,
          bg = "white",
          units = "mm")

}

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# FIGURE 3; Verification, Supps
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
if(run == "v1" | run == "v2") {

  a = plot_grid(
    kinetic.pl.identity + theme(legend.position = "none"),
    NULL,
    kinetic.pl.shared + theme(legend.position = "none"),
    NULL,
    nrow = 1,
    align = "h",
    rel_widths = c(.43, .05, .65, .02)
  )

  a.legend = get_legend(kinetic.pl.shared + theme(legend.position = "bottom") + guides(fill = guide_legend(nrow = 1)))
  fig.a.all = plot_grid(a, a.legend,  nrow = 2, rel_heights = c(.9, .1))

  i.gep = title.lab[1:4]
  identity.gep = lapply(1:length(i.gep), function(x){ consensus_tops(i.gep[x], 25) })
  identity.gep = do.call("rbind", identity.gep)
  identity.gep$CLUSTER_NAME = factor(identity.gep$CLUSTER_NAME, levels = i.gep)

  s.gep = title.lab[5:length(title.lab)]
  shared.gep = lapply(1:length(s.gep), function(x){ consensus_tops(s.gep[x], 25) })
  shared.gep = do.call("rbind", shared.gep)
  shared.gep$CLUSTER_NAME = factor(shared.gep$CLUSTER_NAME, levels = s.gep)


  b = plot_grid(
    gep_tops(identity.gep, "Identity Metagenes"), NULL, gep_tops(shared.gep, "Shared Metagenes"), NULL,
    nrow = 1, rel_widths = c(.4, .05, .6, .02))

  a.b.c = plot_grid(
    fig.a.all, NULL, b, NULL , ora,  nrow = 5, rel_heights = c(.6, 0.01, .7, 0.04, .5),
    labels = c("a", "", "b", "", "c"), label_fontface = "bold", label_size = 11, vjust  = 1.2
  )

  ggsave2(filename = paste0("../figures/additional/", supp.fig, ".pdf"),
          plot = a.b.c,
          width = 297,
          height = 370,
          dpi = 150,
          bg = "white",
          units = "mm")

  ggsave2(filename = paste0("../figures/additional/", supp.fig, "_reactome.pdf"),
          plot = ora_plot(ora.obj.reactome, .go = F),
          width = 297,
          height = 120,
          dpi = 90,
          bg = "white",
          units = "mm")
}

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# SUPPLEMENT: Concatenated coef Heatmap
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
pd = data.frame(colData(meta.prep))
pd = pd[rownames(pd) %in% colnames(h.merge), ]
pd$SAMPLE_ID = rownames(pd)
pd$HOURS_NUM = as.numeric(gsub("h", "", pd$HOURS))
pd = pd[naturalorder(pd$HOURS_NUM, pd$SUBTYPE, decreasing = F), ]

hm.features = c("HOURS", "SUBTYPE", "STUDY", "TISSUE")
anno.col = as.data.frame(pd[,hm.features, drop = F])
anno.col = anno.col %>%  mutate(
  SUBTYPE = case_when(
    SUBTYPE == "thp" ~ "unactivated",
    SUBTYPE == "th0" ~ "Th0",
    SUBTYPE == "th17" ~ "Th17",
    SUBTYPE == "itreg" ~ "iTreg",
    SUBTYPE == "th1" ~ "Th1",
    SUBTYPE == "th2" ~ "Th2"
  )
)
colnames(anno.col) = c("Time point (hours)", "T-cell population", "Study", "Sample source")

hours = c("#BBBBBB", tol9qualitative)[1:length(unique(anno.col[[1]]))]
names(hours) <- unique(anno.col[[1]])

phenotype = c("#4477AA", "#88CCEE", "#44AA99","#117733", "#999933", "#DDCC77")
names(phenotype) = c("unactivated", "Th0", "iTreg", "Th17", "Th2", "Th1")
unique(anno.col[[2]])

study = c("#77AADD", "#99DDFF", "#44BB99", "#BBCC33", "#AAAA00", "#EEDD88")
names(study) = unique(anno.col[[3]])

sample.source = c("#332288", "#88CCEE")
names(sample.source) = unique(anno.col[[4]])

anno.col.colors = list(
  "Time point (hours)" = hours,
  "T-cell population" = phenotype,
  "Study" = study,
  "Sample source" = sample.source
)

mat = h.merge[, rownames(pd)]
for (i in 1:ncol(mat)) {
  mat[ ,i] = mat[ ,i] / sum(mat[ ,i], na.rm = T)
}

ht.fonsize = 6
ht.anno.col = HeatmapAnnotation(
  df = anno.col,
  col = anno.col.colors,
  gap = unit(rep(1,length(anno.col)), "mm"),
  annotation_name_side = "left",
  annotation_legend_param = list(title_gp = gpar(fontsize = ht.fonsize),
                                 labels_gp = gpar(fontsize = ht.fonsize),
                                 ncol = 1),
  annotation_height = unit(rep(2,length(anno.col)), "mm"),
  simple_anno_size_adjust = TRUE,
  show_annotation_name = T,
  annotation_name_gp = gpar(fontsize = ht.fonsize),
  na_col = "black"
)

ht1 =
  ComplexHeatmap::Heatmap(
    mat, name = "mat",
    show_row_names = T,
    show_column_names = FALSE,
    column_names_gp = gpar(fontsize = ht.fonsize),
    row_names_gp = gpar(fontsize = ht.fonsize),
    col = scico(30, palette = "lajolla"),
    column_dend_gp = gpar(lwd = .3),
    row_dend_gp = gpar(lwd = .3),
    cluster_rows = F,
    cluster_columns = F,
    heatmap_legend_param = list(color_bar = "continuous", ncol = 1, direction = "horizontal",
                                grid_height  = unit(.1, "cm"),
                                legend_width = unit(2, "cm"),
                                title = NULL,
                                labels_gp = gpar(fontsize = ht.fonsize),
                                title_gp = gpar(fontsize = ht.fonsize)),
    top_annotation = ht.anno.col
  )

png("../figures/additional/nmf_coef_merged_heatmap.png", width = 210, height = 50, res = 300, units = "mm")
ComplexHeatmap::draw(ht1, heatmap_legend_side = "bottom",
                     annotation_legend_side = "right",
                     row_dend_side = "left")
dev.off()