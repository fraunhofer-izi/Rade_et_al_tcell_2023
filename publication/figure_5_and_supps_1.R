# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Libraries & Functions
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
.cran_packages = c("devtools","reshape2", "ggrepel", "stringr", "dplyr", "scales", "data.table",
                   "yaml", "cowplot", "scico", "ggtext", "circlize", "pdftools", "magick")
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

source("code/R/eda-plots.R")
source("paper-helper.R")
theme_set(mytheme())

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
manifest = yaml.load_file("assets/manifest.yaml")
work = manifest$workdata

ftrs.work = readRDS(paste0(work, manifest$pipeline_discovery))
no.input.nmf = nrow(ftrs.work$nmf.th0$se)

signtrs = readRDS(paste0(work, manifest$signatures_all))
signtrs$ENSEMBL_ID = gsub("\\..+", "", signtrs$ENSEMBL_ID)

v2.signtrs = subset(signtrs, v1_pass == TRUE & v2_pass == TRUE  & Housekeeping == TRUE)

cutoff.v2 = readRDS(paste0(work, manifest$timepoints_cutoff_v2))

se.izi = readRDS(paste0(work, manifest$limma_izi_thp_th0_act))$se

dge.act = readRDS(paste0(work, manifest$limma_izi_thp_th0_act))
dge.ctr = readRDS(paste0(work, manifest$limma_izi_ctrl))

dge.act.res = dge.act$dge.res
dge.ctr.res = dge.ctr$dge.res
rem = names(dge.ctr.res)[!grepl("^ACT_vs_CTR", names(dge.ctr.res))]
dge.ctr.res = dge.ctr.res[rem]

gc.anno = readRDS(paste0(work, manifest$gencode_v29_features))

# Results form metafor
obj.meta = readRDS(paste0(work, manifest$meta_all))
meta.res = obj.meta$meta.res

font.size = 8
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Heatmap
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
dge.res.act = data.frame(rbindlist(dge.act$dge.res, idcol = "CONTRAST"))
# dge.res.act = do.call("rbind", dge.act$dge.res)
dge.res.act = dge.res.act %>% .[!is.na(dge.res.act$confect), ]
dge.res.act.ftrs = unique(dge.res.act$ENSEMBL_ID_ABBR)

mean(table(dge.res.act$CONTRAST))

f = names(dge.ctr$dge.res)[!grepl("^ACT_vs_CTR", names(dge.ctr$dge.res))]
dge.res.ctrl = data.frame(rbindlist(dge.ctr$dge.res[f], idcol = "CONTRAST"))
# dge.res.ctrl = do.call("rbind", dge.ctr$dge.res[f])
dge.res.ctrl = dge.res.ctrl %>% .[!is.na(dge.res.ctrl$confect), ]
dge.res.ctrl.ftrs = unique(dge.res.ctrl$ENSEMBL_ID_ABBR)

mean(table(dge.res.ctrl$CONTRAST))

dge.shared = intersect(dge.res.ctrl.ftrs, dge.res.act.ftrs)
dge.shared = dge.shared[dge.shared %in% v2.signtrs$ENSEMBL_ID]
length(dge.shared)

mat = z_score(assays(se.izi)$cpm[dge.shared, ])

hm.features = c("GROUP_BASE", "HOURS")
anno.col = as.data.frame(colData(se.izi)[,hm.features, drop = F])
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
colnames(anno.col) = c("Condition", "Time point")

group = c("#000000", "#959595", "#009988")
names(group) = c("Unactivated (0h)", "Negative control", "Activated")

time.point = c("#BBBBBB", "#999933", "#DDCC77", "#CC6677", "#882255", "#AA4499")
names(time.point) = c("0h", "6h", "12h", "24h", "48h", "72h")

anno.col.colors = list(
  "Condition" = group,
  "Time point" = time.point
)

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
    use_raster = TRUE, raster_by_magick = TRUE,
    show_row_names = FALSE,
    show_column_names = FALSE,
    column_names_gp = gpar(fontsize = ht.fonsize),
    col = col_fun,
    column_dend_gp = gpar(lwd = .3),
    column_dend_height = unit(5, "mm"),
    row_dend_gp = gpar(lwd = .3),
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

hm = draw(ht1, merge_legend = TRUE)

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Boxplots
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
dge.res.act = do.call("rbind", dge.act$dge.res)
dge.res.act = dge.res.act[dge.res.act$ENSEMBL_ID_ABBR %in% dge.shared, ]
dge.res.act$sign = !is.na(dge.res.act$confect)
dge.res.act = dge.res.act %>% select(logFC, confect, ENSEMBL_ID_ABBR, sign) %>% mutate(Group = "Activated")

f = names(dge.ctr$dge.res)[!grepl("^ACT_vs_CTR", names(dge.ctr$dge.res))]
dge.res.ctrl = do.call("rbind", dge.ctr$dge.res[f])
dge.res.ctrl = dge.res.ctrl[dge.res.ctrl$ENSEMBL_ID_ABBR %in% dge.shared, ]
dge.res.ctrl$sign = !is.na(dge.res.ctrl$confect)
dge.res.ctrl = dge.res.ctrl %>% select(logFC, confect, ENSEMBL_ID_ABBR, sign) %>% mutate(Group = "Negative control")

df = rbind(dge.res.ctrl, dge.res.act)

df$HOURS = gsub("_vs_0h.+", "", rownames(df))
df$HOURS = gsub(".+\\_", "", df$HOURS)
df$HOURS = gsub("h", "", df$HOURS)
df$HOURS = factor(df$HOURS, levels = naturalsort(unique(df$HOURS)))
df$logFC_bool = df$logFC < 0

t.labs <- c("Log2 fold change < 0", "Log2 fold change > 0")
names(t.labs) <- c(TRUE, FALSE)

tmp = table(abs(subset(df, Group == "Negative control" & sign == T)$logFC) > 1)
tmp[2] / sum(tmp)

box.pl =
  # ggplot(df, aes(x = HOURS, y = ifelse(is.na(confect),0,confect), fill = Group)) +
  ggplot(df, aes(x = HOURS, y = logFC, fill = Group)) +
  geom_boxplot(fatten=1.3, lwd = .3, outlier.size = .1) +
  guides(fill = guide_legend(override.aes = list(size = .3))) +
  ylab("Log2 fold change") + xlab("Contrast (6-72 hours vs. 0 hours)") + labs(fill = NULL) +
  scale_fill_manual(values = c("#009988", "#959595")) +
  facet_wrap(~ logFC_bool, scales = "free", labeller = labeller(logFC_bool = t.labs)) +
  mytheme() +
  theme(
    legend.position = "right",
    panel.spacing = unit(.5, "lines"),
    axis.ticks.x=element_blank()
  )

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Venn as barplots (act vs 0 and ctrl vs 0)
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
dge.act.res.cp = dge.act.res
dge.ctr.res.cp = dge.ctr.res
names(dge.act.res.cp) = gsub("HOURS_|CTR_", "", names(dge.act.res.cp))
names(dge.ctr.res.cp) = gsub("HOURS_|CTR_", "", names(dge.ctr.res.cp))

l = list()
for (i in names(dge.act.res.cp)) {

  act = dge.act.res.cp[[i]]
  act = act[!is.na(act$confect), ]
  ctrl = dge.ctr.res.cp[[i]]
  ctrl = ctrl[!is.na(ctrl$confect), ]

  a = act[setdiff(act$ENSEMBL_ID_ABBR, ctrl$ENSEMBL_ID_ABBR), ]
  a.up = table(a$logFC > 0)[2]
  a.dwn = table(a$logFC > 0)[1]

  c = ctrl[setdiff(ctrl$ENSEMBL_ID_ABBR, act$ENSEMBL_ID_ABBR), ]
  c.up = table(c$logFC > 0)[2]
  c.dwn = table(c$logFC > 0)[1]

  act.up = act[act$logFC > 0, ]$ENSEMBL_ID_ABBR
  act.dwn = act[act$logFC < 0, ]$ENSEMBL_ID_ABBR

  ctrl.up = ctrl[ctrl$logFC > 0, ]$ENSEMBL_ID_ABBR
  ctrl.dwn = ctrl[ctrl$logFC < 0, ]$ENSEMBL_ID_ABBR

  i.up = length(intersect(act.up, ctrl.up))
  i.diff = length(intersect(act.dwn, ctrl.up)) + length(intersect(act.up, ctrl.dwn))
  i.dwn = length(intersect(act.dwn, ctrl.dwn))

  l[[i]] = data.frame(
    act = c(a.up, NA, a.dwn),
    ctrl = c(c.up, NA, c.dwn),
    act..ctrl = c(i.up, i.diff, i.dwn),
    row.names = c("up", "diff", "down")
  )
}
df.v = do.call("rbind", l)
df.v = reshape2::melt(as.matrix(df.v))
df.v$dir = gsub(".+\\.", "", df.v$Var1)
df.v$dir = factor(df.v$dir, levels = c("up", "down", "diff"))
df.v$ctrst = gsub("\\..+", "", df.v$Var1)
df.v$ctrst_label = gsub("\\_.+", "", df.v$ctrst)
df.v$ctrst_label = gsub("h", "", df.v$ctrst_label)
df.v$ctrst_label = factor(df.v$ctrst_label, levels = naturalsort(unique(df.v$ctrst_label)))
df.v$fill = paste0(df.v$Var2, "_", df.v$dir)
lvls = c("act_up", "act..ctrl_up", "ctrl_up", "act_down", "act..ctrl_down", "ctrl_down", "act..ctrl_diff")
df.v$fill = factor(df.v$fill, levels = lvls)

# fill.col = c("#B2182B", "#D6604D", "#F4A582", "#2166AC", "#4393C3", "#92C5DE", "#BBBBBB")
fill.col = c("#532C01", "#8D5003", "#C08229", "#2166AC", "#4393C3", "#92C5DE", "#BBBBBB")
names(fill.col) = lvls

leg.label = c(
  "up:\nact",
  # expression(paste("up: act",intersect("neg.ctrl"))),
  "up:\nact & neg. ctrl",
  "up:\nneg. ctrl",
  "down:\nact",
  # expression(paste("down: act",intersect("neg.ctrl"))),
  "down:\nact & neg. ctrl",
  "down:\nneg. ctrl",
  # expression(paste("up/down: act",intersect("neg.ctrl")))
  "diff. LFC:\nact & neg. ctrl"
)
names(leg.label) = lvls

mean(c(
  subset(df.v, fill == "act..ctrl_up")$value,
  subset(df.v, fill == "act..ctrl_down")$value
))

mean(subset(df.v, fill == "act..ctrl_diff")$value)

v.pl =
  ggplot(df.v, aes(x = ctrst_label, y = value, fill = fill)) +
  geom_bar(position="stack", width = .8, stat="identity") +
  scale_fill_manual(values = fill.col, labels = leg.label, name = expression("FDR"<0.05)) +
  guides(fill = guide_legend(label.hjust = 0, label.position = "right", override.aes = list(size=3))) +
  facet_wrap(~dir, ncol = 3) +
  mytheme(base_size = 8) +
  theme(
    aspect.ratio = 1.618,
    panel.spacing = unit(.3, "lines"),
    plot.margin = margin(20, 5, 5, 5, "pt"),
    legend.position = "right",
    legend.text = element_text(margin = margin(b = 2, t = 2)),
    axis.ticks.x = element_blank(),
    strip.text.x = element_blank(),
    plot.title = element_text(hjust = 0.5, face = "bold", colour = "black", size = rel(1))
  ) +
  ylab("Frequency") + xlab("Contrast (6-72 hours vs. 0 hours)")

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Expression profiles (activated and control)
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
ftrs = subset(signtrs, v1_pass == TRUE & v2_pass == TRUE  & Housekeeping == TRUE)
ftrs.fin = subset(signtrs, v1_pass == TRUE & v2_pass == TRUE & v2_ctr_pass == TRUE  & Housekeeping == TRUE)

exprs.cpm =z_score(cpm(assays(se.izi)$raw, log = F))
exprs.cpm = exprs.cpm[ftrs$ENSEMBL_ID, ]
colnames(exprs.cpm) = se.izi$GROUP

a = exprs.cpm[, grepl("naive_0h", colnames(exprs.cpm))]
colnames(a) = rep("control_0h", 4)
b = exprs.cpm[, grepl("naive_0h", colnames(exprs.cpm))]
colnames(b) = rep("activated_0h", 4)
exprs.cpm = cbind(a, b, exprs.cpm[, !grepl("naive_0h", colnames(exprs.cpm))])

exprs.ave = aveByGrp(exprs.cpm)
df = reshape2::melt(as.matrix(exprs.ave))

colnames(df) = c("ENSEMBL_ID", "GROUP", "EXPRS")
df$HOURS = gsub(".*\\_", "", df$GROUP)
df$PHENOTYPE = gsub("\\_.+", "", df$GROUP)
df$CLUSTER = ftrs$CLUSTER[match(df$ENSEMBL_ID, ftrs$ENSEMBL_ID)]
df$HOURS = gsub("h", "", df$HOURS)
lvls = c("0", "6", "12", "24", "48", "72")
df$HOURS = factor(df$HOURS, levels = lvls)
df$ID = paste0(df$ENSEMBL_ID,"_",df$PHENOTYPE)

title.lab = mg.labels[names(mg.labels) %in% unique(ftrs$CLUSTER)]
df$LABEL = title.lab[match(df$CLUSTER, names(title.lab))]
df$LABEL = factor(df$LABEL, levels = title.lab)
df$COL = mg.col[match(df$CLUSTER, names(mg.col))]

df$CLUSTER_HOURS = paste0(df$CLUSTER,"_",df$HOURS)

cutoff.v2["0100"] = NULL
cutoff.v2 = lapply(cutoff.v2, function(x){
  x[1]
})

lfc.thres = paste0(names(unlist(cutoff.v2)), "_", unlist(cutoff.v2))
# lfc.thres = gsub(".\\_", "\\_" , lfc.thres)
lfc.thres = gsub("h$", "", lfc.thres)

df = df %>% mutate(
  LFC_THRES = case_when(
    CLUSTER_HOURS %in% lfc.thres ~ "yes",
    TRUE ~ "no"
  )
)

df = df %>%
  mutate(x.label = paste("<span style = 'color: ",
                         ifelse(LFC_THRES == "yes", "#CC3311", "black"),
                         ";'>",
                         HOURS,
                         "</span>", sep = ""))

df$x.label = factor(df$x.label, levels = unique(df$x.label))

pl_f = function(.obj, .title = NULL) {
  ggplot(.obj, aes(x = x.label, y = EXPRS, group = ID, colour = PHENOTYPE)) +
    geom_line(lwd = .5, alpha = .7) +
    xlab("Hours after activation") +
    ylab("Standardized expression") +
    mytheme(base_size = font.size) +
    theme(
      panel.spacing = unit(.7, "lines"),
      legend.position = "bottom",
      # axis.ticks.x = element_blank(),
      axis.text.x = element_markdown(),
      legend.box.margin=margin(t=-5, b=-5),
      plot.title = element_text(hjust = 0.5, face = "bold", colour = "black", size = rel(1))
    ) +
    scale_x_discrete(expand = c(0, 0)) +
    guides(color = guide_legend(title = NULL, override.aes = list(size = 1, alpha = 1))) +
    scale_colour_manual(values = c("#009988", "#959595"), labels = c("Activated","Negative control")) +
    facet_wrap(~ LABEL, nrow = 2, scales = "free") +
    # stat_summary(data = obj[obj$PHENOTYPE == "control", ], fun=median, colour="#3b3b3b", geom="line", lwd = 1, aes(group = LABEL))  +
    # stat_summary(data = obj[obj$PHENOTYPE == "activated", ], fun=median, colour="#3b3b3b", geom="line", lwd = 1, aes(group = LABEL)) +
    geom_hline(aes(yintercept = Inf), color = factor(.obj$COL), size=3) +
    ggtitle(.title)
}



df.pass = df %>% .[df$ENSEMBL_ID %in% ftrs.fin$ENSEMBL_ID, ]

sel = (as.data.frame.matrix(table(df.pass$HOURS, df.pass$CLUSTER))[1, ] / 2) >=5

p1 = pl_f(subset(df.pass, CLUSTER != "0100" & CLUSTER %in% colnames(sel)[sel]), NULL)
# p1 = pl_f(subset(df.pass, CLUSTER == "1000" | CLUSTER == "0010" | CLUSTER == "0001" | CLUSTER == "0011"), NULL)
p1_i = pl_f(subset(df.pass, CLUSTER == "1000" | CLUSTER == "0010" | CLUSTER == "0001"), "Identity Metagenes")
p1_s = pl_f(subset(df.pass, CLUSTER != "1000" & CLUSTER != "0100" & CLUSTER != "0010" & CLUSTER != "0001"), "Shared Metagenes")

df.fail = df %>% .[!df$ENSEMBL_ID %in% ftrs.fin$ENSEMBL_ID, ]
p1_i.fail = pl_f(subset(df.fail, CLUSTER == "1000" | CLUSTER == "0100" | CLUSTER == "0010" | CLUSTER == "0001"), "Identity Metagenes")
p1_s.fail = pl_f(subset(df.fail, CLUSTER != "1000" & CLUSTER != "0100" & CLUSTER != "0010" & CLUSTER != "0001"), "Shared Metagenes")

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Piechart with fract
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
count.data <- data.frame(
  pheno = c("CD4+ naive", "CD4+ memory", "other CD4+ T-cells", "CD8+ naive", "CD8+ memory", "other CD8+ T-cells", "non T-cells"),
  class = c("38.89% CD4+ naive", "34.52% CD4+ memory", "2.61% other CD4+ T-cells", "11.47% CD8+ naive", "6.88% CD8+ memory", "3.34% other CD8+ T-cells", "2.29% non T-cells"),
  prop = c(38.89, 34.52, 2.61, 11.47, 6.88, 3.34, 2.29)
)

count.data$class = factor(count.data$class, levels = rev(naturalsort(count.data$class)))
count.data$pheno = factor(count.data$pheno, levels = rev(count.data[order(count.data$prop, decreasing = T), ]$pheno))

celltypes.pl =
  ggplot(count.data, aes(x = pheno, y = prop)) +
  geom_bar(width = .8, stat = "identity", color = "white", lwd = 0, fill = "#555555") +
  coord_flip() +
  geom_text(aes(label=prop), hjust= -0.2, size = 2.7) +
  ylim(0, 60) +
  mytheme(base_size = font.size) +
  theme(
    aspect.ratio = 1.618,
    panel.grid.minor.y = element_blank(),
    panel.grid.major.y = element_blank(),
    legend.position = "none",
    axis.ticks.y = element_blank(),
    legend.title = element_text(colour="red", size=8, face="bold")
  ) +
  labs(fill = "FAKE (Daten kommen\nnoch von Sebastian)") +
  guides(fill = guide_legend(label.hjust = 0, label.position = "right", override.aes = list(size = 1))) +
  ylab("Percentage") + xlab(NULL)

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
v1 = subset(signtrs, v1_pass == TRUE & Housekeeping == TRUE)
v1v2 = subset(signtrs, v1_pass == TRUE & v2_pass == TRUE & Housekeeping == TRUE)
v1v2ctr = subset(signtrs, v1_pass == TRUE & v2_pass == TRUE & v2_ctr_pass == TRUE & Housekeeping == TRUE)
v1v2ctr.lymph = subset(signtrs, v1_pass == TRUE & v2_pass == TRUE & v2_ctr_pass == TRUE & Lymphoid_Genes == TRUE & Housekeeping == TRUE)
v1v2ctr.tcell = subset(signtrs, v1_pass == TRUE & v2_pass == TRUE & v2_ctr_pass == TRUE & T_cell_Genes == TRUE & Housekeeping == TRUE)

fc =
  data.frame(x= 1:100, y= 1:100) %>%
  ggplot(aes(x, y)) +
  # ylim(84, 100) +
  scale_x_continuous(expand = c(0, 0),  limits = c(0, 55)) +
  scale_y_continuous(expand = c(0, 0),  limits = c(94, 101.7)) +
  # scale_y_continuous(expand = c(80, 100)) +
  theme_linedraw() +

  geom_rect(xmin = 0, xmax=10, ymin=94, ymax=100, color='black', fill='white', size=0.25) +
  annotate('text', x= 5, y=97,label= paste0(nrow(signtrs), " genes passed\nthe Discovery Set", "\n(of ", no.input.nmf ," genes\nused for NMF)"), size=2.8) +
  annotate(
    'text', x= 12.5, y=101.3,
    label= paste0(nrow(subset(signtrs, Housekeeping == FALSE)), " Housekeeping genes removed"), size=2.8) +

  geom_rect(xmin = 15, xmax=25, ymin=94, ymax=100, color='black', fill='white', size=0.25) +
  annotate('text', x= 20, y=97,label= paste0(nrow(v1), " genes verified by\nthe Memory T-cell\nVerification Set"), size=2.8) +

  geom_rect(xmin = 30, xmax=40, ymin=94, ymax=100, color='black', fill='white', size=0.25) +
  annotate('text', x= 35, y=97,label= paste0(nrow(v1v2), " genes verified by\nthe Pan T-cell\nVerification Set\n(Activated)"), size=2.8) +

  geom_rect(xmin = 45, xmax=55, ymin=94, ymax=100, color='black', fill='white', size=0.25) +
  annotate('text', x= 50, y=97,label= paste0(nrow(v1v2ctr), " genes verified by\nthe Pan T-cell\nVerification Set\n(Negative control)"), size=2.8) +

  geom_segment(
    x=10, xend=15, y=97, yend=97, size=0.15, linejoin = "mitre", lineend = "butt", arrow = arrow(length = unit(1, "mm"), type= "closed")
  ) +
  geom_segment(
      x=12.5, xend=12.5, y=97, yend=100.2, size=0.15, linejoin = "mitre", lineend = "butt", arrow = arrow(length = unit(1, "mm"), type= "closed")
  ) +
  geom_segment(
    x=25, xend=30, y=97, yend=97, size=0.15, linejoin = "mitre", lineend = "butt", arrow = arrow(length = unit(1, "mm"), type= "closed")
  ) +
  geom_segment(
    x=40, xend=45, y=97, yend=97, size=0.15, linejoin = "mitre", lineend = "butt", arrow = arrow(length = unit(1, "mm"), type= "closed")
  ) +
  mytheme(base_size = font.size) +
  theme(
    axis.line=element_blank(),
    axis.text=element_blank(),
    axis.ticks=element_blank(),
    axis.title=element_blank(),
    panel.border=element_blank()
  )


# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
wf = ggdraw() + draw_image(magick::image_read_pdf("../figures/main/izi_seq.pdf", density = 600))
wf = wf + theme(plot.margin = margin(0, 0, 0, 0, "cm"))

a.b = plot_grid(
  wf,
  NULL,
  celltypes.pl + theme(plot.margin = margin(20, 4, 20, 4, "pt")),
  NULL,
  v.pl + theme(plot.margin = margin(20, 4, 20, 4, "pt")),
  NULL,
  rel_widths = c(.35, .05, .22, .02, .46, .01),
  nrow = 1, labels = c("A", "", "B", "", "C"), label_fontface = "bold", label_size = 10
)

c.d.e = plot_grid(
  plot_grid(
    grid.grabExpr(draw(ht1, merge_legend = TRUE)),
    NULL,
    box.pl + theme(plot.margin = margin(5, 0, 5, 5, "pt")),
    nrow = 3, rel_heights = c(.96, .01, .4)
  ),
  NULL,
  plot_grid(
    p1 + theme(legend.position = "bottom"),
    fc + theme(plot.margin = margin(10, 4, 4, 26, "pt")),
    rel_heights = c(.8, .2),
    nrow = 2,
    labels = c("E", "F"), label_fontface = "bold", label_size = 10
  ),
  NULL,
  rel_widths = c(.379, 0.05, .617, 0.01), nrow = 1, labels = c("D", "", "", ""), label_fontface = "bold", label_size = 10
)

ggsave2(filename="../figures/main/fig_5.pdf",
        plot = plot_grid(a.b, c.d.e, rel_heights = c(.25, .55), nrow = 2),
        width = 297,
        height = 200,
        dpi = 300,
        bg = "white",
        units = "mm")

ggsave2(filename="../figures/main/fig_5.png",
        plot = plot_grid(a.b, c.d.e, rel_heights = c(.25, .55), nrow = 2),
        width = 297,
        height = 200,
        # type = "cairo",
        dpi = 300,
        bg = "white",
        units = "mm")

ggsave2(filename="../figures/additional/v2_passed_shared_identity.pdf",
        plot = plot_grid(p1_i, NULL, p1_s, nrow = 1, rel_widths = c(.4, .05, .6)),
        width = 250,
        height = 110,
        dpi = 100,
        bg = "white",
        units = "mm")

ggsave2(filename="../figures/additional/v2_passed_shared_identity_fail.pdf",
        plot = plot_grid(p1_i.fail, NULL, p1_s.fail, NULL, nrow = 1, rel_widths = c(.5, .05, .5, .25)),
        width = 250,
        height = 110,
        dpi = 100,
        bg = "white",
        units = "mm")

