grp_ave_sd = function(.obj, .ave = "mean") {
  grp.ave.l = list()
  grp.sd.l = list()
  for (i in paste0("^", unique(colnames(.obj)), "$")) {

    grp = .obj[, grepl(i , colnames(.obj))]

    grp.sd = data.frame(SD = matrixStats::rowSds(grp), row.names =  rownames(grp))
    colnames(grp.sd) = unique(colnames(grp))

    if (.ave == "mean") {
      grp.ave = data.frame(rowMeans(grp), row.names = rownames(grp))
    } else {
      grp.ave = data.frame(rowMedians(grp), row.names = rownames(grp))
    }
    colnames(grp.ave) = unique(colnames(grp))


    grp.ave.l[[unique(colnames(grp))]] = grp.ave
    grp.sd.l[[unique(colnames(grp))]] = grp.sd
  }
  grp.ave.df = do.call("cbind", grp.ave.l)
  grp.sd.df = do.call("cbind", grp.sd.l)
  return(list(ave = grp.ave.df, sd = grp.sd.df))
}

z_score = function(z) {
  rowmean = apply(z, 1, mean, na.rm=TRUE)
  rowsd = apply(z, 1, sd, na.rm=TRUE)
  rv = sweep(z, 1, rowmean,"-")
  rv = sweep(rv, 1, rowsd, "/")
  return(rv)
}

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Calc combined effect size
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
comb_effect = function(.obj = list(), .geneID = "ENSEMBL_ID",
                    .effectSize = "hedges.eb",  .se = "se.eb", .ncores = 1,
                    .confectsFDR = 0.05) {

  # subsetting obj
  obj.s = lapply(.obj, function(x) {
    dplyr::select(x, all_of(.geneID), all_of(.effectSize), all_of(.se))
  })

  # merging obj
  for (st in names(obj.s)) {
    colnames(obj.s[[st]]) = paste(colnames(obj.s[[st]]), st, sep = "_")
    colnames(obj.s[[st]])[grep(.geneID, colnames(obj.s[[st]]))] = .geneID
  }
  obj.merged = Reduce(function(x, y)
    merge(x = x, y = y, by = .geneID, all = TRUE), obj.s)

  # split by .geneID
  obj.merged.spl = split(obj.merged, obj.merged[[.geneID]])

  # salculating the combined effect size
  res.meta.l = mclapply(obj.merged.spl, function(gene) {
    meta_run(gene, .effectSize, .se)
  }, mc.cores = .ncores)

  res.meta.obj = lapply(res.meta.l, function(x) {x$obj})
  res.meta = lapply(res.meta.l, function(x) {x$df})

  res.meta = cbind(obj.merged, do.call("rbind", res.meta))

  # Topconfects ranking
  res.meta = res.meta %>%
    dplyr::mutate(est_se = (est_ci_r - est_ci_l) / 3.92) %>% # 95% CI
    dplyr::mutate(index = seq(nrow(res.meta)))

  confects = normal_confects(
    effect = res.meta$est_effect,
    se = res.meta$est_se,
    signed = TRUE,
    fdr = .confectsFDR,
    full = TRUE)

  res.meta = merge(
    res.meta,
    dplyr::select(confects$table, c(index, `rank`, confect, confect_fdr_zero = fdr_zero)),
    by = 'index', all = TRUE)
  res.meta$index = NULL

  res.meta = dplyr::arrange(res.meta, `rank`)
  return(list(df = res.meta, obj = res.meta.obj))
}

meta_run = function(gene, .effectSize, .se) {

  g = gene[which(!is.na(gene))]

  ef = as.numeric(dplyr::select(g, matches(.effectSize)))
  se = as.numeric(dplyr::select(g, matches(.se)))
  st = g[, which(grepl(paste0("^", .effectSize), colnames(g)))]
  st = gsub(paste0(.effectSize, "_"), "", colnames(st))

  # compute random effect model for the fold-changes and its variace
  run = tryCatch({ metafor::rma(yi = ef,  vi = se, slab = st, method="REML") },
                 error = function(e){ return(e) })

  # Increase iterations in case Fisher scoring algorithm doesn't converge
  if(any(is(run) == 'simpleError')) {
    run = tryCatch({metafor::rma(yi = ef,  vi = se, slab = st, method="REML",
                                 control = list(maxiter = 5000, stepadj = 0.5))},
                   error = function(e){ print(e); return(e) })
  }

  # If metafor is still returning error, give up and register line for gene
  if(any(is(run) == 'simpleError')) {
    df_res <- data.frame(signcon = length(which(ef>0))-length(which(ef<0)),
                         ntimes = length(ef != 0),
                         est_effect = NA,
                         est_ci_l = NA,
                         est_ci_r = NA,
                         est_pvalue = NA,
                         het_QE = NA,
                         het_QEp = NA,
                         het_QM = NA,
                         het_QMp = NA,
                         het_I2 = NA,
                         het_tau2 = NA,
                         error = TRUE)
  } else {
    df_res <- data.frame(signcon = length(which(ef>0))-length(which(ef<0)),
                         ntimes = length(ef != 0),
                         est_effect = as.numeric(run$beta),
                         est_ci_l = run$ci.lb,
                         est_ci_r = run$ci.ub,
                         est_pvalue = run$pval,
                         het_QE = run$QE,
                         het_QEp = run$QEp,
                         het_QM = run$QM,
                         het_QMp = run$QMp,
                         het_I2 = run$I2,
                         het_tau2 = run$tau2,
                         error = FALSE)
  }
  return(list(obj = run, df = df_res))
}

