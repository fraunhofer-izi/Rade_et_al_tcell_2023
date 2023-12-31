library(igraph)
library(ggraph)
library(scatterpie)
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Emapplot
# Alle nötigen Funktion von "https://github.com/YuLab-SMU/enrichplot/" kopiert
# Grund: Um die font size der node Labels zu ändern
# Änderung: geom_node_text() Zeile 260
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>


##' Get graph.data.frame() result
##'
##' @param y a data.frame of clusterProfiler result
##' @param geneSets a list gene sets with the names of enrichment IDs
##' @param color a string, the column name of y for nodes colours
##' @param line_scale scale of line width
##' @return result of graph.data.frame()
##' @noRd
emap_graph_build <- function(y,geneSets,color,line_scale) {
  if (is.null(dim(y)) | nrow(y) == 1) {
    g <- graph.empty(0, directed=FALSE)
    g <- add_vertices(g, nv = 1)
    V(g)$name <- as.character(y$Description)
    V(g)$color <- "red"
    ##return(ggraph(g))
  } else {
    id <- y[,"ID"]
    geneSets <- geneSets[id]
    n <- nrow(y) #
    w <- matrix(NA, nrow=n, ncol=n)
    colnames(w) <- rownames(w) <- y$Description

    for (i in seq_len(n-1)) {
      for (j in (i+1):n) {
        w[i,j] <- overlap_ratio(geneSets[id[i]], geneSets[id[j]])
      }
    }

    wd <- melt(w)
    wd <- wd[wd[,1] != wd[,2],]
    wd <- wd[!is.na(wd[,3]),]
    g <- graph.data.frame(wd[,-3], directed=FALSE)
    E(g)$width=sqrt(wd[,3] * 5) * line_scale
    g <- delete.edges(g, E(g)[wd[,3] < 0.2])
    ## g <- delete.edges(g, E(g)[wd[,3] < 0.05])
    idx <- unlist(sapply(V(g)$name, function(x) which(x == y$Description)))

    cnt <- sapply(geneSets[idx], length)
    V(g)$size <- cnt

    colVar <- y[idx, color]
    V(g)$color <- colVar
  }

  return(g)
}




##' @rdname emapplot
##' @importFrom igraph graph.empty
##' @importFrom igraph add_vertices
##' @importFrom igraph graph.data.frame
##' @importFrom igraph delete.edges
##' @importFrom igraph V "V<-"
##' @importFrom igraph E "E<-"
##' @importFrom reshape2 melt
##' @importFrom ggplot2 aes_
##' @importFrom ggplot2 scale_color_gradientn
##' @importFrom ggplot2 guide_colorbar
##' @importFrom ggplot2 scale_size
##' @importFrom ggplot2 theme_void
##' @importFrom ggraph ggraph
##' @importFrom ggraph geom_node_point
##' @importFrom ggraph geom_node_text
##' @importFrom ggraph geom_edge_link
##' @importFrom DOSE geneInCategory
##' @param pie_scale scale of pie plot
##' @param line_scale scale of line width
##' @author Guangchuang Yu
emapplot.enrichResult <- function(x, showCategory = 30, color="p.adjust", layout = "nicely",
                                  pie_scale = 1, line_scale = 1, ...) {
  n <- update_n(x, showCategory)
  geneSets <- geneInCategory(x) ## use core gene for gsea result
  y <- as.data.frame(x)
  if (is.numeric(n)) {
    y <- y[1:n,]
  } else {
    y <- y[match(n, y$Description),]
    n <- length(n)
  }


  if (n == 0) {
    stop("no enriched term found...")
  }

  g <- emap_graph_build(y=y,geneSets=geneSets,color=color, line_scale=line_scale)
  if(n == 1) {
    return(ggraph(g) + geom_node_point(color="red", size=5) + geom_node_text(aes_(label=~name)))
  }
  ##} else {
  p <- ggraph(g, layout=layout)
  if (length(E(g)$width) > 0) {
    p <- p + geom_edge_link(alpha=.8, aes_(width=~I(width)), colour='darkgrey')
  }
  p + geom_node_point(aes_(color=~color, size=~size)) +
    geom_node_text(aes_(label=~name), repel=TRUE) + theme_void() +
    scale_color_continuous(low="red", high="blue", name = color, guide=guide_colorbar(reverse=TRUE)) +
    scale_size(range=c(3, 8) * pie_scale)

}



##' Merge the compareClusterResult file
##'
##' @param yy a data.frame of clusterProfiler result
##'
##' @return a data.frame
##' @noRd
merge_compareClusterResult <- function(yy) {
  yy_union<- yy[!duplicated(yy$ID),]
  yy_ids <- lapply(split(yy, yy$ID), function(x) {
    ids <- unique(unlist(strsplit(x$geneID, "/")))
    cnt <- length(ids)
    list(ID=paste0(ids, collapse="/"), cnt=cnt)
  })

  ids <- vapply(yy_ids, function(x) x$ID, character(1))
  cnt <- vapply(yy_ids, function(x) x$cnt, numeric(1))

  yy_union$geneID = ids[yy_union$ID]
  yy_union$Count = cnt[yy_union$ID]
  yy_union$Cluster = NULL
  yy_union
}


##' @importFrom igraph graph.empty
##' @importFrom igraph add_vertices
##' @importFrom igraph graph.data.frame
##' @importFrom igraph delete.edges
##' @importFrom igraph V "V<-"
##' @importFrom igraph E "E<-"
##' @importFrom reshape2 melt
##' @importFrom ggplot2 aes_
##' @importFrom ggplot2 scale_color_gradientn
##' @importFrom ggplot2 guide_colorbar
##' @importFrom ggplot2 scale_size
##' @importFrom ggplot2 theme_void
##' @importFrom ggplot2 geom_text
##' @importFrom ggplot2 coord_equal
##' @importFrom ggplot2 labs
##' @importFrom ggplot2 xlim
##' @importFrom ggplot2 ylim
##' @importFrom ggraph ggraph
##' @importFrom ggraph geom_node_point
##' @importFrom ggraph geom_node_text
##' @importFrom ggraph geom_edge_link
##' @importFrom scatterpie geom_scatterpie
##' @importFrom scatterpie geom_scatterpie_legend
##' @importFrom DOSE geneInCategory
##' @importFrom DOSE geneID
##' @importClassesFrom DOSE compareClusterResult
##' @param split separate result by 'category' variable
##' @param pie proportion of clusters in the pie chart, one of 'equal' (default) or 'Count'
##' @param legend_n number of circle in legend
##' @param pie_scale scale of pie plot
##' @param line_scale scale of line width
##' @method fortify compareClusterResult
##' @importFrom scatterpie geom_scatterpie
##' @importFrom stats setNames
##' @noRd
emapplot_mod <- function(x, showCategory = 5, color = "p.adjust", node.font.size = 6,
                                          layout = "nicely", split=NULL, pie = "equal",
                                          legend_n = 5, pie_scale = 1, line_scale = 1, ...) {

  region <- radius <- NULL

  ## pretreatment of x, just like dotplot do
  y <- fortify.compareClusterResult(x, showCategory=showCategory,
                                    includeAll=TRUE, split=split)

  y$Cluster = sub("\n.*", "", y$Cluster)


  ## geneSets <- geneInCategory(x) ## use core gene for gsea result

  ## Data structure transformation, combining the same ID (Description) genes
  n <- update_n(x, showCategory)

  y_union <- merge_compareClusterResult(y)

  if (n == 0) {
    stop("no enriched term found...")
  }
  geneSets <- setNames(strsplit(as.character(y_union$geneID), "/", fixed = TRUE), y_union$ID)
  g <- emap_graph_build(y=y_union,geneSets=geneSets,color=color, line_scale=line_scale)
  ## when y just have one line
  if(is.null(dim(y)) | nrow(y) == 1) {
    title <- y$Cluster
    p <- ggraph(g) + geom_node_point(color="red", size=5 * pie_scale) +
      geom_node_text(aes_(label=~name)) + theme_void() +
      labs(title=title)
    return(p)
  }

  if(is.null(dim(y_union)) | nrow(y_union) == 1) {
    ##return(ggraph(g) + geom_node_point(color="red", size=5) + geom_node_text(aes_(label=~name)))
    p <- ggraph(g)
    ID_Cluster_mat <- prepare_pie_category(y, pie=pie)

    ID_Cluster_mat <- cbind(ID_Cluster_mat,1,1,0.1*pie_scale)
    colnames(ID_Cluster_mat) <- c(colnames(ID_Cluster_mat)[1:(ncol(ID_Cluster_mat)-3)],"x","y","radius")


    p <- p + geom_scatterpie(aes_(x=~x,y=~y,r=~radius), data=ID_Cluster_mat,
                             cols=names(ID_Cluster_mat)[1:(ncol(ID_Cluster_mat)-3)],color=NA)+
      xlim(-3,3) + ylim(-3,3) + coord_equal()+ geom_node_text(aes_(label=~name), repel=TRUE) +
      theme_void()+labs(fill = "Cluster")
    return(p)

  }
  p <- ggraph(g, layout=layout)
  if (length(E(g)$width) > 0) {
    p <- p + geom_edge_link(alpha=.8, aes_(width=~I(width)), colour='darkgrey')
  }

  ## then add the pie plot
  ## Get the matrix data for the pie plot
  ID_Cluster_mat <- prepare_pie_category(y,pie=pie)


  #plot the edge
  #get the X-coordinate and y-coordinate of pies
  aa <- p$data

  desc <- y_union$Description[match(rownames(ID_Cluster_mat), y_union$Description)]
  i <- match(desc, aa$name)

  ID_Cluster_mat$x <- aa$x[i]
  ID_Cluster_mat$y <- aa$y[i]

  #Change the radius value to fit the pie plot
  ID_Cluster_mat$radius <- sqrt(aa$size[i] / sum(aa$size)) * pie_scale
  #ID_Cluster_mat$radius <- sqrt(aa$size / pi)

  x_loc1 <- min(ID_Cluster_mat$x)
  y_loc1 <- min(ID_Cluster_mat$y)
  ## x_loc2 <- min(ID_Cluster_mat$x)
  ## y_loc2 <- min(ID_Cluster_mat$y)+0.1*(max(ID_Cluster_mat$y)-min(ID_Cluster_mat$y))
  if(ncol(ID_Cluster_mat) > 4) {
    p <- p + geom_scatterpie(aes_(x=~x,y=~y,r=~radius), data=ID_Cluster_mat,
                             cols=colnames(ID_Cluster_mat)[1:(ncol(ID_Cluster_mat)-3)],color=NA) +

      coord_equal()+
      geom_node_text(aes_(label=~name), size = node.font.size, repel=TRUE) + theme_void() +
      geom_scatterpie_legend(ID_Cluster_mat$radius, x=x_loc1, y=y_loc1, n = legend_n,
                             labeller=function(x) {round(sum(aa$size)*((x/pie_scale)^2))}) +
      labs(fill = "Cluster")
    return(p)
  }
  ## annotate("text", label = "gene number", x = x_loc2, y = y_loc2, size = 4, colour = "red")
  title <- colnames(ID_Cluster_mat)[1]
  p + geom_node_point(aes_(color=~color, size=~size)) +
    geom_node_text(aes_(label=~name), repel=TRUE) + theme_void() +
    scale_color_continuous(low="red", high="blue", name = color, guide=guide_colorbar(reverse=TRUE)) +
    scale_size(range=c(3, 8) * pie_scale)  + labs(title= title)
}

##' Prepare pie data for genes in cnetplot.
##' The function only works for compareClusterResult
##'
##' @param y a data.frame converted from compareClusterResult
##' @return a data.frame
##' @noRd
prepare_pie_gene <- function(y) {
  gene_pie <- tibble::as_tibble(y[,c("Cluster", "Description", "geneID")])
  gene_pie$geneID <- strsplit(gene_pie$geneID, '/')
  gene_pie2 <- as.data.frame(tidyr::unnest(gene_pie, cols=geneID))
  gene_pie2 <- unique(gene_pie2)
  prepare_pie_data(gene_pie2, pie =  "equal", type = "gene")
}


##' Prepare pie data for categories in cnetplot/emapplot.
##' The function only works for compareClusterResult
##'
##' @param y a data.frame converted from compareClusterResult
##' @param pie proportion of clusters in the pie chart, one of 'equal' (default) or 'Count'
##' @return a data.frame
##' @noRd
prepare_pie_category <- function(y, pie = "equal") {
  pie <- match.arg(pie, c("equal", "count", "Count"))
  if (pie == "count") pie <- "Count"

  pie_data <- y[,c("Cluster", "Description", "Count")]
  pie_data[,"Description"] <- as.character(pie_data[,"Description"])
  prepare_pie_data(pie_data, pie = pie)
}




prepare_pie_data <- function(pie_data, pie = "equal",type = "category") {
  if(type == "category"){
    ID_unique <- unique(pie_data[,2])
  } else {
    ID_unique <- unique(pie_data[,3])
  }

  Cluster_unique <- unique(pie_data[,1])
  ID_Cluster_mat <- matrix(0, nrow = length(ID_unique), ncol = length(Cluster_unique))
  rownames(ID_Cluster_mat) <- ID_unique
  colnames(ID_Cluster_mat) <- Cluster_unique
  ID_Cluster_mat <- as.data.frame(ID_Cluster_mat, stringAsFactors = FALSE)
  if(pie == "Count") {
    for(i in seq_len(nrow(pie_data))) {
      ID_Cluster_mat[pie_data[i,2],pie_data[i,1]] <- pie_data[i,3]
    }
    for(kk in seq_len(ncol(ID_Cluster_mat))) {
      ID_Cluster_mat[,kk] <- as.numeric(ID_Cluster_mat[,kk])
    }
    return(ID_Cluster_mat)
  }
  for(i in seq_len(nrow(pie_data))) {
    if(type == "category"){
      ID_Cluster_mat[pie_data[i,2],pie_data[i,1]] <- 1
    } else {
      ID_Cluster_mat[pie_data[i,3],pie_data[i,1]] <- 1
    }

  }
  return(ID_Cluster_mat)
}


##' create color palette for continuous data
##'
##'
##' @title color_palette
##' @param colors colors of length >=2
##' @return color vector
##' @importFrom grDevices colorRampPalette
##' @export
##' @examples
##' color_palette(c("red", "yellow", "green"))
##' @author guangchuang yu
color_palette <- function(colors) colorRampPalette(colors)(n = 299)

sig_palette <- color_palette(c("red", "yellow", "blue"))

heatmap_palette <- color_palette(c("red", "yellow", "green"))

overlap_ratio <- function(x, y) {
  x <- unlist(x)
  y <- unlist(y)
  length(intersect(x, y))/length(unique(c(x,y)))
}

fc_readable <- function(x, foldChange = NULL) {
  if (is.null(foldChange))
    return(NULL)

  if(x@readable) {
    gid <- names(foldChange)
    if (is(x, 'gseaResult')) {
      ii <- gid %in% names(x@geneList)
    } else {
      ii <- gid %in% x@gene
    }
    gid[ii] <- x@gene2Symbol[gid[ii]]
    names(foldChange) <- gid
  }
  return(foldChange)
}

# fc_palette <- function(fc) {
# if (all(fc > 0, na.rm=TRUE)) {
# palette <- color_palette(c("blue", "red"))
# } else if (all(fc < 0, na.rm=TRUE)) {
# palette <- color_palette(c("green", "blue"))
# } else {
## palette <- color_palette(c("darkgreen", "#0AFF34", "#B3B3B3", "#FF6347", "red"))
# }
# return(palette)
# }

update_n <- function(x, showCategory) {
  if (!is.numeric(showCategory)) {
    return(showCategory)
  }

  ## geneSets <- geneInCategory(x) ## use core gene for gsea result
  n <- showCategory
  if (nrow(x) < n) {
    n <- nrow(x)
  }

  return(n)
}

extract_geneSets <- function(x, n) {
  n <- update_n(x, n)
  geneSets <- geneInCategory(x) ## use core gene for gsea result
  y <- as.data.frame(x)
  geneSets <- geneSets[y$ID]
  names(geneSets) <- y$Description
  if (is.numeric(n)) {
    return(geneSets[1:n])
  }
  return(geneSets[n]) ## if n is a vector of Description
}

##' Internal plot function for plotting compareClusterResult
##'
##'
##' @title plotting-clusterProfile
##' @param clProf.reshape.df data frame of compareCluster result
##' @param x x variable
##' @param type one of dot and bar
##' @param by one of percentage and count
##' @param title graph title
##' @param font.size graph font size
##' @param colorBy one of pvalue or p.adjust
##' @return ggplot object
##' @importFrom ggplot2 ggplot
##' @importFrom ggplot2 aes
##' @importFrom ggplot2 aes_
##' @importFrom ggplot2 aes_string
##' @importFrom ggplot2 geom_bar
##' @importFrom ggplot2 coord_flip
##' @importFrom ggplot2 geom_point
##' @importFrom ggplot2 %+%
##' @importFrom ggplot2 theme
##' @importFrom ggplot2 xlab
##' @importFrom ggplot2 ylab
##' @importFrom ggplot2 theme_bw
##' @importFrom ggplot2 element_text
##' @importFrom ggplot2 ggtitle
##' @importFrom ggplot2 scale_color_continuous
##' @importFrom ggplot2 guide_colorbar
##' @importFrom DOSE theme_dose
##' @author Guangchuang Yu \url{https://guangchuangyu.github.io}
plotting.clusterProfile <- function(clProf.reshape.df,
                                    x = ~Cluster,
                                    type = "dot",
                                    colorBy = "p.adjust",
                                    by = "geneRatio",
                                    title="",
                                    font.size=12) {
  Description <- Percentage <- Count <- Cluster <- GeneRatio <- p.adjust <- pvalue <- NULL # to satisfy codetools
  if (type == "bar") {
    if (by == "percentage") {
      p <- ggplot(clProf.reshape.df,
                  aes(x=Description, y = Percentage, fill=Cluster))
    } else if (by == "count") {
      p <- ggplot(clProf.reshape.df,
                  aes(x=Description, y = Count, fill=Cluster))
    } else {

    }
    p <- p +
      geom_bar() +
      coord_flip()
  }
  if (type == "dot") {
    if (by == "rowPercentage") {
      p <- ggplot(clProf.reshape.df,
                  aes_(x = x, y = ~Description, size = ~Percentage))
    } else if (by == "count") {
      p <- ggplot(clProf.reshape.df,
                  aes_(x = x, y = ~Description, size = ~Count))
    } else if (by == "geneRatio") {
      p <- ggplot(clProf.reshape.df,
                  aes_(x = x, y = ~Description, size = ~GeneRatio))
    } else {
      ## nothing here
    }
    if (any(colnames(clProf.reshape.df) == colorBy)) {
      p <- p +
        geom_point() +
        aes_string(color=colorBy) +
        scale_color_continuous(low="red", high="blue", guide=guide_colorbar(reverse=TRUE))
      ## scale_color_gradientn(guide=guide_colorbar(reverse=TRUE), colors = sig_palette)
    } else {
      p <- p + geom_point(colour="steelblue")
    }
  }
  p <- p + xlab("") + ylab("") + ggtitle(title) +
    theme_dose(font.size)
  ## theme(axis.text.x = element_text(colour="black", size=font.size, vjust = 1)) +
  ##     theme(axis.text.y = element_text(colour="black",
  ##           size=font.size, hjust = 1)) +
  ##               ggtitle(title)+theme_bw()
  ## p <- p + theme(axis.text.x = element_text(angle=angle.axis.x,
  ##                    hjust=hjust.axis.x,
  ##                    vjust=vjust.axis.x))
  return(p)
}



##' convert compareClusterResult to a data.frame that ready for plot
##'
##'
##' @rdname fortify
##' @title fortify
##' @param includeAll logical
##' @return data.frame
##' @importFrom ggplot2 fortify
##' @importFrom plyr ddply
##' @importFrom plyr mdply
##' @importFrom plyr .
##' @method fortify compareClusterResult
##' @export
##' @author Guangchuang Yu
fortify.compareClusterResult <- function(model, data, showCategory=5, by="geneRatio",
                                         split=NULL, includeAll=TRUE) {
  clProf.df <- as.data.frame(model)
  .split <- split

  ## get top 5 (default) categories of each gene cluster.
  if (is.null(showCategory)) {
    result <- clProf.df
  } else {
    Cluster <- NULL # to satisfy codetools

    topN <- function(res, showCategory) {
      ddply(.data = res,
            .variables = .(Cluster),
            .fun = function(df, N) {
              if (length(df$Count) > N) {
                if (any(colnames(df) == "pvalue")) {
                  idx <- order(df$pvalue, decreasing=FALSE)[1:N]
                } else {
                  ## for groupGO
                  idx <- order(df$Count, decreasing=T)[1:N]
                }
                return(df[idx,])
              } else {
                return(df)
              }
            },
            N=showCategory
      )

    }

    if (!is.null(.split) && .split %in% colnames(clProf.df)) {
      lres <- split(clProf.df, as.character(clProf.df[, .split]))
      lres <- lapply(lres, topN, showCategory = showCategory)
      result <- do.call('rbind', lres)
    } else {
      result <- topN(clProf.df, showCategory)
    }

  }

  ID <- NULL
  if (includeAll == TRUE) {
    result = subset(clProf.df, ID %in% result$ID)
  }

  ## remove zero count
  result$Description <- as.character(result$Description) ## un-factor
  GOlevel <- result[,c("ID", "Description")] ## GO ID and Term
  GOlevel <- unique(GOlevel)

  result <- result[result$Count != 0, ]
  result$Description <- factor(result$Description,
                               levels=rev(GOlevel[,2]))


  if (by=="rowPercentage") {
    Description <- Count <- NULL # to satisfy codetools
    result <- ddply(result,
                    .(Description),
                    transform,
                    Percentage = Count/sum(Count),
                    Total = sum(Count))

    ## label GO Description with gene counts.
    x <- mdply(result[, c("Description", "Total")], paste, sep=" (")
    y <- sapply(x[,3], paste, ")", sep="")
    result$Description <- y

    ## restore the original order of GO Description
    xx <- result[,c(2,3)]
    xx <- unique(xx)
    rownames(xx) <- xx[,1]
    Termlevel <- xx[as.character(GOlevel[,1]),2]

    ##drop the *Total* column
    result <- result[, colnames(result) != "Total"]

    result$Description <- factor(result$Description,
                                 levels=rev(Termlevel))

  } else if (by == "count") {
    ## nothing
  } else if (by == "geneRatio") {
    gsize <- as.numeric(sub("/\\d+$", "", as.character(result$GeneRatio)))
    gcsize <- as.numeric(sub("^\\d+/", "", as.character(result$GeneRatio)))
    result$GeneRatio = gsize/gcsize
    cluster <- paste(as.character(result$Cluster),"\n", "(", gcsize, ")", sep="")
    lv <- unique(cluster)[order(as.numeric(unique(result$Cluster)))]
    result$Cluster <- factor(cluster, levels = lv)
  } else {
    ## nothing
  }
  return(result)
}


##' convert enrichResult object for ggplot2
##'
##'
##' @title fortify
##' @rdname fortify
##' @param model 'enrichResult' or 'compareClusterResult' object
##' @param data not use here
##' @param showCategory Category numbers to show
##' @param by one of Count and GeneRatio
##' @param order logical
##' @param drop logical
##' @param split separate result by 'split' variable
##' @param ... additional parameter
##' @return data.frame
##' @importFrom ggplot2 fortify
##' @method fortify enrichResult
##' @export
fortify.enrichResult <- function(model, data, showCategory=5, by = "Count", order=FALSE,
                                 drop=FALSE, split=NULL, ...) {
  fortify.internal(model, data, showCategory, by, order, drop, split, ...)
}

##' @method fortify gseaResult
##' @export
fortify.gseaResult <- function(model, data, showCategory=5, by = "Count", order=FALSE,
                               drop=FALSE, split=NULL, ...) {
  fortify.internal(model, data, showCategory, by, order, drop, split, ...)
}


fortify.internal <- function(model, data, showCategory=5, by = "Count", order=FALSE,
                             drop=FALSE, split=NULL, ...) {
  res <- as.data.frame(model)
  res <- res[!is.na(res$Description), ]
  if (inherits(model, "gseaResult")) {
    res$Count <- str_count(res$core_enrichment, "/") + 1
    res$.sign <- "activated"
    res$.sign[res$NES < 0] <- "suppressed"
  }
  if (drop) {
    res <- res[res$Count != 0, ]
  }
  if (inherits(model, "gseaResult")) {
    res$GeneRatio <- res$Count / res$setSize
  } else if (inherits(model, "enrichResult")) {
    res$GeneRatio <- parse_ratio(res$GeneRatio)
    if ("BgRatio" %in% colnames(res)) {
      ## groupGO output doesn't have this column
      res$BgRatio <- parse_ratio(res$BgRatio)
    }
  }

  if (order) {
    if (by == "Count") {
      idx <- order(res$Count, decreasing=TRUE)
    } else {
      idx <- order(res$GeneRatio, decreasing=TRUE)
    }
    res <- res[idx,]
  }

  topN <- function(res, showCategory) {
    if ( is.numeric(showCategory) ) {
      if ( showCategory <= nrow(res) ) {
        res <- res[1:showCategory,]
      }
    } else { ## selected categories
      res <- res[res$ID %in% showCategory,]
    }
    return(res)
  }

  if (is.null(split)) {
    res <- topN(res, showCategory)
  } else {
    lres <- split(res, as.character(res[, split]))
    lres <- lapply(lres, topN, showCategory = showCategory)
    res <- do.call('rbind', lres)
  }

  res$Description <- factor(res$Description,
                            levels=rev(unique(res$Description)))

  return(res)
}

str_count <- function(string, pattern="") {
  sapply(string, str_count_item, pattern=pattern)
}

str_count_item <- function(string, pattern = "") {
  length(gregexpr(pattern, string)[[1]])
}

parse_ratio <- function(ratio) {
  gsize <- as.numeric(sub("/\\d+$", "", as.character(ratio)))
  gcsize <- as.numeric(sub("^\\d+/", "", as.character(ratio)))
  return(gsize/gcsize)
}


