

#' CNV_htmp_glist function
#'
#' @description Generates a heatmap of single cell CNV profiles plotted by gene list
#'
#' @param CNVmat copy number variation matrix
#' @param Gen.Loc Genomic Location matrix with list of genes, their associated chromosome numbers, starts and ends
#' @param clustering a TRUE/FALSE variable specifying whether to cluster cells based on their CNV similarities. Default is "FALSE"
#' @param clustering.type: variable specifying the method that will be used to cluster sciCNV profiles (cells), if enabled. Possible options are "pearson",
#'        "euclidean", " spearman", ... "original" (retain the same cell order/clusters as original without further clustering).
#'         Default is "pearson". However this is only enabled when clustering = "TRUE"
#' @param sorting: a TRUE/FALSE variable that enables sorting of cells based on their tumor CNV score from the largest to smallest tumor scores. Default is FASLE.
#' @param CNVscore: the CNV score matrix for all cells (optionally divided by pre-determined clusters). Used only when sorting.clusters = TRUE.
#' @param cluster.lines: is an optional list of values demarcating clusters of cells within the population; only used if multiple
#'        cell clusters exist beyond the simple division between test and control populations
#' @param breakGlist: is a set of values each defines a vertical line that separates chromosomes
#' @param No.test: the number of test cells included in the data; can be used to delineate distinct populations of
#'       of test annd control cells in the heatmap
#'
#' @return The output is the heatmap of sciCNV matrix for test and control cells against list of genes
#'
#' @examples
#' CNVmat <- system.file("extdata", "Sample_CNV_matrix.txt", package="sciCNV")
#' breakGlist <- system.file("extdata", "Sample_breakGlist.txt", package ="sciCNV")
#' CNV_htmp_glist(CNVmat, breakGlist=breakGlist, sorting = FALSE,  No.test=20)
#'
#' @import stats
#' @import robustbase
#' @import dichromat
#' @import graphics
#'
#' @export


CNV_htmp_glist <- function(CNVmat,
                           clustering = FALSE,            # TRUE or FALSE
                           clustering.type = c("pearson", "kendall", "spearman"),   # "pearson", "kendalln", " spearman", defualt: "pearson"
                           sorting = FALSE,               # TRUE or FALSE
                           CNVscore = NULL,               # Only used when sorting = TRUE
                           cluster.lines = NULL,
                           breakGlist,     # separation lines for chromosomes
                           No.test
){

  
  
  ## argument validation
  if ( missing(sorting) ){
    sorting <- FALSE
  }

  if  ( is.null(clustering) ){
    clustering <- FALSE
  }
  if ( clustering == "TRUE" ){
    if ( ! clustering.type %in% c("pearson", "kendall", "spearman")  ){
      stop("Please choose a proper method for unsupervised clustering.")
    } else if ( is.null(clustering.type) ){
      clustering.type <- c("pearson")
    }
    if ( sorting == "TRUE" ){
      stop("Sorting can occur when clustering = FALSE.")
    }
  } else if (clustering == "FALSE"){
    if ( (sorting == "TRUE") & Reduce("|", is.null(CNVscore)) ){
      stop("Please insert a list of CNV-scores")
    }
    if ( (sorting == "FALSE") & Reduce("|", ! is.null(CNVscore)) ){
      stop("Please change sorting status to TRUE to sort data based on CNVscore vector.")
    }
  }

  if  ( (is.null(No.test)) & Reduce("|", is.null(cluster.lines)) ){
    stop("Please insert the number of test cells (No.test).")
  }

  if ( Reduce("|", is.null(cluster.lines)) ){
    cluster.lines <- c(0, nrow(CNVmat) - No.test, nrow(CNVmat))
  }
  if ( Reduce("|", is.null(breakGlist)) ){
    breakGlist <- c(0, ncol(CNVmat))
  }

  ##### sorting of cells within each cluster by CNV-score, from the largest to the smallest (if applicable)
  nr <- 40 #nrow(CNVmat)
  seq1 <- seq_len(No.test)
  seq2 <- No.test + seq_len(nr - No.test)
  if ( sorting == TRUE ){
    tst.score <- base::sort(CNVscore[1, seq1] , decreasing=TRUE)     #MMPCs
    ctrl.score <- base::sort(CNVscore[1, seq2] , decreasing=TRUE)  #NBCs
    ranked.col <- as.matrix( c(colnames(t(as.matrix(ctrl.score))), colnames(t(as.matrix(tst.score))  )) )
    CNV.mat1 <- as.matrix( CNVmat[match(ranked.col, rownames(CNVmat)), ])
    rownames(CNV.mat1) <-  ranked.col

  } else if ( clustering == TRUE ){
    
    if ( is.na(clustering.type) ){
      CNV.mat.tst <- as.matrix(CNVmat[seq1, ])
      hclst <- stats::hclust(stats::as.dist(1-stats::cor( t(CNV.mat.tst), method =  "pearson")), method = "ward.D2")
      hclst.lables <- hclst$labels[hclst$order]
      CNV.mat.clustered <- CNV.mat.tst[hclst.lables , ]

    } else if ( ! missing(clustering.type) ){
      CNV.mat.tst <- as.matrix(CNVmat[seq1, ])
      hclst <- stats::hclust(stats::as.dist(1-stats::cor( t(CNV.mat.tst), method = clustering.type)), method = "ward.D2")
      hclst.lables <- hclst$labels[hclst$order]
      CNV.mat.clustered <- CNV.mat.tst[hclst.lables , ]
    }

    CNVmat.seq2 <- as.matrix(CNVmat[seq2, ])
    CNV.mat1 <- rbind(CNVmat.seq2,as.matrix(CNV.mat.clustered))
    rownames(CNV.mat1) <-  c(rownames(CNVmat.seq2), hclst.lables)
    
  } else if ( (clustering == "FALSE" ) & ( sorting == "FALSE")){
    CNVmat.seq1 <- as.matrix(CNVmat[seq1, ])
    CNVmat.seq2 <- as.matrix(CNVmat[seq2, ])
    CNV.mat1 <- rbind(CNVmat.seq2,CNVmat.seq1)
    rownames(CNV.mat1) <-  c(rownames(CNVmat.seq2),rownames(CNVmat.seq1))
    colnames(CNV.mat1) <- colnames(CNVmat)
  }


  ROWlist <- rownames(CNV.mat1)
  COLlist <- colnames(CNV.mat1)

  ## converging single cell copy number values towards integers
  LL1 <- 0.5
  LL2 <- 1.5
  LL3 <- 2.5
  CNV.mat11 <- matrix(0, ncol = ncol(CNV.mat1), nrow = nrow(CNV.mat1))

  for (w in seq_len(ncol(CNV.mat1))){
    for(l in seq_len(nrow(CNV.mat1))){

      if( abs(CNV.mat1[ l, w]) <  LL1 ){
        CNV.mat11[ l, w] <- sign(CNV.mat1[ l, w])*( CNV.mat1[ l, w] )^2
      } else if((abs(CNV.mat1[ l, w]) >=  LL1) & (abs(CNV.mat1[ l, w]) <=  LL2)){
        CNV.mat11[ l, w] <- sign(CNV.mat1[ l, w])*sqrt( abs(CNV.mat1[ l, w]) )
      } else if(abs(CNV.mat1[ l, w]) >  LL2 ){
        CNV.mat11[ l, w] <- 2*sign(CNV.mat1[ l, w])*sqrt( abs(CNV.mat1[ l, w]) )
      }

    }
  }
  #-------
  CNV.mat1 <- as.matrix(CNV.mat11)
  LL <- 1
  TT <- 0.5
  CNV.mat1[which(CNV.mat1 >  LL)] <-  LL
  CNV.mat1[which( (CNV.mat1 <  TT) & (CNV.mat1 >  -TT) )] <- 0.0
  CNV.mat1[which(CNV.mat1 <  -LL)] <-  -LL

  CNV.mat <- as.matrix(CNV.mat1)
  rownames(CNV.mat) <- matrix(NA, ncol=1, nrow=nrow(CNV.mat1))

  ##################################################
  ## Heatmap of sciCNV profiles plotted by gene list
  ##################################################

  ## Using k-nearest neighbors (kNN)-method
  CNV.mat31 <- as.matrix(CNV.mat)
  CNV.mat3 <- matrix(0, ncol=ncol(CNV.mat31), nrow=nrow(CNV.mat31))
  Orig <- 0.5
  LLength <- c(0, cluster.lines[ 2:( length(cluster.lines) ) ])

  for(j in seq_len( length(LLength)-1)){
    NeighborNo <- min(20, floor((LLength[j+1]-LLength[j])/2))

    for(i in LLength[j]+seq_len(LLength[j]+floor(NeighborNo)-LLength[j])){
      CNV.mat3[i, ] <- CNV.mat31[i, ]*Orig + (1-Orig)*robustbase::colMedians(CNV.mat31[ setdiff(seq(LLength[j]+1,(i+NeighborNo),1),i),  ])
    }
    if( (LLength[j]+NeighborNo+1) <= (LLength[j+1]-(NeighborNo) )){
      for(i in LLength[j]+NeighborNo+seq_len(LLength[j+1]-2*NeighborNo -LLength[j])){
        CNV.mat3[i, ] <- CNV.mat31[i, ]*Orig + (1-Orig)*robustbase::colMedians( CNV.mat31[ setdiff(seq(i-NeighborNo,i+NeighborNo,1),i),  ])
      }
    }
    for(i in LLength[j+1]-NeighborNo+seq_len(NeighborNo)){
      CNV.mat3[i, ] <- CNV.mat31[i, ]*Orig + (1-Orig)*robustbase::colMedians(CNV.mat31[setdiff(seq(i-NeighborNo,LLength[j+1], 1),i),  ])
    }

  }

  CNV.mat3 <- as.matrix(CNV.mat3)

  ## Assigning x-axis location to gene-centered data based on genomic location rank

  label.glist <- t(as.matrix(c(paste("Chr", 1:22, sep = ""),"ChrX", "ChrY")))
  label.call.glist <- rep(NA, ncol(CNVmat))
  breakGlist <- unlist(as.list(t(breakGlist)))

  for(i in seq_len(22)){
    if( breakGlist[i] <  ncol(CNVmat) - 10){
      label.call.glist[breakGlist[i]+10 ] <- as.matrix(label.glist[i])
    }
  }

  if( breakGlist[23] <  ncol(CNVmat)  ){
    label.call.glist[breakGlist[23] ] <- as.matrix(label.glist[23])
    label.call.glist[breakGlist[24] ] <- as.matrix(label.glist[24])
  } else {
    label.call.glist[breakGlist[23] ] <- as.matrix(label.glist[23])
  }


  rownames(CNV.mat3) <- ROWlist
  final.mat <- CNV.mat3

  ## Sketching the heatmap
  COL_vec <- c( rep("steelblue",1), rep("white",2),rep("firebrick",1) )

   if(clustering == TRUE){

     Separns <- c(1,nrow(final.mat)-No.test,nrow(final.mat))

     graphics::plot.new()
     graphics::par(mar=c(5,5,4,2)+1,mgp=c(3,1,0))
     heatmap.3( CNV.mat.clustered ,
           main = "Heatmap of sciCNV profiles of test and control cells
           Thr 0.5 of 1",
           xlab="Genomic location of expressed genes",
           ylab= "Cells",
           breaks = seq(-LL, LL, length.out =16),
           col = colorRampPalette(COL_vec, space = "rgb")(15),
           Colv = "NA",
           trace="none",
           treeheight_row = 0.2,
           sepwidth=c(0.2,0.2),
           sepcolor = "black",
           scale= "none",
           labRaw = NA,
           labCol = label.call.glist,
           Rowv = TRUE,
           dendrogram = "row",
           cluster.by.row = TRUE,
           hclust.FUN = hclst,
           cexCol = 1,
           srtCol=90,
           cutree_rows = 2,
           hieght=50, width = 400,
           legend = TRUE,
           margins = c(4,2),
           key.xlab = "Transcription level",
           denscol = "grey",
           density.info = "density",
           rowsep= Separns,
           add.expr =  graphics::abline(v = breakGlist)
      )

     } else {



     graphics::plot.new()
     graphics::par(mar=c(5,5,4,2)+1,mgp=c(3,1,0))
     heatmap.3( final.mat ,
             main = paste("Heatmap of sciCNV profiles of test and control cells
                          Thr 0.5 of 1, ", NeighborNo ," nearst neighbors", sep="" ),
             xlab = "Genomic location of expressed genes",
             ylab= "Cells",
             breaks = seq(-LL, LL, length.out =16),
             col = colorRampPalette(COL_vec, space = "rgb")(15),
             Rowv = FALSE,
             Colv = FALSE,
             trace ="none",
             sepwidth = c(0.2,0.2),
             sepcolor = "black",
             scale = "none",
             labRow = NA,
             labCol = label.call.glist,
             cexCol = 1,
             srtCol = 90,
             hieght = 50,
             width = 400,
             legend = TRUE,
             margins = c(4,2),
             key.xlab = "Transcription level",
             denscol = "grey",
             density.info = "density",
             rowsep = cluster.lines ,
             add.expr = graphics::abline(v = breakGlist))

  }

}










