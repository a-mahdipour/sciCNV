

#' CNV_htmp_gloc function
#'
#' @description Generates heatmap of sciCNV profiles, plotted by Genomic location (gloc)
#'
#' @note Please refer to the reference and supplemental materials described in the README for additional details.
#
#' @param clustering TRUE/FALSE variable specifying whether to cluster cells based on their CNV similarities. Default is "FALSE"
#' @param clustering.type Variable specifying the clustering method to be used in generating the heatmap. Possible options are "pearson",
#'        "euclidean", " spearman", ... "original" (retains the original cell order without unsupervised clustering).
#'         Default is "pearson". Only enabled when clustering = "TRUE"
#' @param sorting TRUE/FALSE variable specifying whether to sort cells based on their tumor CNV score from the largest to smallest tumor scores. Default is FALSE.
#' @param CNVmat copy number variation matrix
#' @param CNVscore is the tumor CNV score matrix for all cells (possibly ranked within clusters). Only used when sorting.clusters = TRUE.
#' @param cluster.lines is a list of values which can be used to separate cell clusters within the population; only used if multiple clusters present
#' @param breakGloc is a set of values each defines a vertical line that separates chromosomes
#' @param No.test number of test cells included in the data; can be used to delineate the populations of test and control cells in the heatmap
#'
#' @author Ali Mahdipour-Shirayeh, Princess Margaret Cancer Centre, University of Toronto
#'
#' @return The output is the heatmap of sciCNV matrix for test and control cells against genomic location
#'
#' @examples
#' 
#' data(breakGloc)
#' file.path <-  system.file("extdata", "CNV_matrix.txt", package="sciCNV")
#' CNVmat <- read.table(file.path, sep="\t", header=TRUE)
#' 
#' CNV_htmp_gloc(CNVmat=CNVmat, breakGloc=breakGloc, sorting = FALSE,  No.test=20,No.normal=20)
#'
#' @import stats
#' @import robustbase
#' @import dichromat
#' @import graphics
#' @import utils
#'
#' @export



CNV_htmp_gloc <- function(CNVmat,
                          clustering = FALSE,        # TRUE or FALSE option
                          clustering.type = NA,     #c("pearson", "kendall", "spearman"),   # "pearson", "kendalln", "spearman", default: "pearson"
                          sorting = TRUE,            # TRUE or FALSE
                          CNVscore = NULL,                # Only exists when sorting = TRUE
                          cluster.lines = NULL,
                          breakGloc,
                          No.test,
                          No.normal
){

  
  genelist <- CNVmat[,1]
  CNVmat <- as.matrix(CNVmat[,-1])
  rownames(CNVmat) <- genelist
  

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
  if ( Reduce("|", is.null(breakGloc)) ){
    breakGloc <- c(0, ncol(CNVmat))
  }

  ##### sorting of cells within clusters, based on tumor CNV scores, from the largest to the smallest (if applicable)
  #nr <- dim(CNVmat)[1]
  seq1 <- seq_len(No.test)
  seq2 <- No.test + seq_len(No.normal)
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
      CNV.mat.clustered <- CNV.mat.tst[hclst.lables, ]

    } else if ( ! missing(clustering.type) ){
      CNV.mat.tst <- CNVmat[seq1, ]
      hclst <- stats::hclust(stats::as.dist(1-stats::cor( t(CNV.mat.tst), method = clustering.type)), method = "ward.D2")
      hclst.lables <- hclst$labels[hclst$order]
      CNV.mat.clustered <- CNV.mat.tst[hclst.lables , ]
    }
    CNVmat.seq2 <- CNVmat[seq2, ]
    CNV.mat1 <- rbind(CNVmat.seq2,as.matrix(CNV.mat.clustered)  )
    rownames(CNV.mat1) <-  c(rownames(CNVmat.seq2), hclst.lables)

  } else if ( (clustering == "FALSE" ) & ( sorting == "FALSE")){
    CNVmat.seq1 <- CNVmat[seq1, ]
    CNVmat.seq2 <- CNVmat[seq2, ]
    CNV.mat1 <- rbind(CNVmat.seq2, CNVmat.seq1 )
    rownames(CNV.mat1) <-  c(rownames(CNVmat.seq2),rownames(CNVmat.seq1))
    colnames(CNV.mat1) <- colnames(CNVmat)

  }

  ROWlist <- rownames(CNV.mat1)
  COLlist <- colnames(CNV.mat1)

  
  

  ## expanding expressions towrds 0 or 1 providing their expression is less than or bigger to 0.5
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
  LL <- 1
  TT <- 0.2
  CNV.mat11[which(CNV.mat11 >  LL)] <-  LL
  CNV.mat11[which( (CNV.mat11 <  TT) & (CNV.mat11 >  -TT) )] <- 0.0
  CNV.mat11[which(CNV.mat11 <  -LL)] <-  -LL

  CNV.mat11 <- as.matrix(CNV.mat11)
  rownames(CNV.mat11) <- matrix(NA, ncol=1, nrow=nrow(CNV.mat1))

  ##################################################
  ## Heatmap of CNV-curves against genomic locations
  ##################################################

  # Uploading the largest list of genes with chr number, start and end
  path <-  system.file("extdata", "10XGenomics_gen_pos_GRCh38-1.2.0.txt", package="sciCNV")
  genLoc <- read.table(path, header=TRUE, sep="\t")
  Specific_genes <- which( as.matrix(genLoc)[, 1]   %in% colnames(CNVmat))
  M_origin <- genLoc[Specific_genes, ]
  M_sample <-  genLoc[Specific_genes, ]

  ## number of segments on the genome
  No_Intrvl <- 1000

  ############ Finding the length of each Chromosome

  minn <- matrix(0, ncol=24, nrow=1)
  maxx <- matrix(0, ncol=24, nrow=1)

  for(i in seq_len(22)){
    chrinfo <- M_origin[which(as.matrix(M_origin[, 2]) == i) , 3]
    if(length(chrinfo) > 0){
      minn[1,i] <- as.numeric(min(chrinfo ) )
      maxx[1,i] <- as.numeric(max(chrinfo ) )
    }
  }
  
  chr23 <- M_origin[which(M_origin[, 2] == "X") , 3]
  if(length(chr23) > 0){
  minn[1,23] <- as.numeric(min(chr23) )
  maxx[1,23] <- as.numeric(max(chr23) )
  }
    
  chr24 <- M_origin[which(as.matrix(M_origin[, 2]) == "Y") , 3]
  if(length(chr24)>0){
  minn[1,24] <- as.numeric(min(chr24) )
  maxx[1,24] <- as.numeric(max(chr24) )
  }

  Total_length <- sum( (maxx[1,1:24]-minn[1,1:24]) + 1 )
  Seg_size <- ceiling( Total_length/No_Intrvl )
  Start_End_Chr <- rbind(minn, maxx)

  ###########
  nonzero <- function(x){x!=0}
  MaxNo_intervals <- ceiling(max((maxx[1,1:24]-minn[1,1:24])+ 1)/Seg_size)
  POINTS <- matrix(0, nrow = 24, ncol = MaxNo_intervals)
  Min_chr <- rep(0,25)

  Min_chr[1] <- 0
  for(i in 1:24){
    StartEnd <- seq(minn[1,i] , maxx[1,i] , Seg_size)
    POINTS[i, 1:length(StartEnd) ]   <-  t(as.matrix(StartEnd  ))
    Min_chr[i+1]  <-  maxx[1,i]
  }

  Length_POINTS <- rep(0,24)
  for(j in 1:24){
    Length_POINTS[j] <- length(POINTS[j,nonzero(POINTS[j, ])]) #nonzero(POINTS[j, ])
  }

  ALL_POINTS <- length(which(POINTS != 0 ))

  ###########

  CNV.mat4 <- matrix(0, nrow = nrow(CNV.mat11), ncol = ALL_POINTS)
  Minnn <- matrix(0, nrow=1,ncol=24)
  Maxxx <- matrix(0, nrow=1,ncol=24)

  for(i in seq_len(nrow(CNV.mat4))){
    for(z in seq_len(22)){
      for( k in seq_len(length(which(POINTS[z, ]>0))-1) ){
        X <- min(which(as.matrix(M_sample[, 2]) == z)):max(which(as.matrix(M_sample[, 2]) == z))
        ExistExpr <- which( (M_sample[ X ,3] >= POINTS[z, k])
                            & (M_sample[ X ,3] <=  POINTS[z, k+1]) )
        if( z == 1){
          if( length(ExistExpr) == 0){
            CNV.mat4[i,k] <- NA
          }else if( length(ExistExpr) == 1){
            CNV.mat4[i,k] <- CNV.mat11[i, ExistExpr]
          }else{
            CNV.mat4[i,k] <- mean( CNV.mat11[i, ExistExpr])
          }

        }else {
          PrePoint_No <-  sum(Length_POINTS[1:(z-1)])
          if( length(ExistExpr) == 0){
            CNV.mat4[i,k + PrePoint_No] <- NA
          }else if( length(ExistExpr) == 1){
            CNV.mat4[i,k + PrePoint_No] <- CNV.mat11[i, ExistExpr + min(X) - 1 ]
          }else{
            CNV.mat4[i,k + PrePoint_No] <- mean( CNV.mat11[i, ExistExpr + min(X) - 1 ])
          }
        }
      }
    }

    ##########
    
    if( length(which(as.matrix(M_sample[, 2]) == "X")) != 0 ){
      X <- min(which(as.matrix(M_sample[, 2]) == "X")):max(which(as.matrix(M_sample[, 2]) == "X"))
      for( k in seq_len(length(which(POINTS[23, ]>0))-1) ){
        
        ExistExpr <- which( (M_sample[ X ,3] >= POINTS[23, k])
                            & (M_sample[ X ,3] <=  POINTS[23, k+1]) )
        PrePoint_No <- sum(Length_POINTS[1:22])
        if( length(ExistExpr) == 0){
          CNV.mat4[i,k + PrePoint_No] <- NA
        }else if( length(ExistExpr) == 1){
          CNV.mat4[i,k + PrePoint_No] <- CNV.mat11[i, ExistExpr + min(X) - 1 ]
        }else{
          CNV.mat4[i,k + PrePoint_No] <- mean( CNV.mat11[i, ExistExpr + min(X) - 1 ])
        }
      }
    }

    ##########
    if( length(which(as.matrix(M_sample[, 2]) == "Y")) != 0 ){
      X <- min(which(as.matrix(M_sample[, 2]) == "Y")):max(which(as.matrix(M_sample[, 2]) == "Y"))
      for( k in seq_len(length(which(POINTS[24, ]>0))-1) ){
        
        ExistExpr <- which( (M_sample[ X ,3] >= POINTS[24, k])
                            & (M_sample[ X ,3] <=  POINTS[24, k+1]) )
        PrePoint_No <- sum(Length_POINTS[1:23])
        if( length(ExistExpr) == 0){
          CNV.mat4[i,k + PrePoint_No] <- NA
        }else if( length(ExistExpr) == 1){
          CNV.mat4[i,k + PrePoint_No] <- CNV.mat11[i, ExistExpr + min(X) - 1  ]
        }else{
          CNV.mat4[i,k + PrePoint_No] <- mean( CNV.mat11[i, ExistExpr + min(X) - 1 ])
        }
      }
    }
  }



  ######### Interpolations

  for(k in seq_len(ALL_POINTS)){

    if( (k == 1) && is.na(CNV.mat4[1,k])){
      BEF <- 1
      AFT <- 1 + min(which( CNV.mat4[ 1 , 2:ncol(CNV.mat4)] != "NA"))
      BEF_Val <- CNV.mat4[, BEF]
      AFT_Val <- CNV.mat4[, AFT]
      #
      CNV.mat4[ , k] <- AFT_Val    #(AFT_Val + BEF_Val)/(AFT-BEF)
      #
    } else if( (k == ncol(CNV.mat4)) && is.na(CNV.mat4[1,k]) ){
      BEF <- max( which(  CNV.mat4[ 1 ,1:(k-1)] != "NA"))
      AFT <- ncol(CNV.mat4) #min(which( MB4[i,(k+1):ncol(MB4)] != "NA" )  )
      BEF_Val <- CNV.mat4[, BEF]
      AFT_Val <- CNV.mat4[, AFT]
      #
      CNV.mat4[ , (k-1)+l] <- BEF_Val    #(AFT_Val + BEF_Val)/(AFT-BEF)
      #
    } else if( is.na(CNV.mat4[1,k] )){
      BEF <-  max( which( CNV.mat4[1,1:(k-1)] != "NA"))
      AFT <-  k + min(which( CNV.mat4[1,(k+1):ncol(CNV.mat4)] != "NA"))
      BEF_Val <- CNV.mat4[, BEF]
      AFT_Val <- CNV.mat4[, AFT]
      BEF_Val[is.na(BEF_Val)] <- 0
      AFT_Val[is.na(AFT_Val)] <- 0
      #
      CNV.mat4[ , k] <- (AFT_Val + BEF_Val)/2 #(AFT_Val + BEF_Val)/(AFT-BEF)
    }

  }


  ###########################################

  CNV.mat51 <- as.matrix(CNV.mat4)
  CNV.mat5 <- matrix(0, ncol=ncol(CNV.mat51), nrow=nrow(CNV.mat51)) # (ceiling(LLL/10)))
  Orig <- 0.5
  LLength <- c(0, cluster.lines[2:(length(cluster.lines)) ])  # c(0,205,733,834,854)

  for(j in seq_len(length(LLength)-1)){
    NeighborNo <- min(20,floor((LLength[j+1]-LLength[j])/2))
    for(i in LLength[j]+seq_len(floor(NeighborNo))){
      CNV.mat5[i, ] <- CNV.mat51[i, ]*Orig + (1-Orig)*robustbase::colMedians( CNV.mat51[ setdiff(seq(LLength[j]+1,(i+NeighborNo),1),i),  ])
    }
    if( (LLength[j]+NeighborNo+1) <= (LLength[j+1]-(NeighborNo))  ){
      for(i in LLength[j]+NeighborNo+seq_len(LLength[j+1]-2*NeighborNo-LLength[j])){
        CNV.mat5[i, ] <- CNV.mat51[i, ]*Orig + (1-Orig)*robustbase::colMedians( CNV.mat51[ setdiff(seq(i-NeighborNo,i+NeighborNo,1),i),  ])
      }
    }
    if( (LLength[j+1]-NeighborNo+1) <= LLength[j+1] ){
      for(i in LLength[j+1]-NeighborNo+seq_len(NeighborNo)){
        CNV.mat5[i, ] <- CNV.mat51[i, ]*Orig + (1-Orig)*robustbase::colMedians( CNV.mat51[setdiff(seq(i-NeighborNo,LLength[j+1], 1),i),  ])
      }
    }

  }

  ########## Defining separating lines
  labels.gloc <- t(as.matrix(c(paste("Chr", 1:22, sep = ""),"ChrX", "ChrY")))
  labels.call.gloc <- rep(NA, ncol(CNV.mat5))
  breakgloc <- unlist(as.list(t(breakGloc)))

  for(i in seq_len(22)){
    labels.call.gloc[breakgloc[i]+10 ] <- as.matrix(labels.gloc[i])
  }
  labels.call.gloc[breakgloc[23] ] <- as.matrix(labels.gloc[23])
  labels.call.gloc[breakgloc[24] ] <- as.matrix(labels.gloc[24])


  rownames(CNV.mat5) <- ROWlist
  final.mat <- CNV.mat5

  ## Sketching the heatmap

  COL_vec <- c( rep("steelblue",1), rep("white",2),rep("firebrick",1) )


  if(clustering == TRUE){
     Separns <- c(1,nrow(final.mat)-No.test,nrow(final.mat))

     graphics::plot.new()
     graphics::par(mar=c(5,5,4,2)+1,mgp=c(3,1,0))
     heatmap.3( final.mat,
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
           labCol = labels.call.gloc,
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
           add.expr =  graphics::abline(v= c(breakgloc))
      )

    } else {

      graphics::plot.new()
      graphics::par(mar=c(5,5,4,2)+1,mgp=c(3,1,0))
      heatmap.3( final.mat ,
             main = paste("Heatmap of sciCNV profiles of test and control cells
                          Thr 0.5 of 1 -",NeighborNo," nearst neighbors", sep="" ),
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
             labCol = labels.call.gloc,
             cexCol = 1,
             srtCol = 90,
             hieght = 50,
             width = 400,
             legend = TRUE,
             margins = c(4,2),
             key.xlab = "Transcription level",
             denscol = "grey",
             density.info = "density",
             rowsep = cluster.lines,
             add.expr = graphics::abline(v= c(breakgloc)))

  }


}












