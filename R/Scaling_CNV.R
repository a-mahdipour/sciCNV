

#' Scaling_sciCNV function
#'
#' @description To scale preliminary sciCNV-curves derived from scRNA-seq data
#'
#' @note Please refer to the reference and supplemental materials described in the README for additional details.
#'
#' @param V7Alt preliminary CNV matrix egnerated by sciCNV function
#' @param Vn.TestCells number of test cells
#' @param Vscaling.factor scaling factor to re-scale the preliminary CNV estimate
#'
#' @author Ali Mahdipour-Shirayeh, Princess Margaret Cancer Centre, University of Toronto
#'
#' @return Scales sciCNV results to be matched with benchmark (e.g. WES result of the same test sample)
#'
#' @examples
#' scaled_sciCNV <- Scaling_CNV(V7Alt=sciCNV_mat, n.TestCells=100, scaling.factor=0.5)
#'
#' @import utils
#' @import graphics
#'
#' @export


Scaling_CNV <-  function(V7Alt,
                         n.TestCells,
                         scaling.factor
){

  # argument validation
  if ( Reduce("|",is.na(V7Alt)) ){
    stop( "The CNV-matrix has not been properly prepared.")
  }

  if ( is.na( n.TestCells) ){
    stop( "Number of tumor cells needs to be inserted.")
  }

  if ( is.na(scaling.factor) ){
    scaling.factor <- 1.0
  }

  Ave_control <- matrix(0, ncol=1, nrow= nrow(V7Alt ))
  Ave_test <- matrix(0, ncol=1, nrow= nrow(V7Alt ))
  Z2 <- seq( n.TestCells + 1, ncol(V7Alt), 1)   #Normal cells
  Z3 <- seq(1, n.TestCells , 1)     #tumor/test cells

 ##
  for(i in 1:nrow(V7Alt)){
    if( !sum(V7Alt[ i, Z2]) == 0 ){
      Ave_control[i, 1] <- mean(as.matrix(V7Alt[ i, Z2][which(!V7Alt[ i, Z2] == 0 ) ]))
    } else {
      Ave_control[i, 1] <- 0
    }
    Ave_test[i, 1] <- mean(as.matrix(V7Alt[ i, Z3]))
  }

  ## Scaling data if needed
    V7Alt <- V7Alt*(1/scaling.factor)
    Ave_control <- Ave_control*(1/scaling.factor)
    Ave_test <- Ave_test*(1/scaling.factor)

    Ave_control[is.na(Ave_control)] <- 0
    Ave_test[is.na(Ave_test)] <- 0

    ##

    Gen.Loc <- utils::read.table( "../data/10XGenomics_gen_pos_GRCh38-1.2.0.txt", sep='\t', header=TRUE)
    Specific_genes <- which( as.matrix(Gen.Loc)[, 1]   %in% rownames(V7Alt))
    M_sample <-  Gen.Loc[Specific_genes, ]

    Chr_end <- matrix(0, ncol=24, nrow=1)
    Chr_begin <- matrix(0, ncol=24, nrow=1)

    for(l in 1:22){
      Chr_begin[1,l] <- min(which(M_sample [,2]==l))
      Chr_end[1,l] <- max(which(M_sample [,2]==l))
    }
    Chr_begin[1,23] <- min(which(M_sample [,2]=="X"))
    Chr_end[1,23] <- max(which(M_sample [,2]=="X"))
    ##
    Chr_begin[1,24] <- min(which(M_sample [,2]=="Y"))
    Chr_end[1,24] <- max(which(M_sample [,2]=="Y"))

  ## To scale data we used some genes that are correlated to specific number of CNV
  graphics::plot.new()
  graphics::par(mar=c(5,5,4,2)+1,mgp=c(3,1,0))
  graphics::plot(Ave_control[, 1],
       col= "blue",
       pch=16,
       cex=0.5,
       ylim=c(-2,2),
       xlab="Genomic location",
       ylab="Preliminary CNV estimate"
  )
  graphics::points(Ave_test[, 1], col= "red", pch=16, cex=0.5 )
  graphics::abline(h=0, col="black")
  graphics::abline(h=c(-scaling.factor, scaling.factor), lty=1, col="orange")

  graphics::legend("bottomleft",
         legend=c( "Ave control","Ave test"),
         inset=.02, col=c("blue","red"),
         fill=c("blue", "red"),
         horiz=TRUE,
         cex=1.2)

  graphics::abline(v=c(Chr_begin,Chr_end), col="gray65", type="l", lty=2)
  graphics::points(Chr_begin+100, matrix(c(  1.8, 1.7 ),ncol=24, nrow=1), pch=16, col="royalblue1", cex=4)
  graphics::text(Chr_begin+100, matrix(c(  1.8, 1.7 ),ncol=24, nrow=1), c(seq(1,22,1),"X","Y"),col="white", cex=1.2)

  graphics::title(paste("Average pre-sciCNV profile of test/control cells - Threshold=",
              round(scaling.factor,
              digits=2)),
              col.main="brown",
              cex.main=2)

  Final_Mat <- cbind(V7Alt, Ave_test)
  rownames(Final_Mat) <- rownames(V7Alt)
  colnames(Final_Mat) <- c(colnames(V7Alt),c("AveTest"))

  ## this function returnes the scaled CNV-matrix
  return( Final_Mat )



}




