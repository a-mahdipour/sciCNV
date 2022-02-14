

#' Sketch_AveCNV function
#'
#' @description Plots the average pre-sciCNV profile of test cells by genomic location
#'
#' @param Ave.mat the average of sciCNV-curves of test cells
#'
#' @author Ali Mahdipour-Shirayeh, Princess Margaret Cancer Centre, University of Toronto
#'
#' @return Sketches the average sciCNV curve of all test cells; showing all CN gains/losses across entire genome
#'
#' @examples
#' \dontrun{
#' sketch_mean_sciCNV <- Sketch_AveCNV(Ave.mat=mean_sciCNV_mat)
#' }
#'
#' @import utils
#' @import graphics
#'
#' @export


Sketch_AveCNV <- function(Ave.mat, Assoc.Chr){


  graphics::plot.new()
  graphics::par(mar=c(5,5,4,2)+1,mgp=c(3,1,0))
  graphics::plot( Ave.mat,
        col="brown1",
        type="l",
        lty=1,
        lwd=3,
        xlim=c(-nrow(as.matrix(Ave.mat))*.02,nrow(as.matrix(Ave.mat))*1.02),
        ylim=c(-2,2),
        cex=1,
        xlab="Genomic location",
        ylab="Preliminary CNV estimate")
  graphics::abline(h=0, col="black")

  M <- as.matrix(Assoc.Chr) # Chromosome numbers
  Break <- matrix(0, ncol = 24, nrow = 1)
  for(i in 1: 22){
    Break[1,i] <- apply(as.matrix(Assoc.Chr) == i, 2, which.max)
  }
  if( length(which( as.matrix(Assoc.Chr) == "X") > 0) ){
    Break[1,23] <- apply(as.matrix(Assoc.Chr) == "X", 2, which.max)
    if( length(which( as.matrix(Assoc.Chr) == "Y") > 0) ){
      Break[1,24] <- apply(as.matrix(Assoc.Chr) == "Y", 2, which.max)
    } else{
      Break[1,24] <- nrow(as.matrix(Ave.mat))
    }
  } else{
    Break[1,23] <- nrow(as.matrix(Ave.mat))-1
    Break[1,24] <- nrow(as.matrix(Ave.mat))
  }


  Break_lines <- as.matrix(c(Break, nrow(Assoc.Chr)))
  graphics::abline(v=Break_lines, col="gray65", lty=2)
  graphics::par(new=TRUE);
  graphics::points( Ave.mat,
          col = "brown1",
          type = "l",
          lty = 1,
          lwd = 3
          )
  graphics::points(Break , matrix(c(1.8, 1.7 ), ncol=24, nrow=1), pch=16, col="royalblue1", cex=4)
  graphics::text(Break , matrix(c(1.8, 1.7 ), ncol = 24, nrow = 1), c(seq(1, 22, 1),"X","Y"),
        col = "white", cex = 1.2)


  graphics::title("Average pre-sciCNV profile of test cells relative to control cells")

}



