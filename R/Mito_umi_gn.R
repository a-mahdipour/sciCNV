

#' Mito_umi_gn function
#'
#' @description Function to sketch mitochondrial transcript content, nUMI and nGene data at cellular level to detect damaged cells within the population
#'
#' @note Please refer to the reference and supplemental materials described in the README for additional details.
#'
#' @param mat raw data matrix
#' @param percent.mito.G percentage of mitochondrial transcripts relative to total number of reads per cell
#' @param nUMI number of unified molecular identification (UMI)
#' @param nGene number of expressed genes
#' @param No.test number of test cells included in the data; can be used to delineate the populations of test and control cells in the heatmap
#' @param drop.mads number of the median absolute deviation (MAD) to drop damaged cells using the percentage of mitochondrial expression level
#'
#' @author Ali Mahdipour-Shirayeh, Princess Margaret Cancer Centre, University of Toronto
#'
#' @return The output represents damaged cells due to the given  for mitochondrial percentage cutoff
#'
#' @examples
#' \dontrun{
#' damaged_cells <- Mito_umi_gn(mat=raw_data, percent.mito.G, nUMI, nGene, No.test=100, drop.mads = 3)
#' }
#'
#' @import stats
#' @import graphics
#' @import scales
#'
#' @export


Mito_umi_gn <- function( mat,
                         percent.mito.G,
                         nUMI,
                         nGene,
                         No.test,
                         drop.mads = 3
){

  if( missing(percent.mito.G)){
    stop("Please insert proper percentage of mitochondrial matrix (percent.mito.G).")
  }
  if( missing(nUMI)){
    stop("Please insert proper number of UMI per cell matrix (nUMI).")
  }
  if( missing(nGene)){
    stop("Please insert proper number of gene per cell matrix (nGene).")
  }
  if( is.null(drop.mads)){
    drop.mads <- 3
  }

  threshold <- max(0.05, mean(percent.mito.G) + (drop.mads)*(stats::mad(percent.mito.G)) )

  graphics::layout(matrix(c(2,1,0,3,5,4,0,6),2,4) ,c(4.5,1,4.5,1),c(1,5), respect = TRUE)
  graphics::par(mar=c(3,3,0,0),mgp=2:0)
  graphics::plot(percent.mito.G ~ t(as.matrix(nUMI[1:No.test])), col=scales::alpha("black",0.2),  pch=16, cex=1.2,
       xlab="nUMI", ylab="Mitochondrial expression (%)"
  )
  with(mat, graphics::abline(h = threshold, lwd = 2, lty = 2, col = scales::alpha("red", 0.8)))
  graphics::legend("topright", bty = "n", lty = 2, col = scales::alpha("red", 0.8), pt.bg = scales::alpha("red", 0.8),
         legend=paste(drop.mads, "MADs above mean :", threshold))

  graphics::par(mar=c(0,3,1,0))
  HST <- graphics::hist(t(as.matrix(nUMI[1:No.test])) ,breaks=100,col="grey",main=NULL,xaxt="n")
  MTPLR <- HST$counts / HST$density
  Dnsty <- stats::density(nUMI)
  Dnsty$y <- Dnsty$y * MTPLR[1]
  graphics::lines(Dnsty,col=scales::alpha("blue",0.7))

  graphics::par(mar=c(3,0,0,1))
  Dnsty <-  stats::density(as.matrix(percent.mito.G))
  HST <- graphics::hist(percent.mito.G ,breaks=100,plot=FALSE)
  BAR <- graphics::barplot(HST$density,horiz=TRUE,space=0,col="grey",main=NULL,xlab="Density")
  SLOPE <- (max(BAR) - min(BAR)) / (max(HST$mids) - min(HST$mids))
  graphics::lines(y=Dnsty$x * SLOPE + (min(BAR) - min(HST$mids) * SLOPE),
        x=Dnsty$y,lwd=2,col=scales::alpha("blue",0.7))


  #----- Selecting damged cells
  damaged_cells <- NULL
  damaged_cells <- as.matrix(which(percent.mito.G > threshold) )

  graphics::par(mar=c(3,3,0,0),mgp=2:0)
  graphics::plot( t(as.matrix(nGene))[1:No.test] ~ t(as.matrix(nUMI))[1:No.test],
        col=scales::alpha("black",0.2),
        pch=16, cex=1.2, xlab="nUMI", ylab="nGene")
  graphics::points( t(as.matrix(nGene))[damaged_cells] ~ t(as.matrix(nUMI))[damaged_cells] ,
          pch=21,cex=1.2,col=scales::alpha("red",0.5),bg=scales::alpha("red",0.3))
  graphics::legend("topleft",bty="n",pch=21,col=scales::alpha("red",0.8),pt.bg=scales::alpha("red",0.8),
         legend="Damaged cells")


  graphics::par(mar=c(0,3,1,0))
  HST <- graphics::hist(t(as.matrix(nUMI))[1:No.test],breaks=100,col="grey",main=NULL,xaxt="n")
  MTPLR <- HST$counts / HST$density
  Dnsty <- stats::density(nUMI)
  Dnsty$y <- Dnsty$y * MTPLR[1]
  graphics::lines(Dnsty,col=scales::alpha("blue",0.7))


  graphics::par(mar=c(3,0,0,1))
  Dnsty <-  stats::density(as.matrix(nGene[1:No.test]))
  HST <- graphics::hist( nGene ,breaks=100,plot=FALSE)
  BAR <- graphics::barplot(HST$density,horiz=TRUE,space=0,col="grey",main=NULL,xlab="Density")
  SLOPE <- (max(BAR) - min(BAR)) / (max(HST$mids) - min(HST$mids))
  graphics::lines(y=Dnsty$x * SLOPE + (min(BAR) - min(HST$mids) * SLOPE),
        x=Dnsty$y,lwd=2,col=scales::alpha("blue",0.7))


  graphics::title("Detecting damaged cells based on mitochonrial expression level",
        outer = TRUE, line = -2,
        cex.main = 2, col.main ="brown")


  return(damaged_cells)


}





