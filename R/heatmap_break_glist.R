

#' heatmap_break_glist function
#'
#' @description Finding separating lines among diverse chromosomes to be used in sketching heatmap
#'
#' @param CNV.mat2 the sciCNV matrix of test and controcl cells
#'
#' @author Ali Mahdipour-Shirayeh, Princess Margaret Cancer centre, University of Toronto
#'
#' @return Calculates the seperation spots among chromosomes to sketch the heatmap of sciCNV matrix using CNV_htmp_glist function
#'
#' @examples
#' breakpoints_heatmap <- heatmap_break_glist(CNV.mat2 = normalized_data)
#'
#' @import utils
#'
#' @export



heatmap_break_glist <- function(CNV.mat2){

  gen.loc <- utils::read.table( "../data/10XGenomics_gen_pos_GRCh38-1.2.0.txt", sep = '\t', header=TRUE)
  gene.list <- which( as.matrix(gen.loc)[,1] %in% as.matrix(colnames(CNV.mat2)))
  assoc.chr <-  as.matrix(gen.loc[ gene.list, 2]) # as.matrix(mapply( as.matrix(gen.loc[ gene.list, 2]), FUN = as.numeric) )

  break.glist <- rep(0, 24)

  for(i in 1: 22){
    break.glist[i] <- min( which(assoc.chr == i))
  }
  ##
  if( length(which(rownames(assoc.chr) == "X"))  == 0 ){
    break.glist[23] <- ncol(CNV.mat2)
  } else {
    break.glist[23] <- min( which(rownames(assoc.chr) == "X"))
  }
  ##
  if( length(which(  rownames(assoc.chr) == "Y")) == 0 ){
    break.glist[24] <- ncol(CNV.mat2)
  } else {
    break.glist[24] <- min( which(rownames(assoc.chr) == "Y"))
  }



return(break.glist)

}




