

#' heatmap_break_gloc function
#'
#' @description Finding separating lines among diverse chromosomes to be used in sketching heatmap
#'
#' @param CNV.mat2 copy number variation matrix
#' 
#' @author Ali Mahdipour-Shirayeh, Princess Margaret Cancer centre, University of Toronto
#'
#' @return Calculates the seperation spots among chromosomes to sketch the heatmap of sciCNV matrix using CNV_htmp_gloc function
#'
#' @examples
#' breakpoints_heatmap <- heatmap_break_gloc()
#'
#' @export



heatmap_break_gloc <- function(CNV.mat2 ){
  
  
  Specific_genes <- which( as.matrix(genLoc1)[, 1]   %in% colnames(CNV.mat2))
  M_origin <- genLoc1[Specific_genes, ]
  
  ## number of segments on the genome
  No_Intrvl <- 1000
  
  ############
  minn <- matrix(0, ncol=24, nrow=1)
  maxx <- matrix(0, ncol=24, nrow=1)
  
  for(i in 1: 22){
    minn[1,i] <- as.numeric(min(M_origin[which(as.matrix(M_origin[, 2]) == i) , 3]) )
    maxx[1,i] <- as.numeric(max(M_origin[which(as.matrix(M_origin[, 2]) == i) , 3]) )
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
  
  ########
  
  MaxNo_intervals <- ceiling(max((maxx[1,1:24]-minn[1,1:24])+ 1)/Seg_size)
  Intrvls <- matrix(0, nrow = 24, ncol = MaxNo_intervals)
  Min_chr <- rep(0,25)
  
  Min_chr[1] <- 0
  for(i in 1:24){
    StartEnd <- seq(minn[1,i] , maxx[1,i] , Seg_size)
    Intrvls[i, 1:length(StartEnd) ]  <-  t(as.matrix( StartEnd )) + Min_chr[i]
    Min_chr[i+1]  <-  maxx[1,i]
  }
  
  
  
  break.gloc <- rep(0,24)
  break.gloc[1] <- 1
  
  for(i in 2: 24){
    break.gloc[i] <- length( Intrvls[1:(i-1),][ Intrvls[1:(i-1),] != 0])
  }
  
  return(break.gloc)
  
}




