

#' CNV_score function
#'
#' @description Calculates a score for each cell reflecting the similarity of its sciCNV
#'              profile to the mean CNV profile of all suspected tumor cells
#'
#' @note Please see the reference and supplmental materials described in the README file for additional information.
#
#' @param M_nF is a matrix of sciCNV profiles, with an extra column that reflects the average
#'       sciCNV profile of all suspected tumor cells (test cells)
#'
#' @author Ali Mahdipour-Shirayeh, Princess Margaret Cancer centre, University of Toronto
#'
#' @return Calculates the tumor likeness score for each single-cell which is the similarity of sciCNV curve of each cell to average sciCNV of tet cells
#'
#' @examples
#' tumor_likeness_score <- CNV_score(M_nf=sciCNV_mat)
#'
#' @export


CNV_score <-  function( M_nf ){

  # argument validation
  if ( is.matrix(M_nf) == FALSE ){
    stop( "Please insert a proper matrix of CNV-curves")
  }

  #-------
  LL1 <- 0.5
  LL2 <- 1.5
  LL3 <- 2.5
  M_nf.scaled <- matrix(0, ncol=ncol(M_nf), nrow=nrow(M_nf))

  for (w in 1:(ncol(M_nf)-1)  ){
    for(l in 1: nrow(M_nf)){
      if( abs( M_nf[ l, w] ) <  LL1 ){
        M_nf.scaled[ l, w] <- sign(M_nf[ l, w] )*( M_nf[ l, w] )^2
      } else if( (abs(M_nf[ l, w]) >=  LL1) & ( abs(M_nf[ l, w] ) <=  LL2) ){
        M_nf.scaled[ l, w] <- sign(M_nf[ l, w] )*sqrt( abs(M_nf[ l, w] ) )
      } else if( M_nf[ l, w]  >  LL2 ){
        M_nf.scaled[ l, w] <- 2*sign( M_nf[ l, w] )*sqrt( M_nf[ l, w]  )
      }

    }
  }
  M_nf.scaled[ , ncol(M_nf)] <- t(as.matrix(M_nf[ , ncol(M_nf)]))

  ##################
    
  Score <- matrix(0, nrow = nrow(M_nf.scaled), ncol = ncol(M_nf.scaled) )
  TotScore <- matrix(0, nrow = 1, ncol = (ncol(M_nf.scaled)-1) )
  Wr <- ncol(M_nf.scaled)-1

  for(w in 1:(ncol(M_nf.scaled)-1)){
    for(l in 1:nrow(M_nf.scaled)){
      if(M_nf.scaled[l,ncol(M_nf.scaled)] != 0){
        Score[l,w] <- sign(M_nf.scaled[l,ncol(M_nf.scaled)])*sign(M_nf.scaled[l,w])*sqrt(abs(M_nf.scaled[l,w]))
      } else {
        Score[l,w] <- -sqrt( abs(M_nf.scaled[l,w]) )
      }
    }
  }

  #--------------------------------
  Score[1:nrow(M_nf.scaled), ncol(M_nf.scaled)] <- as.matrix(sign(M_nf.scaled[1:nrow(M_nf.scaled),
                                                   ncol(M_nf.scaled)])*sqrt( abs(M_nf.scaled[,ncol(M_nf.scaled)])))
  TotScore[1, ] <- 100*colSums(Score[, 1:(ncol(M_nf.scaled)-1)])/sum(Score[,ncol(M_nf.scaled)])
  colnames(TotScore) <- colnames( M_nf[, -ncol(M_nf)])

  return( TotScore )

}









