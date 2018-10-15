
#----------------------------------------------------------------------------------
# tridiagonal matrices
#----------------------------------------------------------------------------------
#
# name:         tridiag.R
# description:  creation of tridiagonal matrices
# inputs:  
#               upper <- upper diagonal
#               lower <- lower diagonal
#               main  <- main diagonal  
# output:       a matrix of size length(main)
#
# creadet by:   Carolina Gamboa
#
#----------------------------------------------------------------------------------



tridiag <- function(upper, lower, main){
  if(length(upper) !=  length(lower) & length(upper) != (length(main)-1)){
    cat("error in especification!! \n")
  }
  else{ 
    out <- matrix(0,length(main),length(main))
    diag(out) <- main
    indx <- seq.int(length(upper))
    out[cbind(indx+1,indx)] <- lower
    out[cbind(indx,indx+1)] <- upper
    return(out)
  }
}