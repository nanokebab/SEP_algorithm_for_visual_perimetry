changeSupport <- function(pdfs,support1,support2) {
  
  if ( !identical( round((support1[[3]]-support1[[2]]),3), round((support2[[3]]-support2[[2]]),3) ) )
    stop("bin size of the two supports must be the same!")
  
  if ( !identical(ncol(pdfs),length(support1)))
    pdfs <- t(pdfs)
  if ( !identical(ncol(pdfs),length(support1)))
    stop("support of input distribution is not equal to support provided")
  
  support1 <- round(support1,1)
  support2 <- round(support2,1)
  
  nOutputCols <- length(support2)
  nOutputRows <- nrow(pdfs)
  output <- matrix(0,nrow=nOutputRows, ncol=nOutputCols)
  
  matchIndexes <- which(!is.na(match(support2,support1)))
  pickIndexes <- which(!is.na(match(support1,support2)))
  
  output[,matchIndexes] <- pdfs[,pickIndexes]
  
  return(output)
}