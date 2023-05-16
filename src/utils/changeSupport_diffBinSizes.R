changeSupport <- function(pdfs,support1,support2) {
  
  if ( !identical(ncol(pdfs),length(support1)))
    pdfs <- t(pdfs)
  if ( !identical(ncol(pdfs),length(support1)))
    stop("support of input distribution is not equal to support provided")
  
  # NEW
  ###############
  supportedValues <- which(support1 %in% support2)
  # pdfs <- t(as.matrix(pdfs[,supportedValues]))
  pdfs<-t(apply(pdfs,1,function(x) x[supportedValues]))
  support1 <- support1[supportedValues]
  #############
  
  support1 <- round(support1,1)
  support2 <- round(support2,1)
  
  nOutputCols <- length(support2)
  nOutputRows <- nrow(pdfs)
  output <- matrix(0,nrow=nOutputRows, ncol=nOutputCols)
  
  if ( !identical( round((support1[[3]]-support1[[2]]),3), round((support2[[3]]-support2[[2]]),3) ) ) {
    print("bin size of the two supports must be the same!")
    supportedValues2 <- which(support2 %in% support1)
    emptyValues <- supportedValues2[[2]] - supportedValues2[[1]]
    matchIndexes <- seq(1:length(support2))
    pickIndexes <- rep(seq(1:length(support1)),each=emptyValues)[1:length(support2)]
    
  } else {
    matchIndexes <- which(!is.na(match(support2,support1)))
    pickIndexes <- which(!is.na(match(support1,support2)))
  }
  
  output[,matchIndexes] <- pdfs[,pickIndexes]
  
  return(output)
}