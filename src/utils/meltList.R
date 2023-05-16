meltList <- function (list, id) {
  
  meltedListDf <- data.frame()
  id <- which(names(meanDeviations) == id)
  for (i in 1:length(list)) {
    if (i != id) {
      entries <- length(list[[i]])
      currentRow <- nrow(meltedListDf)
      meltedListDf[(currentRow + 1):(currentRow + entries),1] <- rep(names(list)[i],entries)
      meltedListDf[(currentRow + 1):(currentRow + entries),2] <- list[i]
      meltedListDf[(currentRow + 1):(currentRow + entries),3] <- list[[id]][1:entries]
    }
  }
  colnames(meltedListDf)[1] <- "Method"
  colnames(meltedListDf)[2] <- "value"
  colnames(meltedListDf)[3] <- "steps"
  
  return(meltedListDf)
}