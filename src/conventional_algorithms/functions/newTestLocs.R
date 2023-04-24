newTestLocs <- function(states, finishedLocation) {
  
  quadrants <- list(one = list(first= c(13), second=c(5,6,7,12,14,21,22,23), third = c(1,2,11,19,20)), 
                    two = list(first= c(16), second=c(8,9,10,15,17,24,25,26), third = c(3,4,18,27)),
                    three = list(first= c(39), second=c(30,31,32,38,40,45,46,47), third = c(28,29,37,51,52)),
                    four = list(first= c(42), second=c(33,34,35,41,43,48,49,50), third = c(36,44,53,54)))
  
  quadrant <- which( sapply(quadrants, function(x) !all(!unlist(x) %in% finishedLocation)) )
  # level <- which(quadrants[[quadrant]] %in% finishedLocation)
  level <- which(sapply(quadrants[[quadrant]], function(x) !all(!x %in% finishedLocation)))
  
  if (length(quadrants[[quadrant]][[level]][which(quadrants[[quadrant]][[level]] %in% which(!unlist(lapply(states, dynamic.stop))))]) == 0) {
    if (level != 3) {
      level <- level + 1
    } else {}
  }
  newTestLocs <- quadrants[[quadrant]][[level]][which(quadrants[[quadrant]][[level]] %in% which(!unlist(lapply(states, dynamic.stop))))]
  
  return(newTestLocs)
}