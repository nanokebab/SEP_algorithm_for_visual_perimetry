fourTwo.unfinishedFinal <- function (state) 
{
  if (state$finished != "Not") 
    return(state$stairResult)
  else {
    warning("fourTwo.step: asking for final result of unfinished staircase")
    # return(state$lastSeen)
    return(state$currentLevel)
  }
}

environment(fourTwo.unfinishedFinal) <- asNamespace('OPI')