GESTwrapper.clamp75u7 <- function(configuration, locations, support) {
  
  # load functions
  funPath <- "src/utils/"
  source(paste(funPath, "GESTstart_clamp75u7.R",sep = ""))
  source(paste(funPath, "GESTstep_clamp75u7.R",sep = ""))
  
  output <- GEST.start.clamp75u7(configuration, locations, support)
  output <- GEST.step.clamp75u7(output, configuration, locations)
  
  #write down final stats
  thresholdDeviations <- sapply(locations, function(x) x[3]) - sapply(output$output$states, function(x) x$domain[which.min(abs(cumsum(x$pdf) - 0.5))])
  nrOfSteps <- sum(sapply(output$output$states,function(x) x$numPresentations))
  totDev <- output$output$totDevs
  
  return(list(thresholdDeviations, nrOfSteps, totDev))
}
