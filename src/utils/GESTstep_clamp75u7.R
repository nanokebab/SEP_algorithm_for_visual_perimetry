GEST.step.clamp75u7 <- function(output,configuration,locations) {
  
  # issues:
  # conversion from actual values to domian postion for clampVector neccessary (and subsequent back conversion)?
  # 
  
  #load functions
  funPath <- "src/utils/"
  source(paste(funPath, "step.R",sep = ""))
  source(paste(funPath, "stop.R",sep = ""))
  
  # DEFS
  maxIter=10
  nNodes <- output$CRF$vf$n.nodes
  reduceResBy <- output$CRF$supportCRF_resReductionFac
  
  # allocation
  vfIntermediateStates <- vector()
  
  # isFinished <- unlist(lapply(output$states, GEST.stop)) | !(seq(1,nNodes) %in% output$CRF$subsets[[output$CRF$currSubset]]) | (sapply(output$states, function(x) x$numPresentations) == output$CRF$currNoPres)
  
  #Loop through until all states are "stop"
  while (!all(isFinished <- unlist(lapply(output$states, GEST.stop)) | 
              !(seq(1,nNodes) %in% output$CRF$subsets[[output$CRF$currSubset]]) | 
              (sapply(output$states, function(x) x$numPresentations) == output$CRF$currNoPres)
  )) {
    
    # isFinished <- unlist(lapply(output$states, GEST.stop)) | !(seq(1,nNodes) %in% output$CRF$subsets[[output$CRF$currSubset]]) | (sapply(output$states, function(x) x$numPresentations) == output$CRF$currNoPres)

    #choose a random,
    i <- which(!isFinished)
    j <- i[runif(1, min = 1, max = length(i))] # unstopped state
    #read out configuration
    oSD <- configuration$oSD
    uSD <- configuration$uSD
    #step it
    r <- GEST.step(output$states[[j]],update_psychoSD = uSD, optimize_psychoSD = oSD,psychoFn = configuration$fn, psychoFp = configuration$fp)
    #update the state
    output$states[[j]] <- r$state
    # if last step of subset -> do following steps:
    if (length(i) == 1) {
      nextSubset <- (output$CRF$currSubset %% 4) + 1
      # 1. send all pdf of subset to CRF
      output$CRF$vf$node.pot[which(output$CRF$vf$node.id %in% output$CRF$subsets[[output$CRF$currSubset]]),] <- t(sapply(output$states[output$CRF$subsets[[output$CRF$currSubset]]], function(x) {x$pdf[seq(1,length(x$domain),reduceResBy)]}))
      # 2. clamp all nodes of current subset to newly acquired values
      newVals <- sapply(output$states[output$CRF$subsets[[output$CRF$currSubset]]], function(x) x$domain[which.min(abs(cumsum(x$pdf) - 0.5))])
      output$CRF$clampVector[output$CRF$subsets[[output$CRF$currSubset]]] <- sapply(round(newVals), function(x) which(output$CRF$supportCRF %in% x))
      # (if now all nodes are clamped, remove 'oldest' subset i.e. next subset)
      if (length(which(output$CRF$clampVector %in% 0)) == 0)
        output$CRF$clampVector[output$CRF$subsets[[nextSubset]]] <- 0
      vf <- clamp.crf(output$CRF$vf, output$CRF$clampVector)
      # 3. get new priors from crf for next subset (run loopy belief propagation)
      newPriors <- infer.lbp(vf,max.iter=maxIter,cutoff=1e-04,verbose=0)$node.bel
      newPriorsForNextSubset <- newPriors[which(vf$node.id %in% output$CRF$subsets[[nextSubset]]),]
      # 4. write these new values into states for GEST algorithmus
      for (k in 1:length(output$CRF$subsets[[nextSubset]])) {
        output$states[[output$CRF$subsets[[nextSubset]][k]]]$pdf <- newPriorsForNextSubset[k,sort(rep(seq(1,ncol(newPriorsForNextSubset),1),reduceResBy))[1:length(output$states[[output$CRF$subsets[[nextSubset]][k]]]$pdf)]]
      }
      # 5. change subgroup
      output$CRF$currSubset <- nextSubset
      # 6. go to next round if every location in the subset was scanned once
      if (nextSubset == 4)
        output$CRF$currNoPres <- output$CRF$currNoPres+1
      # 7. write down intermediate values of thresholds of visual field for plotting
      clampedLocations <- which(output$CRF$clampVector != 0)
      newThresholds <- output$CRF$supportCRF[max.col(newPriors)]
      allThresholds <- vector(,length=nNodes)
      allThresholds[vf$node.id] <- newThresholds
      if (all(allThresholds[clampedLocations] != 0))
        stop("something went wrong while inserting clamped nodes back into grid")
      allThresholds[clampedLocations] <- output$CRF$supportCRF[output$CRF$clampVector[clampedLocations]]
      vfIntermediateStates <- rbind(vfIntermediateStates, allThresholds)
    }
    # update total deviations
    output$totDevs <- c(output$totDevs, sum(abs(sapply(locations, function(x) x[[3]]) - sapply(output$states, function(x) {tail(x$pdfMax,1)}) )) )
  }
  
  return(list(output = output, history = vfIntermediateStates))
}
