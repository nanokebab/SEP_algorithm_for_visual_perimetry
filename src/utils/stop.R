GEST.stop <- function (state) 
{
  keepGoing <- ((state$numPresentations < state$maxPresentations) && 
                  (state$currNotSeenLimit < state$minNotSeenLimit) && (state$currSeenLimit < 
                                                                         state$maxSeenLimit) && (((state$stopType == "S") && (ZEST.stdev(state) > 
                                                                                                                                state$stopValue)) || ((state$stopType == "H") && (ZEST.entropy(state) > 
                                                                                                                                                                                    state$stopValue)) || ((state$stopType == "N") && (state$numPresentations < 
                                                                                                                                                                                                                                        state$stopValue))
                                                                                                 || ((state$stopType == "R") && (state$responseCounter[[1]] < state$stopValue || state$responseCounter[[2]] < state$stopValue)) 
                                                                                                 || ((state$stopType == "C") && (state$responseCrossings < state$stopValue)) 
                                                                                                 || ((state$stopType == "D") && ((state$numPresentations < 1) || (abs(tail(state$pdfMean,2)[1]-tail(state$pdfMean,2)[2]) > state$stopValue))) 
                                                                                                 || ((state$stopType == "2D") && ((state$numPresentations < 2) || (abs(state$pdfMean[1]-state$pdfMean[3]) > state$stopValue)))
                                                                                                 ))
  # stopType == "R" is new -> it stops the updating process as soon as 
  # at least 'stopVal' of both answers (TRUE,FALSE) are recieved
  
  # stopType == "C" is new -> it stops after the response has changed 'stopVal' times 
  # i.e. when 'stopVal' crossings happend
  return(!keepGoing)
}
environment(GEST.stop) <- asNamespace('OPI')