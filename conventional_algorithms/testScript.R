#test script for Serife

#DEFS
simType <- "G" # "G" -> patinet with glaucoma, "N" -> normal patient, "C" -> combined

# load packages
require(OPI)
require(visualFields)

#load functions
funPath <- "functions/"
source(paste(funPath, "dynamicWrapper.R",sep = ""))
source(paste(funPath, "TopWrapper.R",sep = ""))

# initialize OPI
chooseOpi("SimHenson")
if (!is.null(opiInitialize(type = simType, cap = 6)))
  stop("opiInitialize failed")

# initialize locations (List of (x, y, true threshold, normative value) triples)
location <- list(c(-9,21,20,25.73017))
locations <- rep(location,54)

result_dynamicStrategy <- dynamicWrapper(locations)

result_topStrategy <- TOPwrapper(locations)
