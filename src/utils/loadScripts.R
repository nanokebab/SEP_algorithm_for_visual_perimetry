loadScripts <- function() {
  
  # load scripts
  funPath <-"src/scripts/"
  
  source(paste(funPath, "oneStep_mutualInformation.R",sep = ""))
  source(paste(funPath, "oneStep_optimalQuestion_varyNormalPriorSD.R",sep = ""))
  source(paste(funPath, "oneStep_optimalQuestion_varyLogNormPriorSD.R",sep = ""))
  source(paste(funPath, "oneStep_optimalQuestion_varyPsychoSD.R",sep = ""))
  source(paste(funPath, "oneStep_optimalQuestion_varyPsychoFn.R",sep = ""))
  source(paste(funPath, "oneStep_optimalQuestion_varySkewNormPriorAlpha.R",sep = ""))
  
  source(paste(funPath, "oneLocation_PDFs_greedy.R",sep = ""))
  source(paste(funPath, "oneLocation_Error&Answers.R",sep = ""))
  source(paste(funPath, "oneLocation_PDFskewness.R",sep = ""))
  source(paste(funPath, "oneLocation_Error_GreedyVsQUESTvs4-2Staircase.R",sep = ""))
  source(paste(funPath, "oneLocation_Error_GreedyVsQUEST.R",sep = ""))
  source(paste(funPath, "oneLocation_Error_GreedyVaryPsychoSD.R",sep = ""))
  source(paste(funPath, "oneLocation_Error_GreedyVaryPsychoSD_vs_QUEST.R",sep = ""))
  source(paste(funPath, "oneLocation_ConvSpeedVsPriorMean_Greedy&QUESTs.R",sep = ""))
  source(paste(funPath, "oneLocation_ConvSpeedVsPriorMean_Greedy&QUESTs&4-2Staircase.R",sep = ""))
  source(paste(funPath, "oneLocation_ConvSpeedVsPriorMean_GreedyVaryPsychoSD.R",sep = ""))
  source(paste(funPath, "oneLocation_ConvSpeedVsGreedyPsychoSD_varyPriorMean.R",sep = ""))
  source(paste(funPath, "oneLocation_ConvSpeedHeatmap_varyUpdatePsychoSD&EstimationPsychoSD.R",sep = ""))
  source(paste(funPath, "oneEye_oneExperiment_Errors_GreedyVsQUESTvs4-2Staircase.R",sep = ""))
  source(paste(funPath, "oneEye_oneExperiment_Errors_GreedyVsQUESTvsTOP.R",sep = ""))
  source(paste(funPath, "oneEye_manyExperiments_Errors_GreedyVsQUESTvs4-2Staircase.R",sep = ""))
  source(paste(funPath, "oneEye_manyExperiments_FinalThresholdSD_vs_InitialDeviation.R",sep = ""))

  source(paste(funPath, "RotterdamDataset_overview.R",sep = ""))
  source(paste(funPath, "RotterdamDataset_neighbourCorrelations.R",sep = ""))  
  source(paste(funPath, "mutipleEyes_parameterParetoOptimization_nonUniformSampling_VaryParamSeparately.R",sep = "")) 
  source(paste(funPath, "mutipleEyes_parameterParetoOptimization_uniformSampling_VaryParamSeparately.R",sep = "")) 
  source(paste(funPath, "mutipleEyes_parameterParetoOptimization_nonUniformSampling_VaryParamTogether.R",sep = "")) 
  source(paste(funPath, "mutipleEyes_parameterParetoOptimization_uniformSampling_VaryParamTogether.R",sep = "")) 
  
  source(paste(funPath, "crfVisualField_Interactive_ClampToTrueThresholds.R",sep = ""))
  source(paste(funPath, "crfVisualField_Interactive_Modify.R",sep = ""))
}
  
  
