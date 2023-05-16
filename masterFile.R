# FAST VISUAL FIELD PERIMETRY ACQUISITION
# 
# Abstract:
# In this project a Fast Perimetrical Visual Field Acquisition Algorithm was implemented using the greedy strategy
# and 'communicating' visual field locations represented by a conditional random field (CRF).
# The goal is to achieve higher acquisition speed and/or higher accuracy than using the well established TOP 
# algorithm which is currently used in the OCTUPUS perimeters of HAAG-STREIT.
# The framework is based on the R package 'OPI' which contains a readily implemented Mean-, Median and Mode-QUEST
# algorithm and an interface for Octupus perimeters. The R package 'VisualFields' was used for patient data. 
# The work is based on the following publications: Jedynak et al. (2012), King-Smith et al. (1994) and Watson (1983)
# 
# This file contains all scripts used for simuation, data analysis and creation of plots.
# The following R-packages are required in order to run this code: OPI, VisualFields, ggplot2, gridExtra, sn, CRF 
#
# By Derk Wild, ARTORG Center for Biomedical Engineering Research, University of Bern, 10.3.2016 

setwd("*add-your-working-directory-here*")

# empty memory
rm(list=ls())

# load scripts
path <- "src/utils/"
source(paste(path, "loadScripts.R",sep = ""))
loadScripts()

###############
# RUN SCRIPTS #
###############
# Chapter 1: The Greedy in a Single Step
########################################
# 1.1 Compute mutual information as a function of 'question' i.e. light intensity.
oneStep_mutualInformation()

# 1.2 Display how question to be asked changes with normal prior SD acc. to the different QUEST strategies.
# and the Greedy
oneStep_optimalQuestion_varyNormalPriorSD()

# 1.3 Display how question to be asked changes with log-normal prior SD acc. to the different QUEST strategies.
# and the Greedy
oneStep_optimalQuestion_varyLogNormPriorSD()

# 1.4 Display how question to be asked according to Greedy changes with SD of psychometric function.
oneStep_optimalQuestion_varyPsychoSD()

# 1.5 Display how question to be asked according to Greedy changes with Fn (frac. false neg.) of psychometric function.
oneStep_optimalQuestion_varyPsychoFn()

# 1.6 Display how question to be asked according to Greedy changes with skewness of prior.
oneStep_optimalQuestion_varySkewNormPriorAlpha()


# Chapter 2: The Greedy on a Single (Visual-) location
######################################################
# 2.1 Display pdf's as they get updated with a psychometric function.
oneLocation_PDFs_greedy()

# 2.2 Show single experiment of one location including the given aswers of the simulated patient.
oneLocation_ErrorAndAnswers()

# 2.3 Display a skewness measure (diff. between mean and median) over the coarse of the exmination of one location.
oneLocation_PDFskewness()

# 2.4 Show convergence of Greedy vs mean-, mode-, median-QUEST and the 4-2 Staircase algorithm.
oneLocation_Error_GreedyVsQUESTvs42Staircase()

# 2.5 Show convergence of Greedy vs mean-, mode- and median-QUEST.
oneLocation_Error_GreedyVsQUEST()

# 2.6 Show convergence of Greedies using different combinations of psychometric functions (different SD's) for 
# updating and computing the greedy, respectively.
oneLocation_Error_GreedyVaryPsychoSD()

# 2.7 Show convergence of Greedies using different combinations of psychometric functions (different SD's) for 
# updating and computing the greedy, respectively vs mean-, mode- and median-QUEST algorithms.
oneLocation_Error_GreedyVaryPsychoSD_vs_QUEST()

# 2.8 Plots convergengence speed (Num. of Steps) vs. initial difference between prior mode and true threshold
# for Greedy and QUEST algorithms.
oneLocation_ConvSpeedVsPriorMean_GreedyAndQUESTs()

# 2.9 Plots convergengence speed (Num. of Steps) vs. initial difference between prior mode and true threshold
# for Greedy,QUEST algorithms and 4-2 Staircase algorithm.
oneLocation_ConvSpeedVsPriorMean_GreedyAndQUESTsAnd42Staircase()

# 2.10 Plots convergengence speed (Num. of Steps) vs. initial difference between prior mode and true threshold
# for Greedy algorithms with psychometric functions with different SD's.
oneLocation_ConvSpeedVsPriorMean_GreedyVaryPsychoSD()

# 2.11 Plots convergengence speed (Num. of Steps) vs. SD of psychometric function used for computing the greedy for 
# different initial differences between prior mode and true threshold.
oneLocation_ConvSpeedVsGreedyPsychoSD_varyPriorMean()

# 2.12 Plots heatmap that shows convergence speed (mean error at step 5) and bias (mean error at step 16) for
# different combinations of psychometric funtions (different SD's) for updating and computing the greedy, respectively.
oneLocation_ConvSpeedHeatmap_varyUpdatePsychoSDandEstimationPsychoSD()


# Chapter 3: The Greedy on a complete visual field
##################################################
# 3.1 Shows the end results and final errors of one experiment performed with the following algorithms: Greedy,
# mean-QUEST and 4-2 Staircase.
# IMPORTANT: If the following error appears: 'Error in grid.Call(L_textBounds, as.graphicsAnnot(x$label), x$x, x$y,  : 
# polygon edge not found' run code again
oneEye_oneExperiment_Errors_GreedyVsQUESTvs42Staircase()

# 3.2 Shows the end results and final errors of one experiment performed with the following algorithms: Greedy,
# mean-QUEST and TOP (only staircase from TOP, no spatial information exchanged).
# IMPORTANT: If the following error appears: 'Error in grid.Call(L_textBounds, as.graphicsAnnot(x$label), x$x, x$y,  : 
# polygon edge not found' run code again
oneEye_oneExperiment_Errors_GreedyVsQUESTvsTOP()

# 3.3 Shows the end results and final errors of multiple experiments performed with the following algorithms: 
# Greedy, mean-QUEST and TOP (only staircase from TOP, no spatial information exchanged).
# IMPORTANT: If the following error appears: 'Error in grid.Call(L_textBounds, as.graphicsAnnot(x$label), x$x, x$y,  : 
# polygon edge not found' run code again
oneEye_manyExperiments_Errors_GreedyVsQUESTvs42Staircase()

# 3.4 Displays SD of final threshold estimates vs. the initial deviation of the prior to the true threshold.
oneEye_manyExperiments_FinalThresholdSD_vs_InitialDeviation()


# Chapter 4: The greedy on multiple visual fields
#################################################
# 4.1 Shows an overview over the Rotterdam 'Longitudinal Glaucoma Visual Fields' data set.
rotterdamDataset_overview()

# 4.2 Shows correlations between neighbouring visual field locations from all visual fields in the Rotterdam 
# data set.
rotterdamDataset_neighbourCorrelations()

# 4.3 Pareto optimization between speed and accuracy for the numerous parameters used in the algorithm using
# visual fields from the rotterdam dataset with a MD-distribution according to the whole Rotterdam data set. 
# Here, only one parameter is changed at a time.
mutipleEyes_parameterParetoOptimization_nonUniformSampling_VaryParamSeparately()

# 4.4 Pareto optimization between speed and accuracy for the numerous parameters used in the algorithm using
# visual fields from the rotterdam dataset with a uniform MD-distribution. 
# Here, only one parameter is changed at a time.
mutipleEyes_parameterParetoOptimization_uniformSampling_VaryParamSeparately()

# 4.5 Pareto optimization between speed and accuracy for the numerous parameters used in the algorithm using
# visual fields from the rotterdam dataset with a MD-distribution according to the whole Rotterdam data set.
# Here, a selection of advantageous parameters is changed together.
mutipleEyes_parameterParetoOptimization_nonUniformSampling_VaryParamTogether()

# 4.6 Pareto optimization between speed and accuracy for the numerous parameters used in the algorithm using
# visual fields from the rotterdam dataset with a uniform MD-distribution. 
# Here, a selection of advantageous parameters is changed together.
mutipleEyes_parameterParetoOptimization_uniformSampling_VaryParamTogether()

# Chapter 5: Adding spatial dependencies between visual field locations using a Conditional Random Field (CRF)
##############################################################################################################
# 5.1 Graphical interface for CRF of visual field to modify node potentials, edge potentials and to clamp nodes
# o a certain value to observe how changes propagate through grid.
crfVisualField_Interactive_ClampToTrueThresholds()

# 5.2 Graphical interface for CRF of visual field to clamp certian locations to true threshold and obseerve 
# how probabilities are modifield troughout the grid.
crfVisualField_Interactive_Modify()

# compare strategies:
oneEye_oneExperiment_TestCRF.R

