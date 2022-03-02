# DisCoNet
Tools for assessing optimal distance correlation threshold for metabolic correlation networks 


R Files in Package:

DisCoNetUtils.R: Contains all functions to be called by other scripts

ComparisonPlots.R: Contains function plotRateofChange(), which generates necessary networks and figures for all analyses except stability analyses and scale free analysis. Returns nothing, generates .png files for all plots in working directory. Arguments with associated defaults are filename, seqmin=0, seqmax=1, seqstep=0.1, method="distance", pthresh=0.05

StabilityAnalysis.R: Contains function stabilityAnalysis(), which generates necessary networks and figures for the “split and compare” stability analysis with resampling. Returns nothing, generates .png files for all plots in working directory. Arguments with associated defaults are filename, seqmin=0, seqmax=1, seqstep=0.1, method="distance", pthresh=0.05, boots=25

TODO add ScaleFreeAnalysis.R
TODO change the way plots are presented/arranged/saved?
