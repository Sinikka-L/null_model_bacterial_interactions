# Null model for bacterial interactions in culture experiment
Consumer-resource model for mono- and co-cultures in incubation experiments of four Roseobacter strains, cultured in a full factorial design (Füssel et al.)

This repository contains the MATLAB code implementing the null hypothesis model described in Füssel et al., "Bacterial interactions shape the molecular composition of dissolved organic matter." Model parameters are fitted independently for each monoculture. The resulting parameter sets are then used to simulate co-cultures comprising all possible combinations of two, three, and four strains. The model accounts exclusively for competition for the sole organic carbon source in the growth medium, glucose. It therefore serves as a null hypothesis, providing a prediction of the expected growth curves in the absence of any interactions beyond resource competition.

The model is a simple consumer-resource model with four parameter for each bacteria strain: maximum uptake rate, substrate affinity (half-saturation constant for glucose), mortality rate and a conversion factor from OD to glucose. The approach is as follows:
1) Inverse modelling of monocultures: Parameters for the four mono-cultures of the strains SH4, SH22, SH24 and SH40 are fitted seperately for each monoculture with a lsqnonlin function in Matlab (model_monoculture_optimisation/addon_opti_glucose.m). Change ispec in line 26 to select one of the cultures 2-5 (1 is control).
2) Forward modelling of co-cultures: The best set of parameters for each mono-culture is then used in a forward mode to simulate the co-cultures (model_monoculture_optimisation/incubationmodel_higherlevels_glucose_1.m).
3) Uncertainty assessment: 10000 simulations with random parameters are performed for each monoculture, and the best 10 are selected based on their RMSE values. These parameter sets are used in forward mode as in 2). Script is in LHS_uncertainty_assessment/incubationmodel_higherlevels_glucose_LHS.m).

contact: Sinikka Lennartz, sinikka.lennartz [at] uni-oldenburg [dot] de
