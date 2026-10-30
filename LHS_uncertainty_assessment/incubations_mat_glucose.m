function [t,od]=incubations_mat_glucose(par,M)
% calculates OD of one culture based on growth rate, substrate specificity,
% and mortality
itreat=M.itreat;
treatmentmatrix=M.treatmentmatrix;
tspan=[01:1:400]; % time in hours
options=[];
load 'ODmat_corr.mat'
load 'gluc.mat'
%setmom_mat % run script to set parameters
n0=(OD(itreat,1)+0.01)/M.nspec;
N0=treatmentmatrix(itreat,:)*n0; % initial conditions of culture and resouce
y0=[N0 gluc(itreat,1)]; % initial conditions of culture 
[t,od] = ode45(@(t,y)incubation_derivs_mat_glucose(t,y,par),tspan,y0);


end