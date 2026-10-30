function [t,od]=incubations_glucose(par)
% calculates OD of one culture based on growth rate, substrate specificity,
% and mortality
global ispec
tspan=[0 400]; % time in hours
options=[];
load 'ODmat_corr.mat'

load 'gluc.mat'
y0=[OD(ispec,1)+0.01 gluc(ispec,1)]; % initial conditions of culture and resouce

[t,od] = ode45(@(t,y)incubation_derivs_glucose(t,y,par),tspan,y0);


end