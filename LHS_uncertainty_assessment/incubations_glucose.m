function [t,od]=incubations_glucose(par,M)
% calculates OD of one culture based on growth rate, substrate specificity,
% mortality, and conversion factor OD -> glucose
ispec=M.itreat;
tspan=[0 400]; % time in hours
options=[];
load 'ODmat_corr.mat'
load 'gluc.mat'
y0=[OD(ispec,1)+0.005 gluc(ispec,1)]; % initial conditions of culture and resouce, set to nonzero to get model started
[t,od] = ode45(@(t,y)incubation_derivs_glucose(t,y,par),tspan,y0);


end