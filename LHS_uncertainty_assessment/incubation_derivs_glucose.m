function [dtdy]=incubation_derivs_glucose(t,y,par)
% derivatives for incubation experiments
mumax=par(1);
K=par(2);
lambda=par(3);
conv=par(4);


mu=mumax*y(2)./(K+y(2));

dtdy(1)=(mu-lambda)*y(1); % bacteria

dtdy(2)=-conv*mu*y(1); % glucose


dtdy=dtdy';
end