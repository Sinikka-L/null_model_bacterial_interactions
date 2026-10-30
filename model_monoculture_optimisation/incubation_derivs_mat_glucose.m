function [dydt]=incubation_derivs_mat_glucose(t,y,par)
% derivatives for incubation experiments
mumax=par(1,:)';
K=par(2,:)';
lambda=par(3,:)';
conv=par(4,:)';

N=y(1:4);
R=y(end);
K(N==0)=0;
lambda(N==0)=0;

mu=mumax*R./(K+R);
mu(N==0)=0;

dydt(1:4)=(mu-lambda).*N;
dydt(5)=sum(-conv.*mu.*N);

dydt=dydt';
end