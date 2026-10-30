%%% optimization of parameters for mono-cultures
clc
clear all
close all

global ispec
% ispec is the number of treatment (fitting averages), i.e. 1=control

treatmentstr={'control','SH4-1','SH22-1','SH24-1b','SH40',...
    'SH4-1 + SH22-1','SH4-1 + SH24-1b','SH4-1 + SH40','SH22-1 + SH24-1b',...
    'SH22-1 + SH40', 'SH24-1b + SH40', 'SH4-1 + SH22-1 + SH24-1b',...
    'SH4-1 + SH22-1 + SH40','SH4-1 + SH24-1b + SH40','SH22-1 + SH24-1b + SH40',...
    'SH4-1 + SH22-1 + SH24-1b + SH40'}; % string describing treatments

% set options for nonlinear optimization routine
options=optimset('lsqnonlin');
options.TolX=1e-20;
options.TolFun=1e-20;
options.Display='iter';
options.Jacobian='off';
options.Algorithm='levenberg-marquardt';
options.MaxFunctionEvaluations=1000;
lb=[0,0,0,0];
ub=[];
%warning('off');
ispec=2; % 2,3,4,5 !!!! CHANGE HERE FOR FITTING DIFFERENT MONOCULTURES: 2, 3, 4, 5 (1=control)


par=[1,50,0.02,1]; % starting values

[phat,resnorm,residual,exitflag,output,lambda,Jacob]=lsqnonlin(@lsq_addon_glucose,par, lb, ub, options);
ci=nlparci(phat, residual, 'jacobian', Jacob);  

%% RESIDUAL CHECK
load 'ODmat_corr.mat'% measured OD
load 'gluc.mat' % measured glucose 
% time points for OD and glucose measurements:
tOD=[0.00	15.33	38.50	63.25	86.50	110.50	135.50	159.43	183.33...
    207.33	231.00	255.50	302.67	374.88];
tgluc=[0.00	63.25	135.50	255.50];

% run model with best parameter set derived from optimization above
[mt,mod]=incubations_glucose(phat);
% determine residuals
res=lsq_addon(phat);
% determine root mean square error
RMSE=sqrt(mean(res.^2));

f=figure;
set(f,'units','centimeters','position',[5 5 20 12])
subplot(2,1,1) % 
hold on
errorbar(tOD, OD(ispec,:),ODstd(ispec,:),'k','LineStyle','none');
p1=plot(mt,mod(:,1),'-','LineWidth',2,'Color',[53,151,143]/255);
hold on
p2=plot(tOD,OD(ispec,:),'o','MarkerFaceColor', [253,174,97]/255, 'MarkerEdgeColor','k');
yl=ylabel('OD');
title(treatmentstr(ispec))
legend([p1 p2],'modelled','measured','location','eastoutside')
xlim([0 405])
ylim([0 0.7])
t1=text(300,0.58,['\mu_{max}=' sprintf('%.2f',phat(1))]);
t2=text(300,0.5,['K=' sprintf('%.2f',phat(2))]);
t3=text(300,0.42,['\lambda=' sprintf('%.2f',phat(3))]);
t4=text(300,0.35,['RMSE=' sprintf('%.4f',RMSE)]);

subplot(2,1,2)
hold on
errorbar(tgluc, gluc(ispec,:),glucstd(ispec,:),'k','LineStyle','none');
p1=plot(mt,mod(:,2),'-','LineWidth',2,'Color',[53,151,143]/255);
p2=plot(tgluc, gluc(ispec,:),'o','MarkerFaceColor', [253,174,97]/255, 'MarkerEdgeColor','k');
xl=xlabel('incubation time [h]');
yl=ylabel('Glucose [mM]');
title(treatmentstr(ispec))
legend([p1 p2],'modelled','measured','location','eastoutside')
ylim([0 3.5])