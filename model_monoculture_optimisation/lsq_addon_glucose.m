% lsq_OCS
function res=lsq_addon_glucose(par)
global ispec
[mt,mod]=incubations_glucose(par); %model time and model od
%ispec=5;
load 'ODmat_corr.mat'% measured OD
load 'gluc.mat' % measured glucose
tOD=[0.00	15.33	38.50	63.25	86.50	110.50	135.50	159.43	183.33...
    207.33	231.00	255.50	302.67	374.88];
tgluc=[0.00	63.25	135.50	255.50];
index=indexclosestto(mt, tOD);
indexg=indexclosestto(mt, tgluc);

resod=abs(mod(index,1)-OD(ispec,:)');
resgluc=abs(mod(indexg,2)-gluc(ispec,:)');

%residuals are scaled to account for different magnitudes of glucose
%concentration and OD
res=[resod./nanmean(OD(ispec,:)); resgluc./nanmean(gluc(ispec,:))];
