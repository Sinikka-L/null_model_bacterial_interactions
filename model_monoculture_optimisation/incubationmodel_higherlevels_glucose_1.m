% script for using the parameter sets from mono-cultures for co-cuture
% simulations

clear all

% best aprameter sets from previous fits for each mono-culture:
mumax=[0.9747, 0.27, 0.3135, 0.6938];
K=[47.3381, 65.01, 23.81, 56.14];
lambda=[0.0056, 0.0053, 0.0242, 0.0235];
conv=[4.2319, 2.8419, 0.8168, 0.3017];
M.treatmentmatrix=[0 0 0 0;...
                1 0 0 0;...
                0 1 0 0;...
                0 0 1 0;...
                0 0 0 1;...
                1 1 0 0;...
                1 0 1 0;...
                1 0 0 1;...
                0 1 1 0;...
                0 1 0 1;...
                0 0 1 1;...
                1 1 1 0;...
                1 1 0 1;...
                1 0 1 1;...
                0 1 1 1;...
                1 1 1 1];
par=[mumax;...
    K;...
    lambda;...
    conv];

mtime=zeros(500,16);
modtot=zeros(500,5,16);
for k=2:16

    M.itreat=k; % treatmentnr.
    M.nspec=sum(M.treatmentmatrix(M.itreat,:),2); % no of species in incubation
    [mt,mod]=incubations_mat_glucose(par,M);
    mtime(1:length(mt),k)=mt;
    modtot(1:length(mod),:,k)=mod;
end

load 'ODmat_corr.mat'
figure
for i=1:16
    if i==1, n=1; end
    if i>1 && i<6, n=i+5; end
    if i>5 && i<16, n=i+7; end
    if i==16, n=25; end
    subplot(5,6,n)
    hold on
    pe=errorbar(t,OD(i,:),ODstd(i,:),'Color', [253,174,97]/255);
    p2=plot(t,OD(i,:),'o','MarkerFaceColor', [253,174,97]/255, 'MarkerEdgeColor','k');
    timed=mtime(:,i);
    data=sum(squeeze(modtot(:,1:4,i)),2);
    timed(data==0)=NaN;
    data(data==0)=NaN;
    plot(timed,data,'-','LineWidth',2,'Color',[53,151,143]/255)
    title(treatmentstr(i))
    ylim([0 0.7])
    if i==1 || i==2 ||i==6 || i==12 || i==16
        ylabel('OD')
    end
end
