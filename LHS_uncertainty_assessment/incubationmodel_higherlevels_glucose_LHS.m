% Script with Latin-Hyper-cube design to sample parameter sets for
% monocultures, select the top10 with respect to RMSE and then use these for
% simulating in-silico CO-cultures.

clear all

% set random seed
rng(42)

% define parameters
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
                1 1 1 1]; % treatments x species, 1=presence, 0=absence

% LHS design
lb = [0.1,  10,  0.001,  0.1];   % Lower bounds
ub = [1.5,  80,  0.1, 10.0];   % Upper bounds

% Number of samples
N = 10000;
% generate parameters with LHS design
X_norm = lhsdesign(N, 4);
% Scale samples to parameter bounds
X = lb + X_norm .* (ub - lb);

nrand=N; % number of randomisations
%preallocation of result arrays
mtime=zeros(500, 16, nrand); % time
modtot=zeros(500, 2, 16, nrand); % total model results

for k=2:5 % iterating through treatments --> needed because of different initial combinations (0 is control --> no fit)

    M.itreat=k; % treatmentnr.
    for m=1:nrand % iterating through random parameter values

        % generate parameter sets
        par=X(m,:)'; % select from LHS design
        [mt,mod]=incubations_glucose(par,M);
        mtime(1:length(mt),k,m)=mt;
        modtot(1:length(mod),:,k,m)=mod;
    end

end

%% evaluate against observations
load 'ODmat_corr.mat' % measured OD
load 'gluc.mat' % measured glucose concentration
tOD=[0.00	15.33	38.50	63.25	86.50	110.50	135.50	159.43	183.33...
    207.33	231.00	255.50	302.67	374.88]; % time points OD measured
tgluc=[0.00	63.25	135.50	255.50]; % time points glucose measured

for k=2:5 % iterate over 4 monocultures
    for m=1:nrand % iterate over randomisations
        ispec=k; % index for treatment in observational data
        datamod=modtot(:,1,k,m);
        datatime=mtime(:,k,m);
        datatime(datamod==0)=[];
        datamod(datamod==0)=[];

        % interpolate for comparison
        ODmod=interp1(datatime, datamod, tOD);
        glucmod=interp1(datatime, datamod, tgluc);

        rmserror(k,m)=rmse([ODmod/nanmean(OD(ispec,:)), ...
            glucmod/nanmean(gluc(ispec,:))],[OD(ispec,:)/nanmean(OD(ispec,:)), gluc(ispec,:)/nanmean(gluc(ispec,:))]); % this is scaled to mean
    end
end

% identify the indices with the 10 lowest rmserrors
for k=2:5
    data=[rmserror(k,:)', [1:1:nrand]'];
    sorted=sortrows(data,1);
    top10(:,:,k-1)=sorted;
end

% take model output from those and plot
figure
for i=1:4
    subplot(2,2,i)
    hold on
    pe=errorbar(t,OD(i+1,:),ODstd(i+1,:),'Color', [253,174,97]/255);
    p2=plot(t,OD(i+1,:),'o','MarkerFaceColor', [253,174,97]/255, 'MarkerEdgeColor','k');

    timed=squeeze(mtime(:,i+1,top10(1:10,2,i)));
    data2plot=squeeze(modtot(:,1,i+1,top10(1:10,2,i)));
    timed(data2plot==0)=NaN;
    data2plot(data2plot==0)=NaN;
    plot(timed,data2plot,'-','LineWidth',2,'Color',[53,151,143]/255)
    title(treatmentstr(i+1))
    ylim([0 0.7])
end
%% full system
% now we use the top 10 parameters sets for simulating the co-cultures
tops=10;
for i=1:16
    for k=1:tops
        indk=squeeze(top10(k,2,:)); % 4 indices that gives us the best parameter sets, sorted
        % combine them in one matrix:
        params=[X(indk(1),:); X(indk(2),:); X(indk(3),:); X(indk(4),:)]'; % sh4, sh22, sh24, sh40 / mumax, K, lambda, conv
        % run the model
        M.itreat=i; % treatmentnr.
        M.nspec=sum(M.treatmentmatrix(M.itreat,:),2); % no of species in incubation
        [mt,mod]=incubations_mat_glucose(params,M);
        %collect it, time is always the same grid
        modall(:,:,k,i)=mod; % time, species, top10 parvars, treatment
    end
end



figure
for i=1:16
    if i==1, n=1; end
    if i>1 && i<6, n=i+5; end
    if i>5 && i<16, n=i+7; end
    if i==16, n=25; end
    subplot(5,6,n)
    hold on

    data=sum(squeeze(modall(:,1:4,:,i)),2);
    data(data==0)=NaN;
    plot(mt,squeeze(data),'-','LineWidth',0.5,'Color',[.7 .7 .7])
    pe=errorbar(t,OD(i,:),ODstd(i,:),'Color', 'k');
    p2=plot(t,OD(i,:),'o','MarkerFaceColor', 'k', 'MarkerEdgeColor','k','MarkerSize', 3);
    title(treatmentstr(i))
    ylim([0 0.7])
    if i==1 || i==2 ||i==6 || i==12 || i==16
        ylabel('OD_{600}')
    end
    if i==10 || i==11 || i==13 || i==14 || i==15 || i==16
        xlabel('time [h]')
    end
end
