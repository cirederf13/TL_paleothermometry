% written by Rabiul Haque Biswas (biswasrabiul@gmail.com)
% This is the code for inverse modelling

clear all; close all; clc;
tic
%% Load TL parameters
sample ='MBTP1';
filename = ['Summary_Parameter_GOK_' sample];
rawdata=xlsread(filename);

% pool=parpool(10);

kparams=rawdata(:,1:2:18); skparams=rawdata(:,2:2:18);
nN_mean=rawdata(:,23); sigmanN=rawdata(:,24);
KarsnN=rawdata(:,25); sigmaKarsnN=rawdata(:,26);

%% Defining Parameters
timemax=1.0;        % Time in Ma
timemin=0;          % Time in Ma
nstep=2000;         % this is to reinterpolate the Ttpath for the nN model
niter=3000;         % number of random paths, fixed to equal the number of random realisations

Dtime   =timemax/(nstep-1);
timeM   =0:Dtime:timemax;

% temp define from delta 18O data
data = xlsread('Greenland_delta18O');
time =data(:,1)/10^6;            % time in Ma
detal18O=data(:,2);

% time_scaled    = linspace(0, 1, 2000);      % In ka
detal18O_f =fliplr(interp1(time,detal18O,timeM));

%scaling
[val,idx]=min(abs(timeM-0.983)); %base point at 17 ka
max_temp=detal18O_f(end);
min_temp=detal18O_f(idx);

detal18O_norm=(-min_temp+detal18O_f)./(max_temp-min_temp);


parfor i=1:niter
residuals=zeros(size(kparams,1),1);
nNf=zeros(nstep,size(kparams,1),niter);
Predicted_nN=zeros(niter,size(kparams,1));
nN=zeros(size(kparams,1),1);


% Generating Random oscillation path 
amp=30*rand;
meanT=40*rand-20;
tempM=meanT+detal18O_norm.*amp;


   
plot(timeM, tempM); hold on
    xlim([max(timeM)-100/1000,max(timeM)]);                  % Plotting of last 200 ka
% ylim([-10,30]);
    xlabel('time (ka)');
    ylabel('temp (^oC)');
    set(gca,'FontSize',15);

%% Calculate nN for each Tt paths
 
    for k=1:size(kparams,1)
        nNf(:,k,i) = TLModel_GOK_MCerror(timeM,tempM,kparams(k,:),skparams(k,:));
        nN(k)=nN_mean(k);
        residuals(k)=((nN(k)/sigmanN(k))*0.5.*log(nN(k)/nNf(end,k,i))).^2;
        Predicted_nN(:,k)=nNf(end,k,i);    
    end
    sum_residual=sum(residuals)/length(nN);
    
    misfit(i,1)=sum_residual;
    TPath(i,:)=tempM;
    nNPred(i,:) = Predicted_nN(i,:);
    
    Tbase(i,1)=meanT;
    Tamp(i,1)=amp;
    
    

fprintf('Path%i   \n ',i);
end

Out_misfit = misfit; 
Out_TPath = TPath;
Out_nNPred = nNPred;
Out_Tbase = Tbase;
Out_Tamp = Tamp;


% Sorting the data
disp('sorting table')
[sortedmisOUT,IX] = sort(Out_misfit(:,1));
sortedTt = Out_TPath(IX,:);
sortedPrednN = Out_nNPred(IX,:);
sortedTbase=Out_Tbase(IX,1);
sortedTamp=Out_Tamp(IX,1);


% Saving the data
save([sample '_Tt.mat'], 'sortedTt', '-v7.3');
save ([sample '_misfit.mat'], 'sortedmisOUT', '-v7.3');
save ([sample '_PrednN.mat'], 'sortedPrednN', '-v7.3');
save ([sample '_time.mat'], 'timeM', '-v7.3')
save ([sample '_Tbase.mat'], 'sortedTbase', '-v7.3')
save ([sample '_Tamp.mat'], 'sortedTamp', '-v7.3')

% delete(pool);
toc
