% written by Rabiul Haque Biswas (biswasrabiul@gmail.com)
% this code estimates the atherthermal decay parameter (rhoprime), and
% Growth parameter (D0 and a) and natural TL (n/N)


clear all; close all; clc
tic

sample = 'MBTP1'
ddot=7.39;                      %dose rate in Gy/ka
SystemDoseRate=0.236;           % in Gy/s Risoe3

for jr=1:4
%% Fading

% Export data from Excel file
filename = [sample '_Fading'];

rawdata1=xlsread(filename);
TLfading=rawdata1(2:end-1,2:end);
% n=input('enter disc number>');
fad_TL=TLfading(:,1:2:end);
TD_TL=TLfading(:,2:2:end);


% Signal Integration
signal1=210+(jr-1)*10;
signal2=220+(jr-1)*10;


TL1 = zeros(size(fad_TL,2),1);
TL2 = zeros(size(TD_TL,2),1);

for i=1:1:size(fad_TL,2)
    sum1 = 0; sum2=0;
    for j=signal1:signal2          
    sum1 = sum1 + fad_TL(j,i);
    sum2 = sum2 + TD_TL(j,i);
    TL1(i) = sum1;TL2(i)=sum2;
    end
end

fady=TL1./TL2; 

filename = [sample '_Fdelay'];
tdelay1=xlsread(filename);
tdelay = reshape(tdelay1',[],1);      % delay time in hour
fadx=tdelay*3600;


%%nlinfit solution rhop%%
    beta01=[1 -5]; %%initial parameters, a=scaling parameter, b=rhop estimation%%
    [beta1,R,J,Cov,MSE]=nlinfit(fadx,fady,@Ffitrho,beta01);

    %Confidence interval rhop%
    fvec=linspace(1,1e8,1e3);
    [Ypred1,delta1]=nlpredci(@Ffitrho,fvec,beta1,R,'covar',Cov,'alpha',.05,'predopt','curve');
    sbeta1=nlparci(beta1,R,'covar',Cov,'alpha',1-.68); 
    sigma1=beta1'-sbeta1(:,2); %calculate error on rhop
    
    rho10=beta1(2); errrho10=sigma1(2);
    
    
    %%g value calculation%%
    beta02=[1 -0.05]; %%initial parameters, I, m%%
    fy=fady./fady(1); fx=fadx./172800;
    [beta2,R,J,Cov,MSE]=nlinfit(fx,fady,@Ffitg,beta02); g2d=(-100*beta2(2));
    
    %Confidence interval g2days%
    [Ypred2,delta2]=nlpredci(@Ffitg,fvec,beta2,R,'covar',Cov,'alpha',.05,'predopt','curve');
    sbeta2=nlparci(beta2,R,'covar',Cov,'alpha',1-.68); sigma2=beta2'-sbeta2(:,2); %calculate error on m for g2days
    
    g=-beta2(2)*100;errg=-sigma2(2)*100;
  
    %%print figure of fit      
    f(jr)=figure(jr); 
    subplot(1,2,1)
    d=cat(2, Ypred1-delta1, fliplr(Ypred1+delta1));
    xx=cat(2, fvec, fliplr(fvec));
    semilogx(fvec,Ypred1); hold on; axis square; box on;
    fill([xx], [d], [0.5 0.5 0.5], 'edgecolor','none'); 
    semilogx(fvec,Ypred1,'Color',[0.8500    0.3250    0.0980],'LineWidth',1.0);
    plot(fadx,fady,'o','MarkerSize', 8, 'MarkerEdgeColor','none','MarkerFaceColor','k');
    ylim([0.7 1.05]); xlim([5*1e1 1e8]); xlabel('delay time (s)'); ylabel('TL Intensity (a.u.)');
    ax = gca; set(ax,'XTick',[1e2, 1e3, 1e4, 1e5, 1e6, 1e7, 1e8]);
    set(ax,'YTick',[0.7, 0.8, 0.9, 1.0]); box on
    set(gca,'FontSize',12);
    
     %Add params to figures%
    rhoplab=sprintf('log_1_0(\\rho′)= %0.2f ± %0.2f', rho10, errrho10); 
%     glab=sprintf('g_{2days}= %0.2f �� %0.2f', g, errg);
    text(1.5e2,0.745,rhoplab,'FontSize', 10);
%     text(1.5e2,0.765,glab,'FontSize', 10);
    

%% Growth
    
    rhop10=beta1(2);
    
    
 % Export data from Excel file
% filename = [sample '_MAR'];
A = xlsread([sample '_MAR_PG']);
rawdata=A(2:end-1,2:end);

% ncf1=rawdata(1:250,2:8:end); peak1=ncf1(150,:);
% ncf2=rawdata(1:250,5:8:end); peak2=ncf2(150,:);
% ratio =peak1./peak2;

nat_TL =rawdata(:,2:6:end);
reg_TL =rawdata(:,4:6:end);
TD_TL = rawdata(:,6:6:end),


N =sum(nat_TL(signal1:signal2,:));
R =sum(reg_TL(signal1:signal2,:));
TD =sum(TD_TL(signal1:signal2,:));

NTL =N./TD;
NAvg=mean(NTL);
errNavg=std(NTL);

LxTx =R./TD;
% LxTx0=mean(LxTxplusTD(1:3));
Npt=0.2;

% NAvg = NAvgplusTD-LxTx0;
% LxTx = LxTxplusTD -LxTx0;

Ds=[0,0,0,100,100,100, 200,200,200,500,500,500,1000,1000,1000,2000,2000,2000,4000,4000,4000,8000,8000,8000];   %dose in s
D=Ds*SystemDoseRate;                       % dose in Gy

rhop=10.^rhop10;
fcorrecLxTx=LxTx./exp(-rhop*log(1.8*3e15.*(0.5*Ds)).^3);
fcorrecLxTx(isnan(fcorrecLxTx))=0;


% Estimation of a parameter
lb=[0,0, 1.0];
ub=[200,3000,2.0];
x0=[150 1000 1.5];
fun = @(x)x(1).*(1-(1+(D./x(2)).*(x(3)-1)).^(-1/(x(3)-1)))-fcorrecLxTx;
[params,resnorm,residual,exitflag,output,lambda,jacobian] = lsqnonlin(fun,x0,lb,ub);
jacobian = full(jacobian);
N=length(fcorrecLxTx);
Cov = resnorm*inv(jacobian'*jacobian)/N;
sparams = sqrt(diag(Cov));
global a;
a=params(3); erra=sparams(3);

%%nlinfit solution%%
beta0=[50 500]; %%initial parameters, a=scaling parameter, b=D0 estimation%%
[beta,R,J,Cov,MSE]=nlinfit(D,fcorrecLxTx,@SARGOK,beta0);
modvec=linspace(1,1e6,1e5); %xvecLn=ones(size(LnLmax)).*0.2;

%Confidence intervals%
[Ypred,delta]=nlpredci(@SARGOK,modvec,beta,R,'covar',Cov,'alpha',.05,'predopt','curve');
sbeta=nlparci(beta,R,'covar',Cov,'alpha',1-.68);
sigma=sbeta(:,2)-beta'; %calculate error on D0



%%Normalise data and model for figure
rawy=LxTx/beta(1); yplot=fcorrecLxTx./beta(1); Ypredplot=Ypred./beta(1); deltaplot=delta./beta(1); Lnplot=NAvg./beta(1);

D0=beta(2);errD0=sigma(2);
naty=Lnplot; errnaty=std(NTL)/beta(1);

%%figure

subplot(1,2,2);

%plot CI
d=cat(2, Ypredplot-deltaplot, fliplr(Ypredplot+deltaplot)); xx=cat(2, modvec, fliplr(modvec));
d(isnan(d))=0.00001; 
semilogx(modvec,Ypredplot); axis square; box on; hold on;
fill(xx,d,[0.5 0.5 0.5],'edgecolor','none'); 
semilogx(modvec,Ypredplot,'Color',[0.8500    0.3250    0.0980],'LineWidth',1.0); 

%plot data points
fadx(fadx==0)=0.2; %remove zero values to faciliate plotting
X1=plot(D,rawy,'o','MarkerSize', 10,'MarkerEdgeColor','none','MarkerFaceColor',[0.7 0.7 0.7]);
X2=plot(D,yplot,'o','MarkerSize', 10,'MarkerEdgeColor','none','MarkerFaceColor','k');
X3=plot(Npt,Lnplot,'s','MarkerSize', 12,'MarkerEdgeColor','none','MarkerFaceColor','g');
legend ([X1 X2 X3],'Raw data', 'Fading corrected data','Natural point','Location', 'northwest');
legend boxoff
% semilogx(modvec,UFYpredplot,'k-');
ylim([0 1.1]); xlim([0.1 1e4]); ax = gca; set(ax,'XTick',[1e0, 1e1, 1e2, 1e3, 1e4, 1e5]);
xlabel('Dose (Gy)'); ylabel('n'); box on
set(gca,'FontSize',12);


%%Calculate Kars model%%
kb = 8.617343*10^(-5); alpha = 1; Hs=3E15; %s value after Huntley (2006) J. Phys. D.
rprime=linspace(0.01,2.5,100)'; [rpm,rpn]=size(rprime); %create vector of rprime distances

% rhop=0; rhop(rhop<0)=0;
rhop=10.^rho10; rhop(rhop<0)=0;
rhopub=10.^(rho10+errrho10); rhopub(rhopub<0)=0;
rhoplb=10.^(rho10-errrho10); rhoplb(rhoplb<0)=0;
rhopbound=[rhop rhopub rhoplb];
for j=1:length(rhopbound)
rho = 3*alpha^3*rhopbound(j)/(4*3.1415);
r = rprime./(4*3.1415*rho/3)^(1/3);
pr = 3.*rprime.^2.*exp(-rprime.^3.);
tau=((1/Hs)*2.71.^(alpha.*r))/(3600*24*365*1000);
Ls=1./(1+D0./(ddot.*tau));
Lstrap=(pr.*Ls)/sum(pr);
Karsres(j)=sum(Lstrap);
end
% end

KarsAv=(Karsres(1)); errKarsAv=max([(abs(Karsres(1)-Karsres(2))) (abs(Karsres(1)-Karsres(3)))],[],2);


%Add params to figures%
    D0lab=sprintf('D_0= %0.0f ± %0.0f Gy' , D0, errD0);
    growthpara=sprintf('a=%0.2f ± %0.2f', a, erra);
    ntllab=sprintf('n_{obs}= %0.2f ± %0.2f', naty, errnaty);
%     Karslab=sprintf('Nss/N= %0.2f �� %0.2f', KarsAv, errKarsAv);
    text(0.25,0.6,D0lab,'FontSize', 10);
    text(0.25,0.5,growthpara,'FontSize', 10);
    text(0.25,0.4,ntllab,'FontSize', 10);
%     text(0.25,0.45,Karslab,'FontSize', 18);
    
    % Parameter extraction
    param_growth(jr,1)=D0;
    param_growth(jr,2)=errD0;
    param_growth(jr,3)=a;
    param_growth(jr,4)=erra;
    param_fading(jr,1)=rho10;
    param_fading(jr,2)=errrho10;
    param_fading(jr,3)=g;
    param_fading(jr,4)=errg;
    param_natTL(jr,1)=naty;
    param_natTL(jr,2)=errnaty;
    param_natTL(jr,3)=KarsAv;
    param_natTL(jr,4)=errKarsAv;
    
end
%%
save([sample '_param_growth'],'param_growth');
save([sample '_param_fading'],'param_fading');
save([sample '_param_natTL'],'param_natTL');


toc


