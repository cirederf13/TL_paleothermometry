% written by Rabiul Haque Biswas (biswasrabiul@gmail.com)
% this code estimates the thermal decay parameters (E,s,and b)

clear all; close all

% Export data from Excel file
sample='MBTP1';
filename = [sample '_TmTstop'];
rawdata=xlsread(filename);

% Read TL data
Fracglow = rawdata(2:1:end-1,5:2:end);
TD=rawdata(2:1:end-1,6:2:end);

% Glow curve normalization
for i=1:length(TD(1,:))
    maxTD(i)=max(TD(:,i));
    normTL(:,i)=Fracglow(:,i)./max(TD(:,i));
end

k=8.617*10^-5;     % Boltzman constant

% Read temperature
T1=273+1:273+450;
T=T1';

%%%%%% Fitting subpeak %%%%%%%%%
deltaTm=5;
deltaE=0.015;
for i=10:size(normTL,2)- 4                    %%%%% 2 curve less
    subpeak(:,i)=normTL(:,i)-normTL(:,4+i);  %%% 1 minus third
    y(:,i)=subpeak(:,i)./max(subpeak(:,i));
    y(y<=0)=0;
    params(:,i)=[1 475+(i-30)*deltaTm 1.34+(i-30)*deltaE 1.6];
    [beta(:,i),R,J,Cov,MSE]=nlinfit(T,y(:,i),@single_peak,params(:,i));
    sbeta(:,:,i)=nlparci(beta(:,i),R,'covar',Cov,'alpha',1-.95);
    sigma(:,i)=sbeta(:,2,i)-beta(:,i);          %calculate error
end

%%%%%% Representative curve of subpeak and fitting %%%%%%
f1=figure(1);
pbaspect([2 1 1]);
% axis square; 
box on; hold on
for i=1:2:size(normTL,2)-4
    hold on
    A(:,i)=y(:,i)*max(subpeak(:,i));
    B(:,i)=single_peak(beta(:,i),T)*max(subpeak(:,i));
    plot(T-273,A(:,i),'o','markersize',1.0);
    plot(T-273,B(:,i),'LineWidth',0.5);
end
xlim([0 500]);
xlabel('Temperature (\circC)');
ylabel('Subpeak TL Intensity (a.u.)');
% set(gca,'Yscale','log');
% ylim([0.001 2]);
set(gca,'FontSize',11);
% ax.LineWidth = 2.0;


%%%%%%%% Parameter Extraction %%%%%%%
Im=(beta(1,1:end).*max(subpeak))'; errIm=(sigma(1,1:end).*max(subpeak))';
E=beta(3,1:end)';   errE=sigma(3,1:end)';
Tm=beta(2,1:end)';  errTm=sigma(2,1:end)';
b=beta(4,1:end)';   errb=sigma(4,1:end)';
HR=1;     %%%%% Heating rate

%%%%%% Plotting of Activation Energy vs. subpeak temperature %%%%%
f2=figure(2); 
% subplot(2,2,1);
box on; axis square; hold on
errorbar(Tm-273,E,errE,'o','Color',[1 0.5 0],'MarkerSize',8,'MarkerFaceColor',[0.25 0.25 0.25],'MarkerEdgeColor','none');

%%%% Smoothing of curve%%%%
yE = smooth(Tm,E,0.95,'rloess');
[xE,ind] = sort(Tm);
plot(xE-273,yE(ind),'LineWidth', 0.5);

xlabel('Peak Temperature (\circC)');
ylabel('E (eV)');
xlim([0 500]);
ylim([0.5 2.5]);
set(gca,'FontSize',20);
ax.LineWidth = 2.0;
% gtext('a)');
% legend('Subpeak Fitting', 'Smooth Fit', 'Location', 'northwest');
% legend boxoff;

%%%%% Calculation of frequency factor and life time %%%%%
for i=1:length(E)
s(i)=HR*E(i)/k/Tm(i)^2/(1+(b(i)-1)*2*k*Tm(i)/E(i))*exp(E(i)/k/Tm(i));
errs(i)=s(i)*sqrt((errE(i)/E(i))^2+(errTm(i)/Tm(i))^2+(errb(i)/b(i))^2);
s10(i)=log10(s(i)); errs10(i)=s10(i)*errs(i)/s(i);
tau(i)=s(i)^-1*exp(E(i)/k/288)/3600/24/365;
tau10(i)=log10(tau(i));
errtau(i)=tau(i)*sqrt((errs(i)/s(i))^2+(errE(i)/E(i))^2);
errtau10(i)=tau10(i)*errtau(i)/tau(i);
end;

%%%%%% Plotting of Freqency factor vs. subpeak temperature %%%%%
f3=figure(3); 
% subplot(2,2,2);
box on; axis square; hold on;
errorbar(Tm-273,s10,errs10,'o','Color',[1 0.5 0],'MarkerSize',8,'MarkerFaceColor',[0.25 0.25 0.25],'MarkerEdgeColor','none');

%%%% Smoothing of curve%%%%
ys10 = smooth(Tm,s10,0.95,'rloess');
[xs10,ind] = sort(Tm);
plot(xs10-273,ys10(ind),'LineWidth', 0.5);

xlabel('Peak Temperature (\circC)');
ylabel('log_{10}(s (s^{-1})) ');
xlim([0 500]); ylim([0 20]);
set(gca,'FontSize',20);
ax.LineWidth = 2.0;
% gtext('b)');

% legend('Subpeak Fitting', 'Smooth Fit', 'Location', 'northwest');
% legend boxoff;

%%%%%% Plotting of b vs. subpeak temperature %%%%%
f4=figure(4);
% subplot(2,2,3)
axis square; box on; hold on
errorbar(Tm-273, b,errb, 'o','Color',[1 0.5 0],'MarkerSize',8,'MarkerFaceColor',[0.25 0.25 0.25],'MarkerEdgeColor','none');

%%%% Smoothing of curve%%%%
yb = smooth(Tm,b,0.95,'rloess');
[xb,ind] = sort(Tm);
plot(xb-273,yb(ind),'LineWidth', 0.5);

xlabel('Peak Temperature (\circC)');
ylabel('b'); 
xlim([0 500]); ylim([0.0 3.0]);
set(gca,'FontSize',20);
ax.LineWidth = 2.0;
% gtext('c)');

% legend('Subpeak Fitting', 'Smooth Fit', 'Location', 'northwest');
% legend boxoff;


%%%%%% Plotting of Life time vs. subpeak temperature %%%%%
f5=figure(5); 
% subplot(2,2,4)
axis square; box on; hold on
errorbar(Tm-273,tau10,errtau10,'o','Color',[1 0.5 0],'MarkerSize',8,'MarkerFaceColor',[0.25 0.25 0.25],'MarkerEdgeColor','none');
xlabel('Peak Temperature (\circC)');
ylabel('Life time (year)');,
xlabel('Peak Temperature (\circC)');
ylabel('log_{10}(\tau (y))');
xlim([0 500]);
ylim([-5 20]);
set(gca,'FontSize',20);
ax.LineWidth = 2.0;
% gtext('d)');



%%% Literature data %%%%%%
litx=[90 110 210 280 320 350];
lity1=[0.16*10^-3 43*10^-3 3.6*10^3 3.9*10^6 1*10^9 9.2*10^9];
lity=log10(lity1);
plot(litx,lity,'o','MarkerSize',8,'MarkerFaceColor','r','MarkerEdgeColor', 'none')
legend('Subpeak Fitting', 'Literature', 'Location', 'northwest');
legend boxoff;


% Parameter extraction from smooth fitting data
Temp=[215 225 235 245];

for j=1:length(Temp)
for i=1:length(xE)
    L(i) = abs(273+Temp(j) - xE(i));
end
[M,I(j)]=min(L);
avgE(j)=yE(I(j));
avgs10(j)=ys10(I(j));
avgb(j)=yb(I(j));
end

parama_decay=[avgE' avgs10' avgb'];
save ([sample '_param_decay'], 'parama_decay')



print(f2,['E_' sample], '-dpdf', '-r1000');
print(f3,['s_' sample], '-dpdf', '-r1000');
print(f4,['b_' sample], '-dpdf', '-r1000');

