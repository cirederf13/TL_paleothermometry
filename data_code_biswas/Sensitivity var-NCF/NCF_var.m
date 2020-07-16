% written by Rabiul Haque Biswas (biswasrabiul@gmail.com)
% this code estimates the var-NCF (estimate for initial sensitivity change) 

clear all; close all; clc;

file = xlsread('MBTP_NCF');
rawdata = file(2:end-1, 2:end);

for i=1:size(rawdata,2)       
    TL1=rawdata(:,1:2:end);
    TL2=rawdata(:,2:2:end);
    ratio =TL1./TL2;
end

temp=linspace(0,200,250)';
% plot(temp,TL2);

Temp_ncf= 90:10:150;
factor=interp1(temp,ratio,Temp_ncf);

for k=1:length(Temp_ncf)
    for i=1:size(factor,2)/3
        ncf_mean(k,i)=mean(factor(k,3*i-2:3*i));
        ncf_std(k,i)=std(factor(k,3*i-2:3*i));
    end  
end 

TL_temp=[215 225 235 245]';
temp_plot=linspace(90,300,10);
for i=1:size(factor,2)/3
    [P(i,:),S(i,:)]=polyfit(Temp_ncf',ncf_mean(:,i),1);
    ncf_var(:,i)=P(i,1)*TL_temp+P(i,2);
    ncf_plot(:,i)=P(i,1)*temp_plot+P(i,2);
end

f1=figure(1);axis square; box on; hold on
plot(temp,TL1(:,1),'Linewidth',2.0);
plot(temp,TL2(:,1),'Linewidth',2.0);
xlim([0 400]);
% ylim([0 3]);
legend('Before Natural (TL1)','After Natural (TL2)');
legend boxoff;
xlabel('TL temperature (^oC)');
ylabel('TL Intensity (counts/s)');
set(gca,'FontSize',20);


f2=figure(2); axis square; box on; hold on
col=[0    0.4470    0.7410; 0.8500    0.3250    0.0980; 0.9290, 0.6940, 0.1250; 0.4940, 0.1840, 0.5560; 0.4660, 0.6740, 0.1880; 0.3010, 0.7450, 0.9330];

for i=1:6
P(i)=plot(Temp_ncf,ncf_mean(:,i),'o','MarkerSize', 10, 'MarkerFaceColor',col(i,:),'MarkerEdgeColor','none'); 
plot(temp_plot,ncf_plot(:,i),'Color',col(i,:));
plot(TL_temp,ncf_var(:,i),'o','MarkerSize', 10, 'MarkerFaceColor','none','MarkerEdgeColor',col(i,:)); 
end
xlim([90 300]);
ylim([0 3]);
legend([P(1) P(2) P(3) P(4) P(5) P(6)],'MBTP1','MBTP2', 'MBTP5','MBTP6', 'MBTP9','MBTP11');
legend boxoff;
xlabel('TL temperature (^oC)');
ylabel('var-NCF');
set(gca,'FontSize',20);


% print(f1,'MBTP1_NCF', '-dpdf', '-r300');
% print(f2,'MBTP_var_NCF', '-dpdf', '-r300');
