% written by Rabiul Haque Biswas (biswasrabiul@gmail.com)
% This code is for PDF plot and temperature extranction

tic

clear all;clc;clf,close all,
set(gcf,'units','centimeters')
set(gca,'FontSize',10)
set(0,'DefaultAxesFontName','Helvetica');
set(0,'DefaultTextFontName','Helvetica');

% Loading the file
sample='MBTP1';
for l=1:10               % iteration to stabilize the solution 
load([sample '_Tt.mat']); Tt=sortedTt;
load([sample '_misfit.mat']); misOUT=sortedmisOUT;
% load([filename '_nN.mat']); nNs=sortednNsave; 
load([sample '_PrednN.mat']); PrednNs=sortedPrednN; 
load([sample '_time.mat']); time=timeM;

load([sample '_Tbase.mat']); Tbase=sortedTbase;
load([sample '_Tamp.mat']); Tamp=sortedTamp;

rawdata=xlsread(['Summary_Parameter_GOK_' sample]);
nN=rawdata(:,23); sigmanN=rawdata(:,24);
KarsnN=rawdata(:,25); sigmaKarsnN=rawdata(:,26);

MTemp=[215 225 235 245];

% Grid define
nAvM=1000;
[m,nt]=size(PrednNs);

% to compute the PDF for the likelhood, the Zt paths are reinterpolated onto a grid (Av_matrix)
time_max=1; time_min=0; dt=(time_max-time_min)/(nAvM-1);
T_max=500; T_min=-500; dT=(T_max-T_min)/(nAvM-1);
vec_time=time_min:dt:time_max;
Tvec=T_min:dT:T_max;

% Data selection %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
R = rand(m,1); 
prob=exp(-misOUT); scale=max(prob); % scaling for resampling the PDF
% scale=1;
test=prob/scale>R;


idefix=find(test);
idefix=idefix(end:-1:1); % goes from largest misfit to smallest
movea=length(idefix);

max_likelihood=max(prob(idefix)); min_likelihood=min(prob(idefix)); % bounds within the max and min you want the color scheme to be
index_color=max(1,floor(63*(prob-min_likelihood)/(max_likelihood-min_likelihood))+1); % scaling the value between 1 and 64

MTEMP=ones(movea,1)*MTemp(1:nt);
INDEX=index_color(idefix)*ones(1,nt);
PrednNs=PrednNs(idefix,:);

Av_matrix=zeros(nAvM);
% add accepted model to a matrix to compute the PDF, this is based on a rejection algorithm
for k = idefix';
	vec_T=interp1(time,Tt(k,:),vec_time,'linear');
    Tpath=(0:nAvM-1)*nAvM+round((vec_T-T_min)/dT)+1;
    Tpath(Tpath<=0)=[];
	Av_matrix(Tpath)=Av_matrix(Tpath)+1;
end
X=cumsum(Av_matrix/movea); % for computing CIs and median


% Plot figures %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
map=colormap(parula); % creates a table with 64 rows and three columns for RGB
col=round(100*(max(misOUT(idefix)):(min(misOUT)-max(misOUT(idefix)))/6:min(misOUT)))/100; % create colorbar for misfit



% Figure 1: Predicted vs. observed nN
f1=figure(1); axis square; box on; hold on
xlabel('TL temperature (^oC)');
ylabel('n');
c1=colorbar('yticklabel',col);
set(get(c1,'title'),'string','misfit');
% colormap(parula);

median_path=scatter(MTEMP(:),PrednNs(:),60,INDEX(:),'filled');   % predicted nNs
actual=errorbar(MTemp(1:nt),nN(1:nt),sigmanN(1:nt),'ko','MarkerSize', 20);   % observed nNs
legend([actual,median_path],'Observed','Predicted', 'Location', 'NorthWest'); legend boxoff;

xlim([200 250]);
ylim([1e-2 10]);
set(gca,'FontSize',20);
ax.LineWidth = 2.0;
set(gca,'Yscale','log');


% Figure 2: Spaghetti plot of selected paths
f2=figure(2); axis square; box on; hold on
xlabel('Time (Ma)');
ylabel('Temp (^oC)');
c1=colorbar('yticklabel',col);colormap(jet);
set(get(c1,'title'),'string','misfit');
P = plot(time,Tt(idefix,:)); % 'spaghetti' plot
for curve = 1:movea
	set(P(curve),'Color',map(index_color(idefix(curve)),:))
end

xlim([time_max-0.1 time_max]);
ylim([-20 20])
set(gca,'FontSize',20);
ax.LineWidth = 2.0;

% Figure 3: Probability Density Plot
f3=figure(3); axis square; box on; hold on
xlabel('Time (ka)');
ylabel('Temperature (^oC)');
axis([time_max-0.05 time_max -20 20],'square');
cc=othercolor('Reds9',100)
c2=colorbar;
set(get(c2,'title'),'string','PDF');
hold on
contourf(vec_time,Tvec,Av_matrix./movea,100,'edgecolor','none'); % PDF
shading flat
twosig=contour(vec_time,Tvec,X,[0.05,0.95],'k','LineWidth',0.5); % 95 CI
onesig=contour(vec_time,Tvec,X,[0.32,0.68],'w','LineWidth',1.5); % 68 CI
median_path=contour(vec_time,Tvec,X,1,'r','LineWidth',2); % median

caxis([0 0.1]); 
set(gca,'FontSize',20);
ax.LineWidth = 2.0;
set(gca,'XTickLabel', {'50','40','30','20','10','0'});
% set(gca,'xaxisLocation','top');
% set(gca,'yaxisLocation','right');
% colorbar('westoutside')
% set(gca,'caxisLocation','west');


Tbase_final=Tbase(idefix,1);            n_Tbase=round(max(Tbase_final)-min(Tbase_final));
Tamp_final=Tamp(idefix,1);              n_Tamp=round(max(Tamp_final)-min(Tamp_final));
Tsurf_final=Tbase_final+Tamp_final;     n_Tsurf=round(max(Tsurf_final)-min(Tsurf_final));


f4=figure(4); axis square; box on; hold on
histogram(Tbase_final,n_Tbase,'FaceColor',[0.4 0.4 0.4]);
[hi cx]=hist(Tbase_final,n_Tbase);
AA=median(Tbase_final);
bar(AA,max(hi));
xlabel('Temp at 17 ka (^oC)');
ylabel('Frequency');
xlim([-20 30]);
set(gca,'FontSize',20);
ax.LineWidth = 2.0;

f5=figure(5); axis square; box on; hold on
histogram(Tsurf_final,n_Tsurf,'FaceColor',[0.2 0.2 0.2]);
[hi cx]=hist(Tsurf_final,n_Tsurf);
BB=median(Tsurf_final);
bar(BB,max(hi));
xlabel('Temp at present (^oC)');
ylabel('Frequency');
xlim([-20 30]);
set(gca,'FontSize',20);
ax.LineWidth = 2.0;


%% temp extraction from PDF
% Median and one sigma path
outliers = find (median_path(2,:) > 100);
median_path(:,outliers)=[];
[~,idx]=find(median_path(1,:)==1);
median_path=median_path(:,1:idx);

outliers = find (onesig(2,:) > 100);
onesig(:,outliers)=[];

[~,idx]=find(onesig(1,:)==1);
onesig_minus=onesig(:,1:idx(1));
onesig_plus=onesig(:,1+idx(1):end);

% Present day temp
T_present_median(l)=median_path(2,end);
T_present_plus(l)=onesig_plus(2,end);
T_present_minus(l)=onesig_minus(2,end);


% LGM temp 
LGM1=1-17/1000;

[~,idx]=min(abs(median_path(1,:)-LGM1)); T_LGM1_median(l)=median_path(2,idx);
[~,idx]=min(abs(onesig_plus(1,:)-LGM1)); T_LGM1_plus(l)=onesig_plus(2,idx);
[~,idx]=min(abs(onesig_minus(1,:)-LGM1)); T_LGM1_minus(l)=onesig_minus(2,idx);


T_extract(l,1)=T_present_median(l);
T_extract(l,2)=T_present_plus(l);
T_extract(l,3)=T_present_minus(l);
T_extract(l,4)=T_LGM1_median(l);
T_extract(l,5)=T_LGM1_plus(l);
T_extract(l,6)=T_LGM1_minus(l);

end

T_extract(l+2,:)=mean(T_extract);


print(f1,[sample '_nN_prediction'], '-dpdf', '-r300');
print(f3,[sample '_probPlot'], '-dpdf', '-r300');
print(f4,[sample '_Tbase_17ka'], '-dpdf', '-r300');
print(f5,[sample '_Tpresent'], '-dpdf', '-r300');

