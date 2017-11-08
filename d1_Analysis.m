%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   d1_Analysis.m                                                              %
%                                                                              %
%   Leiden University College &                                                %
%   Institute of environmental sciences (CML), Leiden University               %
%                                                                              %
%   Run input-output analysis. CBA of average diet, NRD, NRD isocaloric        %
%   - Plot Stacked bars                                                        %
%   - Plot Maps                                                                %
%   - Plot regression between P and N                                          %
%                                                                              %
%   Paul Behrens: p.a.behrens@luc.leidenuniv.nl                                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 0. Set-Up                                                                    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all; clc;
pre_process_data = true;

% run scripts if needed
if pre_process_data == true
    %a1_ProcessFAO;
    %a2_ImportExio2011;
    a3_ProcessNRD;
    c1_MakeDietVectors;
    c2_MakeDietMatrix;
    c3_MakePrices;
    %b1_CompareNRDaverage;  % note: intermediate analysis - run after the others
end

addpath('./0. article_functions/');
Supplementary_Data_File = '../Data/Supplementary_Data.xlsx';
Characterisation_File = '~/Source/EXIOBASE/characterisation_DESIRE_version3.3.xlsx';

fig_str_LQ = '../Figures_LQ/';
fig_str_PR = '../Figures_PR/';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1. Load data                                                                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% load EXIOBASE
load('~/Source/EXIOBASE/2011/sys_pxp.mat');
exio_idx = xlsread(Supplementary_Data_File,'Diet_Classifications','P3:P19');
grps = xlsread(Supplementary_Data_File,'Diet_Classifications','Q3:Q19');
[~, ~, grps_str] = xlsread(Supplementary_Data_File,'Diet_Classifications','O22:O27');

% load problem classes
GWP = xlsread(Characterisation_File,'Q_emissions','E9:PE9');
Land = xlsread(Characterisation_File,'Q_resources','C4:P4');
Eutrophication = xlsread(Characterisation_File,'Q_emissions','E55:PE55');
Eutrophication_N2 = xlsread(Characterisation_File,'Q_emissions','E56:PE56');

% load ancillary country data at this stage
pop = xlsread(Supplementary_Data_File,'Country_Classifications','D3:D46');

% load diets and prices
load('../Data/Diets.mat');
load('../Data/Prices.mat');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2. Calculate environmental impacts related to diets                          %                                                                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Format data
% multiply diets by prices
Prices = repmat(Prices,49,1)./10^6;
AD = AD.*Prices;
ND = ND.*Prices;
NCD = NCD.*Prices;

%% pre-calculate for each impact
L(1,1:9800) = GWP*io.rc*io.L;
L(2,1:9800) = 1000*Eutrophication*io.rc*io.L;
L(3,1:9800) = 365.*100.*Land*io.kc*io.L;
L_N2 = 1000*Eutrophication_N2*io.rc*io.L;

%% calculate consumption-based impacts for GHG, Eutrophication, and Land
for i = 1:1:44
    impactsAD(1:3,1:9800,i) = L*diag(AD(:,i));
    impactsND(1:3,1:9800,i) = L*diag(ND(:,i));
    impactsNCD(1:3,1:9800,i) = L*diag(NCD(:,i));
    L_N2_AD(1:9800,i) = L_N2*diag(AD(:,i));
    L_N2_NCD(1:9800,i) = L_N2*diag(NCD(:,i));
end

impactsAD = reshape(impactsAD,3,200,49,[]);
impactsND = reshape(impactsND,3,200,49,[]);
impactsNCD = reshape(impactsNCD,3,200,49,[]);
L_N2_AD = reshape(L_N2_AD,200,49,[]);
L_N2_NCD = reshape(L_N2_NCD,200,49,[]);

totAD = squeeze(sum(impactsAD,3));
totND = squeeze(sum(impactsND,3));
totNCD = squeeze(sum(impactsNCD,3));
tot_N2_AD = squeeze(sum(L_N2_AD,2));
tot_N2_NCD = squeeze(sum(L_N2_NCD,2));

totAD = totAD(:,exio_idx,:);
totND = totND(:,exio_idx,:);
totNCD = totNCD(:,exio_idx,:);
tot_N2_AD = tot_N2_AD(exio_idx,:);
tot_N2_NCD = tot_N2_NCD(exio_idx,:);

totAD = funAccumarray(grps,totAD);
totND = funAccumarray(grps,totND);
totNCD = funAccumarray(grps,totNCD);
tot_N2_AD = funAccumarray(grps,tot_N2_AD);
tot_N2_NCD = funAccumarray(grps,tot_N2_NCD);

diff = totND-totAD;
diff = permute(diff,[3,2,1]);
diff2 = totNCD-totAD;
diff2 = permute(diff2,[3,2,1]);
diff3_N2 = tot_N2_NCD-tot_N2_AD;
diff3_N2 = diff3_N2';

X = permute(totAD,[2 3 1]);
Y = permute(totND,[2 3 1]);
Z = permute(totNCD,[2 3 1]);

RM = [3 5 17 18 20 27 41]; % remove nations in EXIOBASE for which we have no NRD
X(:,RM,:) = [];
Y(:,RM,:) = [];
Z(:,RM,:) = [];
diff(RM,:,:) = [];
diff2(RM,:,:) = [];
diff3_N2(RM,:) =[];
tot_N2_AD(:,RM) = [];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 3. Plotting out results                                                      %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%load in shortened names for nations and income cats
income_cat = xlsread(Supplementary_Data_File,'Country_Classifications','F3:F46');
[~, ~, nations_short] = xlsread(Supplementary_Data_File,'Country_Classifications','K3:K46');

% remove nations in EXIOBASE for which we have no NRD
nations_short(RM) = [];
pop(RM) = [];
income_cat(RM) = [];

I = find(income_cat==2);
J = find(income_cat==3);
P = find(income_cat == 4);
Q = find(income_cat == 5);
K = [I' J' P' Q'];
X = X(:,K,:);
Y = Y(:,K,:);
diff = diff(K,:,:);
diff2 = diff2(K,:,:);
diff3_N2 = diff3_N2(K,:);
tot_N2_AD = tot_N2_AD(:,K);
pop = pop(K);
income_cat = income_cat(K);
nations_short = nations_short(K);

%% call text script for writing out results
writeOutResults(X,diff,diff2,income_cat,pop)

%% Relative difference plot - stacked bar
figure(1) % FIGURE 2B main text
y_lab = {'GHGs [kg CO2_{eq}]','Eutroph. [kg PO_{4}^{3-}]','Land [ha]'};
[h, sb, leg, h3, leg2] = funMultiStackedBar({diff2(:,:,1) diff2(:,:,2) diff2(:,:,3)},grps_str,nations_short,y_lab,[-0.05 0.03 -0.05 0]);
[h, sb, leg, h3, leg2] = funMSB_post_diet_impacts(h, sb, leg, h3, leg2, [-2.1 1.5; -50 20; -0.4 0.3]);
hold on
text(-2.5,-0.6,'B','fontangle','Italic','fontsize',22)
set(leg(2),'Visible','on')

export_fig([fig_str_LQ 'F2B.pdf']);
export_fig([fig_str_PR 'F2B.eps'],'-painters');

%% box plot of differences
figure(99)
h = funboxplot(diff(1:10,:,:),diff2(1:10,:,:),X(:,1:10,:));
xlabel('Low/middle income (relative difference in %)')
export_fig('../Text/figures/box_overview_low_mid.pdf');

figure(99)
h = funboxplot(diff(10:end,:,:),diff2(10:end,:,:),X(:,10:end,:));
xlabel('High income (relative difference in %)')
export_fig('../Text/figures/box_overview_high.pdf');

figure(2) % FIGURE 3 main text
h = funboxplot2(diff2(1:9,:,:),diff2(10:end,:,:),X(:,1:9,:),X(:,10:end,:),pop);
xlabel('Relative difference [%]')
export_fig([fig_str_LQ 'F3.pdf']);
export_fig([fig_str_PR 'F3.eps'],'-painters');

%% plot absolute differences
figure(3) % FIGURE S5
[h, sb, leg] = funMultiStackedBar({Y(:,:,1)' Y(:,:,2)' Y(:,:,3)'},grps_str,nations_short,y_lab,[-0.05 0.03 -0.05 0]);
[h, sb, leg, h3, leg2] = funMSB_post_diet_impacts(h, sb, leg, h3, leg2, [0 6; 0 130; 0 1]);
set(leg(2),'Visible','on');
export_fig([fig_str_LQ 'S5.pdf']);
export_fig([fig_str_PR 'S5.eps'],'-painters');

figure(4) % FIGURE 2A main text
[h, sb, leg] = funMultiStackedBar({X(:,:,1)' X(:,:,2)' X(:,:,3)'},grps_str,nations_short,y_lab,[-0.05 0.03 -0.05 0]);
[h, sb, leg, h3, leg2] = funMSB_post_diet_impacts(h, sb, leg, h3, leg2, [0 6; 0 130; 0 1]);
set(leg(2),'Visible','on');
hold on
text(-2.5,-0.29,'A','fontangle','Italic','fontsize',22)

export_fig([fig_str_LQ 'F2A.pdf']);
export_fig([fig_str_PR 'F2A.eps'],'-painters');

%% Plot map of environmental impacts
% FIGURES S6 (a), (b), and (c)
[~, ~, nations] = xlsread(Supplementary_Data_File,'Country_Classifications','C3:C46');
nations(RM) = [];
nations = nations(K);

figure(5) % S6(a)
ratio = (sum(diff2(:,:,1),2)./sum(X(:,:,1)',2)*100);
h = funMapping(ratio,nations,'Relative impact [%]',2,'Exiobase');
str_overlay('(a) GHG emissions',[0.001 0.03],22)
export_fig([fig_str_LQ 'S6a.pdf']);
export_fig([fig_str_PR 'S6a.eps'],'-painters');

figure(6) % S6(b)
ratio = (sum(diff2(:,:,2),2)./sum(X(:,:,2)',2)*100);
h = funMapping(ratio,nations,'Relative impact [%]',2,'Exiobase');
str_overlay('(b) Eutrophication',[0.001 0.03],22)
export_fig([fig_str_LQ 'S6b.pdf']);
export_fig([fig_str_PR 'S6b.eps'],'-painters');

figure(7) % S6(c)
ratio = (sum(diff2(:,:,3),2)./sum(X(:,:,3)',2)*100);
h = funMapping(ratio,nations,'Relative impact [%]',2,'Exiobase');
str_overlay('(c) Land use',[0.001 0.03],22)
export_fig([fig_str_LQ 'S6c.pdf']);
export_fig([fig_str_PR 'S6c.eps'],'-painters');

%% plot Eutro vs N2 regression
% FIGURE S10
figure(8)
diff_Eutro = squeeze(diff2(:,:,2));
diff_Eutro = sum(diff_Eutro,2);
diff3_N2 = sum(diff3_N2,2);
Eutro_x = sum(X(:,:,2),1)';
N2_x = sum(tot_N2_AD,1)';
x = (diff3_N2./N2_x)*100;
y = (diff_Eutro./Eutro_x)*100;
mdl = fitlm(x,y)
coeffs = polyfit(x,y,1);
fittedX = linspace(min(x), max(x),200);
fittedY = polyval(coeffs, fittedX);
plot(x,y,'.','markersize',17);
hold on
plot(fittedX,fittedY,'k-','Linewidth',1.5); %coeffs = twodp(coeffs);
text(-40,40,{'y = 0.89x + 3.3','R^2 = 0.841'},'fontsize',15)
xlim([-50 40])
xlabel('Eutro. [% reduction in NOx eq.]')
ylabel('Eutro. [% reduction in PO^{3-}_{4} eq.]')
grid on
set(gca,'xtick',-50:10:40)
set(gca,'fontsize',15)
export_fig([fig_str_LQ 'S10.pdf']);
export_fig([fig_str_PR 'S10.eps'],'-painters');

%% Plot import export with arrows
% FIGURES S7ab, S8ab, S9ab
[~, ~, nations] = xlsread(Supplementary_Data_File,'Country_Classifications','C3:C46');
nat_cat_map = xlsread(Supplementary_Data_File,'Country_Classifications','M3:M46');
nat_pop_map = xlsread(Supplementary_Data_File,'Country_Classifications','D3:D46');

[Domestic, Import, Export, allBalance] = Exio_trade(squeeze(sum(impactsAD(1,:,:,:),2)));

imp_dep = Import./(Import+Domestic+Export);

[Domestic, Import, Export, allBalance] = Exio_trade(squeeze(sum(impactsAD(1,:,:,:),2)),'Agg',[nat_cat_map nat_pop_map]);
imp_dep2 = Import./(Import+Domestic+Export);

for i = 1:1:max(nat_cat_map)
    imp_dep(nat_cat_map == i) = imp_dep2(i);
end
imp_dep = imp_dep.*100;

figure(9)
[h, cb] = funMapping(imp_dep,nations,'GHG Imports [% of total]',1);
nat_latlon = xlsread(Supplementary_Data_File,'Country_Classifications','Q3:R14');

[M,I] = sort(abs(allBalance(:)));
[I_row, I_col] = ind2sub(size(allBalance),I);
I_row = flipud(I_row);
I_col = flipud(I_col);

for i = 1:1:10;
    funAnnotateMap(nat_latlon(I_row(i),:),nat_latlon(I_col(i),:),allBalance(I_row(i),I_col(i))*1000,max(max(allBalance)));
end

setm(h,'flonlim',[-177 183]);
set(gca,'Position',[0 -0.2 1.15 1.3]);
set(cb,'location','south');
set(cb,'Position',[0.3 0.05 0.5 0.04]);
set(cb,'FontSize',20);
cb.Label.Rotation = 0;
cb.Label.Position = [(max(imp_dep)+min(imp_dep))/2 2.3 0];
str_overlay('(a)',[0.001 0.03],22)

export_fig([fig_str_LQ 'S7a.pdf']);
export_fig([fig_str_PR 'S7a.eps'],'-painters');
clf;

%figure 2 plot of changes due to NRDs
deltaTrade = impactsNCD;
[Domestic, Import, Export, allBalance2] = Exio_trade(squeeze(sum(deltaTrade(1,:,:,:),2)));

imp_dep = Import./(Import+Domestic+Export)-imp_dep;

[Domestic, Import, Export, allBalance2] = Exio_trade(squeeze(sum(deltaTrade(1,:,:,:),2)),'Agg',[nat_cat_map nat_pop_map]);
imp_dep2 = Import./(Import+Domestic+Export)-imp_dep2;

allBalance = allBalance2-allBalance;
for i = 1:1:max(nat_cat_map)
    imp_dep(nat_cat_map == i) = imp_dep2(i);
end
imp_dep = imp_dep.*100;

[h, cb] = funMapping(imp_dep,nations,'\Delta GHG Imports [% of total]',2);
nat_latlon = xlsread(Supplementary_Data_File,'Country_Classifications','Q3:R14');

[M,I] = sort(abs(allBalance(:)));
[I_row, I_col] = ind2sub(size(allBalance),I);
I_row = flipud(I_row);
I_col = flipud(I_col);

for i = 1:1:10;
    funAnnotateMap(nat_latlon(I_row(i),:),nat_latlon(I_col(i),:),allBalance(I_row(i),I_col(i))*1000,max(max(allBalance)));
end

setm(h,'flonlim',[-177 183]);
set(gca,'Position',[0 -0.2 1.15 1.3]);
set(cb,'location','south');
set(cb,'Position',[0.3 0.05 0.5 0.04]);
set(cb,'FontSize',20);
cb.Label.Rotation = 0;
cb.Label.Position = [(max(imp_dep)+min(imp_dep))/2 2.3 0];
str_overlay('(b)',[0.001 0.03],22)

export_fig([fig_str_LQ 'S7b.pdf']);
export_fig([fig_str_PR 'S7b.eps'],'-painters'); clf;

%% water
figure(10)
[Domestic, Import, Export, allBalance] = Exio_trade(squeeze(sum(impactsAD(2,:,:,:),2)));

imp_dep = Import./(Import+Domestic+Export);

[Domestic, Import, Export, allBalance] = Exio_trade(squeeze(sum(impactsAD(2,:,:,:),2)),'Agg',[nat_cat_map nat_pop_map]);
imp_dep2 = Import./(Import+Domestic+Export);

for i = 1:1:max(nat_cat_map)
    imp_dep(nat_cat_map == i) = imp_dep2(i);
end
imp_dep = imp_dep.*100;

[h, cb] = funMapping(imp_dep,nations,'Eutroph. Imports [% of total]',1);
nat_latlon = xlsread(Supplementary_Data_File,'Country_Classifications','Q3:R14');

[M,I] = sort(abs(allBalance(:)));
[I_row, I_col] = ind2sub(size(allBalance),I);
I_row = flipud(I_row);
I_col = flipud(I_col);

for i = 1:1:10;
    funAnnotateMap(nat_latlon(I_row(i),:),nat_latlon(I_col(i),:),allBalance(I_row(i),I_col(i))*1000,max(max(allBalance)));
end

setm(h,'flonlim',[-177 183]);
set(gca,'Position',[0 -0.2 1.15 1.3]);
set(cb,'location','south');
set(cb,'Position',[0.3 0.05 0.5 0.04]);
set(cb,'FontSize',20);
cb.Label.Rotation = 0;
cb.Label.Position = [(max(imp_dep)+min(imp_dep))/2 2.3 0];
str_overlay('(a)',[0.001 0.03],22)

export_fig([fig_str_LQ 'S8a.pdf']);
export_fig([fig_str_PR 'S8a.eps'],'-painters'); clf;

%figure 2 plot of changes due to NRDs
deltaTrade = impactsNCD;
[Domestic, Import, Export, allBalance2] = Exio_trade(squeeze(sum(deltaTrade(2,:,:,:),2)));

imp_dep = Import./(Import+Domestic+Export)-imp_dep;

[Domestic, Import, Export, allBalance2] = Exio_trade(squeeze(sum(deltaTrade(2,:,:,:),2)),'Agg',[nat_cat_map nat_pop_map]);
imp_dep2 = Import./(Import+Domestic+Export)-imp_dep2;

allBalance = allBalance2-allBalance;
for i = 1:1:max(nat_cat_map)
    imp_dep(nat_cat_map == i) = imp_dep2(i);
end
imp_dep = imp_dep.*100;

[h, cb] = funMapping(imp_dep,nations,'\Delta Eutroph. Imports [% of total]',2);
nat_latlon = xlsread(Supplementary_Data_File,'Country_Classifications','Q3:R14');

[M,I] = sort(abs(allBalance(:)));
[I_row, I_col] = ind2sub(size(allBalance),I);
I_row = flipud(I_row);
I_col = flipud(I_col);

for i = 1:1:10;
    funAnnotateMap(nat_latlon(I_row(i),:),nat_latlon(I_col(i),:),allBalance(I_row(i),I_col(i))*1000,max(max(allBalance)));
end

setm(h,'flonlim',[-177 183]);
set(gca,'Position',[0 -0.2 1.15 1.3]);
set(cb,'location','south');
set(cb,'Position',[0.3 0.05 0.5 0.04]);
set(cb,'FontSize',20); cb.Label.Rotation = 0; cb.Label.Position = [(max(imp_dep)+min(imp_dep))/2 2.3 0];
str_overlay('(b)',[0.001 0.03],22)

export_fig([fig_str_LQ 'S8b.pdf']);
export_fig([fig_str_PR 'S8b.eps'],'-painters'); clf;

%% Land
figure(11)
[Domestic, Import, Export, allBalance] = Exio_trade(squeeze(sum(impactsAD(3,:,:,:),2)));

imp_dep = Import./(Import+Domestic+Export);

[Domestic, Import, Export, allBalance] = Exio_trade(squeeze(sum(impactsAD(3,:,:,:),2)),'Agg',[nat_cat_map nat_pop_map]);
imp_dep2 = Import./(Import+Domestic+Export);

for i = 1:1:max(nat_cat_map)
    imp_dep(nat_cat_map == i) = imp_dep2(i);
end
imp_dep = imp_dep.*100;

[h, cb] = funMapping(imp_dep,nations,'Land Use Imports [% of total]',1);
nat_latlon = xlsread(Supplementary_Data_File,'Country_Classifications','Q3:R14');

[M,I] = sort(abs(allBalance(:)));
[I_row, I_col] = ind2sub(size(allBalance),I);
I_row = flipud(I_row); I_col = flipud(I_col);

for i = 1:1:10;
    funAnnotateMap(nat_latlon(I_row(i),:),nat_latlon(I_col(i),:),allBalance(I_row(i),I_col(i))*1000,max(max(allBalance)));
end

setm(h,'flonlim',[-177 183]);
set(gca,'Position',[0 -0.2 1.15 1.3]);
set(cb,'location','south');
set(cb,'Position',[0.3 0.05 0.5 0.04]);
set(cb,'FontSize',20); cb.Label.Rotation = 0;
cb.Label.Position = [(max(imp_dep)+min(imp_dep))/2 2.3 0];
str_overlay('(a)',[0.001 0.03],22)

export_fig([fig_str_LQ 'S9a.pdf']);
export_fig([fig_str_PR 'S9a.eps'],'-painters'); clf;

%figure 2 plot of changes due to NRDs
deltaTrade = impactsNCD;
[Domestic, Import, Export, allBalance2] = Exio_trade(squeeze(sum(deltaTrade(3,:,:,:),2)));

imp_dep = Import./(Import+Domestic+Export)-imp_dep;

[Domestic, Import, Export, allBalance2] = Exio_trade(squeeze(sum(deltaTrade(3,:,:,:),2)),'Agg',[nat_cat_map nat_pop_map]);
imp_dep2 = Import./(Import+Domestic+Export)-imp_dep2;

allBalance = allBalance2-allBalance;
for i = 1:1:max(nat_cat_map)
    imp_dep(nat_cat_map == i) = imp_dep2(i);
end

imp_dep = imp_dep.*100;

[h, cb] = funMapping(imp_dep,nations,'\Delta Land Use Imports [% of total]',2);
nat_latlon = xlsread(Supplementary_Data_File,'Country_Classifications','Q3:R14');

[M,I] = sort(abs(allBalance(:)));
[I_row, I_col] = ind2sub(size(allBalance),I);
I_row = flipud(I_row); I_col = flipud(I_col);

for i = 1:1:10;
    funAnnotateMap(nat_latlon(I_row(i),:),nat_latlon(I_col(i),:),allBalance(I_row(i),I_col(i))*1000,max(max(allBalance)));
end

setm(h,'flonlim',[-177 183]);
set(gca,'Position',[0 -0.2 1.15 1.3]);
set(cb,'location','south');
set(cb,'Position',[0.3 0.05 0.5 0.04]);
set(cb,'FontSize',20);
cb.Label.Rotation = 0;
cb.Label.Position = [(max(imp_dep)+min(imp_dep))/2 2.3 0];
str_overlay('(b)',[0.001 0.03],22)

export_fig([fig_str_LQ 'S9b.pdf']);
export_fig([fig_str_PR 'S9b.eps'],'-painters');