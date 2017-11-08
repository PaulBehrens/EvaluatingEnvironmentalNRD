%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   b2_CompareNRDaverage.m                                                     % 
%                                                                              %
%   Leiden University College &                                                %
%   Institute of environmental sciences (CML), Leiden University               %
%                                                                              %
%   Compare the amounts in different diets for summary plots.                  %
%   Create plots for supplementary information.                                %
%                                                                              %
%                                                                              %
%   Paul Behrens: p.a.behrens@luc.leidenuniv.nl                                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all; clc;
addpath('./0. article_functions/')
fig_str_LQ = '../Figures_LQ/';
fig_str_PR = '../Figures_PR/';

%% First plot Raw NRD data results
load ../Data/unifiedNRD.mat

[cnt_idx] = unique([NRD{:,1}]);  [cat_idx, cat_str] = unique([NRD{:,8}]);
cat_str = NRD(cat_str,3);
[~,~,cnt_str] = xlsread('../Data/Supplementary_Data.xlsx','Country_Classifications','K3:K46');
[income_cat] = xlsread('../Data/Supplementary_Data.xlsx','Country_Classifications','F3:F46');
I = find(income_cat==2); J = find(income_cat==3); P = find(income_cat == 4); Q = find(income_cat ==5); K = [I' J' P' Q'];
RM = [3 14 26 27 30 35 43]; K(RM) =[]; % Remove nations for which we have no NRD

for i = 1:1:length(NRD)
    A(NRD{i,1},NRD{i,8},1) = NRD{i,4};
    A(NRD{i,1},NRD{i,8},2) = NRD{i,9};
    A(NRD{i,1},NRD{i,8},3) = NRD{i,10};
    A(NRD{i,1},NRD{i,8},4) = NRD{i,11};
end

% Use an aggregation to make the plot easier to read.
% FIGURE S2
sum_idx = [1 2 3 4 5 6 7 8 9 10 10 8 8 11 11 11 12];
cat_str{10} = 'Dairy'; cat_str{8} = 'Fats'; cat_str{11} = 'Red Meat'; cat_str{12} = 'Other meat'; cat_str = cat_str(1:12);

for i = 1:1:4
    PD{i} = A(cnt_idx,cat_idx,i);
    PDx{i} = funAccumarray(sum_idx,PD{i}');
end

y_lab = {'Amount [g]','Energy [kcal]','Protein [g]','Fats [g]'};
[h, sb, leg] = funMultiStackedBar({PDx{1}(:,K)' PDx{2}(:,K)' PDx{3}(:,K)' PDx{4}(:,K)'},cat_str,cnt_str(K),y_lab,[-0.05 0.03 -0.05 0]);
[h, sb, leg, h3, leg2] = funMSB_post_diet_impacts(h, sb, leg, 0, 0, [0 2600; 0 2500; 0 150; 0 100]); % call user function
set(leg(1),'Visible','on');
set(gcf,'Position',[94 51 920 800])
export_fig([fig_str_LQ 'S2.pdf'])
export_fig([fig_str_PR 'S2.eps'],'-painters')

%% Plot comparisons with the other diets
Type_code = xlsread('../Data/Supplementary_Data.xlsx','Diet_Classifications','L2:L89');
load ../Data/FAO.mat

FAO_less_waste_g = FAO_less_waste./repmat(FAO_less_waste(:,:,1),1,1,4);
FAO_less_waste_g = zeroinf(FAO_less_waste_g);
I = size(FAO_less_waste_g,3);

load ../Data/UnifiedNRD2.mat
NRD = Data;

kcal_g = squeeze(FAO_less_waste(:,:,2))./squeeze(FAO_less_waste(:,:,1));
kcal_g = zeroinf(kcal_g);
FAO_data = squeeze(FAO_less_waste(:,:,1));

NCD = zeros(size(NRD));
ratio = sum((kcal_g.*FAO_data),2)./sum((kcal_g.*Data),2);
tmp = (kcal_g.*Data).*repmat(ratio,1,88);
NCD = tmp./kcal_g; NCD = zeroinf(NCD);

for i = 1:1:I
   NRD_nut{i} = NRD.*FAO_less_waste_g(:,:,i);
   NCD_nut{i} = NCD.*FAO_less_waste_g(:,:,i);
   NRD_nut{i} = funAccumarray(Type_code,NRD_nut{i}');
   NCD_nut{i} = funAccumarray(Type_code,NCD_nut{i}');
   FAO_nut{i} = funAccumarray(Type_code,squeeze(FAO_less_waste(:,:,i))');
end

NRD_nut{3} = ((NRD_nut{3}(1:6,:).*4)./repmat(sum(NRD_nut{2}),6,1)).*100; % 4 is kcal/gram protien
NRD_nut{4} = ((NRD_nut{4}(1:6,:).*9)./repmat(sum(NRD_nut{2}),6,1)).*100; % 9 is kcal/gram sugars

NCD_nut{3} = ((NCD_nut{3}(1:6,:).*4)./repmat(sum(NCD_nut{2}),6,1)).*100;
NCD_nut{4} = ((NCD_nut{4}(1:6,:).*9)./repmat(sum(NCD_nut{2}),6,1)).*100;

FAO_nut{3} = ((FAO_nut{3}(1:6,:).*4)./repmat(sum(FAO_nut{2}),6,1)).*100;
FAO_nut{4} = ((FAO_nut{4}(1:6,:).*9)./repmat(sum(FAO_nut{2}),6,1)).*100;


%% plot out summary information figures
figure(1) % FIGURE S1
l_lab = {'Meat','Fish','Dairy','Grains','VFN','Other'};
y_lab = {'Energy [kcal]','Protein [% of energy]','Fats [% of energy]'};
[h, sb, leg] = funMultiStackedBar({FAO_nut{2}(:,K)' FAO_nut{3}(:,K)' FAO_nut{4}(:,K)'},l_lab,cnt_str(K),y_lab,[-0.05 0.03 -0.05 0]);
[h, sb, leg, h3, leg2] = funMSB_post_diet_impacts(h, sb, leg, 0, 0, [0 3200; 0 17; 0 60]);
set(leg(2),'Visible','on');
set(gcf,'Position',[94 51 920 800])
export_fig([fig_str_LQ 'S1.pdf'])
export_fig([fig_str_PR 'S1.eps'],'-painters')

figure(2)
l_lab = {'Meat','Fish','Dairy','Grains','VFN','Other'};
y_lab = {'Energy [kcal]','Protein [% of energy]','Fats [% of energy]'};
[~, C] = sort(sum(NRD_nut{2},1),'descend');
[h, sb, leg] = funMultiStackedBar({NRD_nut{2}(:,K)' NRD_nut{3}(:,K)' NRD_nut{4}(:,K)'},l_lab,cnt_str(K),y_lab,[-0.05 0.03 -0.05 0]);
[h, sb, leg, h3, leg2] = funMSB_post_diet_impacts(h, sb, leg, 0, 0, [0 3700; 0 19; 0 60]);
set(leg(2),'Visible','on');
set(gcf,'Position',[94 51 920 754])
%export_fig('../Text/figures/NRD.pdf')

figure(3) % FIGURE S4
l_lab = {'Meat','Fish','Dairy','Grains','VFN','Other'};
y_lab = {'Energy [kcal]','Protein [% change in energy]','Fats [% change in energy]'};
[h, sb, leg] = funMultiStackedBar({NRD_nut{2}(:,K)'-FAO_nut{2}(:,K)' NRD_nut{3}(:,K)'-FAO_nut{3}(:,K)' NRD_nut{4}(:,K)'-FAO_nut{4}(:,K)'},l_lab,cnt_str(K),y_lab,[-0.05 0.03 -0.05 0]);
[h, sb, leg, h3, leg2] = funMSB_post_diet_impacts(h, sb, leg, 0, 0, [-1500 1600; -10 10; -21 15]);
set(leg(2),'Visible','on');
A = get(sb(1),'ylabel');
A.Position = [-0.05 0.5 0];
A = get(sb(2),'ylabel');
A.String = {'Protein';'[% change in energy]'};
set(sb(2),'ylabel',A);
A = get(sb(3),'ylabel');
A.String = {'Fats';'[% change in energy]'};
set(sb(3),'ylabel',A);

set(gcf,'Position',[94 51 920 800])
export_fig([fig_str_LQ 'S4.pdf'])
export_fig([fig_str_PR 'S4.eps'],'-painters')

figure(4) % FIGURE S3
l_lab = {'Meat','Fish','Dairy','Grains','VFN','Other'};
y_lab = {'Energy [kcal]','Protein [% change in energy]','Fats [% change in energy]'};
[h, sb, leg] = funMultiStackedBar({NCD_nut{2}(:,K)'-FAO_nut{2}(:,K)' NCD_nut{3}(:,K)'-FAO_nut{3}(:,K)' NCD_nut{4}(:,K)'-FAO_nut{4}(:,K)'},l_lab,cnt_str(K),y_lab,[-0.05 0.03 -0.05 0]);
[h, sb, leg, h3, leg2] = funMSB_post_diet_impacts(h, sb, leg, 0, 0, [-1500 1600; -10 10; -21 15]);
set(leg(2),'Visible','on');
A = get(sb(1),'ylabel');
A.Position = [-0.05 0.5 0];
A = get(sb(2),'ylabel');
A.String = {'Protein';'[% change in energy]'};
set(sb(2),'ylabel',A);
A = get(sb(3),'ylabel');
A.String = {'Fats';'[% change in energy]'};
set(sb(3),'ylabel',A);

set(gcf,'Position',[94 51 920 800])
export_fig([fig_str_LQ 'S3.pdf'])
export_fig([fig_str_PR 'S3.eps'],'-painters')

%% NOTE: figure 1 main text
figure(5)  
l_lab = {'Meat','Fish','Dairy','Grains','VFN','Other'};
y_lab = {'Average diet [kcal]','NRD [kcal]','Difference [kcal]'};
[h, sb, leg] = funMultiStackedBar({FAO_nut{2}(:,K)' NCD_nut{2}(:,K)' NCD_nut{2}(:,K)'-FAO_nut{2}(:,K)'},l_lab,cnt_str(K),y_lab,[-0.05 0.03 -0.05 0]);
[h, sb, leg, h3, leg2] = funMSB_post_diet_impacts(h, sb, leg, 0, 0, [0 3200; 0 3200 ; -1000 1000]);
set(leg(2),'Visible','on');

for i = 1:3
    A = get(sb(i),'ylabel');
    A.Position = [-0.05 0.5 0];
end

set(gcf,'Position',[94 51 920 800])
export_fig([fig_str_LQ 'F1.pdf'])
export_fig([fig_str_PR 'F1.eps'],'-painters')
