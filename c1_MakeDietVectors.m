%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Create Diets for stimulus vector                                           %
%                                                                              %
%   Leiden University College &                                                %
%   Institute of environmental sciences (CML), Leiden University               %
%                                                                              %
%   1. Average diet                                                            %
%   2. Nationally recommended/WHO modified diet                                %
%   3. Isocaloric diet                                                         %
%                                                                              %
%   Unit is g/year/person                                                      % 
%                                                                              %
%   Paul Behrens: p.a.behrens@luc.leidenuniv.nl                                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all; clc;

Supplementary_Data_File = '../Data/Supplementary_Data.xlsx';
EXIO_code = xlsread(Supplementary_Data_File,'Diet_Classifications','J2:J89');
NRD_code = xlsread(Supplementary_Data_File,'Diet_Classifications','K2:K89');
FAO_data = xlsread(Supplementary_Data_File,'g_per_capita_per_day_less_waste','B5:CL48');
kcal_g = xlsread(Supplementary_Data_File,'kcal per g','B3:CK46');
empty_cals = xlsread(Supplementary_Data_File,'Diet_Classifications','N2:N89');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1. Make average diet                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Data = FAO_data;

M = length(unique(EXIO_code));
N = size(Data,1);

for i = 1:1:N
    D(:,i) = accumarray(EXIO_code,Data(i,:));
end

M = size(D,1); N = size(D,2);
D(M+1:200,:) = zeros((200-M),44);
D = reshape(D,8800,1);
D = [D' zeros(1,200*5)];
D = reshape(D,200,49);

AD = D;
save ../Data/AverageDiet.mat AD
clear D

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2. Make nationally recommended diet                                          %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load ../Data/unifiedNRD.mat

for i = 1:length(NRD)
    I = find(NRD{i,8}==NRD_code);
    Data(NRD{i,1},I) = NRD{i,4}.*Data(NRD{i,1},I)./sum(Data(NRD{i,1},I));
end

% zero out other types of oil, and meat (by implication these are
% recommended to be oilve oil only, and red/white meats (not offal).

K = find(NRD_code==98);
for i = 1:1:44
    I = find([NRD{:,1}] == i);
    J = find([NRD{I,8}] == 111);
    if length(J)>0
        Data(i,K) = 0;
    end
end

% lower the empty calories to 350 a day max

kcal = Data.*kcal_g;
for i = 1:1:size(Data,1)
    empty(:,i) = sum(kcal(i,empty_cals==1));
    if empty(1,i) > 350
       empty_ratio = 350/empty(1,i);
       Data(i,empty_cals==1) = Data(i,empty_cals==1).*empty_ratio;
    end
end

save ../Data/unifiedNRD2.mat Data

% format the vector for processing
for i = 1:1:size(Data,1);
    D(:,i) = accumarray(EXIO_code,Data(i,:));
end

M = size(D,1); N = size(D,2);
D(M+1:200,:) = zeros((200-M),44);
D = reshape(D,8800,1);
D = [D' zeros(1,200*5)];
D = reshape(D,200,49);
ND = D;

save ../Data/NRD_Diet.mat ND

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 3. Make isocaloric diet
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Type_code = xlsread(Supplementary_Data_File,'Diet_Classifications','L2:L89');
load ../Data/FAO.mat

NCD = zeros(size(NRD));
ratio = sum((kcal_g.*FAO_data),2)./sum((kcal_g.*Data),2);
tmp = (kcal_g.*Data).*repmat(ratio,1,88);
NCD = tmp./kcal_g; NCD = zeroinf(NCD);

clear D
for i = 1:1:size(Data,1);
    D(:,i) = accumarray(EXIO_code,NCD(i,:));
end

M = size(D,1); N = size(D,2);
D(M+1:200,:) = zeros((200-M),44);
D = reshape(D,8800,1);
D = [D' zeros(1,200*5)];
D = reshape(D,200,49);
NCD = D;

save ../Data/NRD_isocaloric.mat NCD

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % 4. Make non-isocaloric diet (scaled to 2200 kcal/person/day)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% load ../Data/unifiedNRD2.mat 
% 
% NCD = zeros(size(NRD));
% ratio = 2400./sum((kcal_g.*Data),2);
% tmp = (kcal_g.*Data).*repmat(ratio,1,88);
% NCD = tmp./kcal_g; NCD = zeroinf(NCD);
% Data = NCD;
% save ../Data/unifiedNRD2.mat Data
% 
% clear D
% for i = 1:1:size(Data,1);
%     D(:,i) = accumarray(EXIO_code,NCD(i,:));
% end
% 
% M = size(D,1); N = size(D,2);
% D(M+1:200,:) = zeros((200-M),44);
% D = reshape(D,8800,1);
% D = [D' zeros(1,200*5)];
% D = reshape(D,200,49);
% ND = D;
% 
% save ../Data/NRD_Diet.mat ND
