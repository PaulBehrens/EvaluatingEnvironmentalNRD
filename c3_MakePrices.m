%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   c3_MakePrices.m                                                            %
%                                                                              %
%   Leiden University College &                                                %
%   Institute of environmental sciences (CML), Leiden University               %
%                                                                              %
%   Calculate prices for EXIOBASE products using FAO statistics and EXIO       %  
%   monetary data. Final prices are in Euros/g                                 %
%                                                                              %
%   Paul Behrens: p.a.behrens@luc.leidenuniv.nl                                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all; clc;

Supplementary_Data_File = '../Data/Supplementary_Data.xlsx';
EXIO_code = xlsread(Supplementary_Data_File,'Diet_Classifications','J2:J89');
FAO_data = xlsread(Supplementary_Data_File,'g_per_cap_per_day','B3:DS46');
population = xlsread(Supplementary_Data_File,'Country_Classifications','D3:C46');

Data = FAO_data;

M = length(unique(EXIO_code));
N = size(Data,1);

for i = 1:1:N
    D(:,i) = accumarray(EXIO_code,Data(i,:))*population(i)*365;
end

M = size(D,1); N = size(D,2);
D(M+1:200,:) = zeros((200-M),44);
D = reshape(D,8800,1);
D = [D' zeros(1000,1)']';

FD = dlmread('~/Source/EXIOBASE/2011/mrFinalDemand_3.3_2011.txt','\t',2, 3);
FD = reshape(FD,200,49,[]);
FD = squeeze(sum(FD,2));
FD = reshape(FD,200,7,[]);
FD = squeeze(FD(:,1,:));    
FD = reshape(FD,[],1);
 
Prices = (FD./D);
Prices(isnan(Prices) == 1 | Prices == Inf) = 0;

Prices = reshape(Prices,200,49);
FD = reshape(FD,200,49);
D = reshape(D,200,49);

Prices = Prices.*10^6;  
save ../Data/Prices.mat Prices

