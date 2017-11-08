%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   c2_MakeDietMatrix.m                                                        %
%                                                                              %
%   Leiden University College &                                                %
%   Institute of environmental sciences (CML), Leiden University               %
%                                                                              %
%   Create vectors of diets, specified as fractions of total household         %
%   income spent on different food categories for each country/region.         %
%                                                                              %
%   Paul Behrens: p.a.behrens@luc.leidenuniv.nl                                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc; clear all;

load ../Data/AverageDiet.mat AD
load ../Data/NRD_Diet.mat ND
load ../Data/NRD_isocaloric.mat NCD

FD = dlmread('~/Source/EXIOBASE/2011/mrFinalDemand_3.3_2011.txt','\t',2, 3);

FD = reshape(FD,200,49,7,[]);
FD = squeeze(FD(:,:,1,:)); % select only HH consumption
totals = sum(FD,2); % create totals of consumption for each type

totals = repmat(totals,[1 49 1]); % create a division matrix for the FD matrix

fractions = FD./totals;
fractions(isnan(fractions) == 1) = 0;
fractions(isinf(fractions) == 1) = 0;

% Diet specified as fractions of total household income spent on different
% food categories.
for i = 1:1:49
    for j = 1:1:200
      tmp1(j,:,i) = fractions(j,:,i).*AD(j,i);
      tmp2(j,:,i) = fractions(j,:,i).*ND(j,i);
      tmp4(j,:,i) = fractions(j,:,i).*NCD(j,i);
    end
end

AD2 = reshape(tmp1,9800,49);  % Average Diet
ND2 = reshape(tmp2,9800,49);  % National Recommended Diet
NCD2 = reshape(tmp4,9800,49); % Iso Caloric Diet

AD = AD2; ND = ND2; NCD = NCD2;
save ../Data/Diets.mat AD ND NCD
