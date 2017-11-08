%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Script: a3_ProcessNRD.m                                                    %
%                                                                              %
%   Leiden University College &                                                %
%   Institute of environmental sciences (CML), Leiden University               %
%                                                                              %
%   Nationally recommended diets come in all sorts of units and                %
%   choices, this script unifies the diets for units, and choices. Script      %
%   writes out.                                                                %
%                                                                              %
%   Further steps:                                                             %
%   1. Splits vegetable/fruit and meat/fish advice along FAO proportions of    %
%      nationally consumed amounts.                                            %
%   2. Splits red meat, lean meat, and meat categories proportions in FAO for  %
%   meat categories                                                            %
%                                                                              %
%   Paul Behrens: p.a.behrens@luc.leidenuniv.nl                                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all; clc; N = [];


%% load in data

% Load NRD data
%[~, ~, N] = xlsread('../Data/Supplementary_Data.xlsx','Diet_Vector','A2:H386');
[~, ~, N] = xlsread('../Data/Supplementary_Data.xlsx','Diet_Vector2','A2:H387');

% Load NRD codes
[~, ~, NC] = xlsread('../Data/Supplementary_Data.xlsx','Diet_Classifications','R3:S23');

% Load default portions when grams are not available
[~, ~, E] = xlsread('../Data/Supplementary_Data.xlsx','Default_portions','A2:F22');

% Load fruit/veg and meat/fish splits
[~, ~, S] = xlsread('../Data/Supplementary_Data.xlsx','g_per_capita_per_day_less_waste','B5:CL48');

% Load FAO codes
[FAO_codes] = xlsread('../Data/Supplementary_Data.xlsx','Diet_Classifications','H2:H89');

cnt_str = unique([N{:,1}]);
cnt_length = length(cnt_str);

%% Catagorise fresh and aged cheese as just cheese!
I = find(strcmp(N(:,3),'Aged Cheese') == 1 | strcmp(N(:,3),'Fresh Cheese') == 1);
for i = 1:length(I)
   N{I(i),3} = 'Cheese';
end

%% Remove Fiber and Salt (not of interest in this work)
I = find([N{:,8}] == 190 | [N{:,8}] == 191);
N(I,:) = [];

%% Make Lean Meat, Non-Specified Meat
I = find(strcmp(N(:,3),'Lean Meat') == 1);
for i = 1:length(I)
   N{I(i),3} = 'Non-specified Meat';
   N{I(i),8} = 106; 
end

%% catagorise wholegrains as grains
I = find(strcmp(N(:,3),'Wholegrains') == 1);
for i = 1:length(I)
   N{I(i),3} = 'Grains';
   N{I(i),8} = 101; 
end

%% check codes
for i = 1:length(N)
    I = find(strcmp(NC(:,1),N{i,3}) == 1);
    if NC{I,2} - N{i,8} == 0
        T(i) = 0;
    else
        T(i) = 1;
    end
end

if sum(find(T==1)) > 0
    disp('ERROR')
end

%% Split 'choices' within nationally recommended diets evenly
for i = 1:cnt_length
    I = find([N{:,1}] == cnt_str(i));
    I2 = find(isnan([N{I,6}]) == 0);
    if isempty(I2)==0
        j_max = max([N{I,6}]);
        for j = 1:j_max
            I3 = find([N{I,6}] == j);
            for k = 1:length(I3)
                N{I(I3(k)),4} = [N{I(I3(k)),4}]./length(I3);
            end
        end
    end
end

%% Split fruit/veg & combined catagories
I = find(strcmp(N(:,3),'Fruit/Veg') == 1);
for i = 1:length(I)
    N(I(i)+1:end+1,:) = N(I(i):end,:);
    N{I(i),3} = 'Fruit'; N{I(i),8} = 104;
    N{I(i)+1,3} = 'Vegetables'; N{I(i)+1,8} = 105;
    J = find(strcmp(S(:,1),N(I(i),2)) == 1);
    if isempty(J) == 0
       N{I(i),4} =  N{I(i),4}.*(sum([S{J,45:55}])./(sum([S{J,42:44}])+sum([S{J,45:55}])));
       N{I(i)+1,4} = N{I(i)+1,4}-N{I(i),4};
    end
    I = I+1;
end

%% Split meat/fish catagories
I = find(strcmp(N(:,3),'Meat/Fish') == 1);
for i = 1:length(I)
    N(I(i)+1:end+1,:) = N(I(i):end,:);
    N{I(i),3} = 'Non-specified Meat'; N{I(i),8} = 106;
    N{I(i)+1,3} = 'Fish'; N{I(i)+1,8} = 109;
    J = find(strcmp(S(:,1),N(I(i),2)) ==1);
    if isempty(J) == 0
       N{I(i),4} =  N{I(i),4}.*(sum([S{J,67:71}])./(sum([S{J,78:89}])+sum([S{J,67:71}])));
       N{I(i)+1,4} = N{I(i)+1,4}-N{I(i),4};
    end
    I = I+1;
end

%% Convert servings and portions to grams numbers, portions, to grams
servings2gram = xlsread('../Data/Supplementary_Data.xlsx','Default_portions','B2:C22');
I = find(strcmp(N(:,5),'servings')==1 | strcmp(N(:,5),'portions')==1);
for i = 1:length(I)
    J = find(N{I(i),8} == servings2gram(:,1));
    N{I(i),4} = N{I(i),4}.*servings2gram(J,2);
    N{I(i),5} = 'grams';
end

%% Convert remaining units to grams
egg2gram = 50;  %(US medium Egg)
tablespoon2gram = 13.7; 
UScupveg2gram = 150;
UScupfruit2gram = 175;
UScupgrain2gram = 129.60;
cupmilk2gram = 236.588;
handful2gram = 28.3;
wholegrains2grains = 3.7;
teaspoon2gram = 15;

% change eggs to grams
I = find(strcmp(N(:,5),'number')==1 | strcmp(N(:,5),'units')==1 & strcmp(N(:,3),'Eggs')==1);
for i = 1:length(I)
   N{I(i),4} = N{I(i),4}.*egg2gram;
   N{I(i),5} = 'grams';
end

I = find(strcmp(N(:,5),'tablespoons')==1);
for i = 1:length(I)
    N{I(i),4} = N{I(i),4}.*tablespoon2gram;
    N{I(i),5} = 'grams';
end

I = find(strcmp(N(:,5),'teaspoon')==1);
for i = 1:length(I)
    N{I(i),4} = N{I(i),4}.*teaspoon2gram;
    N{I(i),5} = 'grams';
end

I = find(strcmp(N(:,5),'cups')==1);
for i = 1:length(I)
    if strcmp(N{I(i),3},'Grains') ==  1
        N{I(i),4} = N{I(i),4}.*UScupgrain2gram;
    elseif strcmp(N{I(i),3},'Vegetables') == 1
        N{I(i),4} = N{I(i),4}.*UScupveg2gram;
    elseif strcmp(N{I(i),3},'Fruit') == 1
        N{I(i),4} = N{I(i),4}.*UScupfruit2gram;
    elseif strcmp(N{I(i),3},'Milk') == 1
        N{I(i),4} = N{I(i),4}.*cupmilk2gram;
    end
    N{I(i),5} = 'grams';
end

I = find(strcmp(N(:,5),'handful')==1);
for i = 1:length(I)
    N{I(i),4} = N{I(i),4}.*handful2gram;
    N{I(i),5} = 'grams';
end

I = find(strcmp(N(:,5),'ml')==1);
for i = 1:length(I)
    N{I(i),5} = 'grams';
end

%% Combine US non-specified meat
I = find(strcmp(N(:,3),'Non-specified Meat') ==1 & strcmp(N(:,2),'United States') == 1);
N{I(1),4} = N{I(1),4} + N{I(2),4};
N(I(2),:) = [];

%% Split non-specified and red meat into broader catagories
meat_codes = [118 119 120 121];
[~,meat_FAO_codes] = intersect(FAO_codes,[2731 2732 2733 2734]);
meat_FAO_codes = meat_FAO_codes+1;
C = {'Bovine Meat','Mutton & Goat Meat','Pigmeat','Poultry Meat'};
for i = 1:cnt_length
        I = find([N{:,1}] == cnt_str(i));
        I2 = find(strcmp(N(I,3),'Non-specified Meat') == 1);
        if isempty(I2) == 0
            total = N{I(I2),4};
            for j = 1:length(meat_codes)
                N(I(I2)+1:end+1,:) = N(I(I2):end,:);
                N{I(I2),3} = C{j}; N{I(I2),8} = meat_codes(j);
                N{I(I2),4} = total.*(S{cnt_str(i),meat_FAO_codes(j)}/sum([S{cnt_str(i),meat_FAO_codes}]));
            end
        end
        
        I2 = find(strcmp(N(I,3),'Red Meat') == 1);
        if isempty(I2) == 0
            total = N{I(I2),4};
            for j = 1:length(meat_codes)-1
                N(I(I2)+1:end+1,:) = N(I(I2):end,:);
                N{I(I2),3} = C{j}; N{I(I2),8} = meat_codes(j);
                N{I(I2),4} = total.*(S{cnt_str(i),meat_FAO_codes(j)}/sum([S{cnt_str(i),meat_FAO_codes(1:3)}]));
            end
        end              
end

for i = 1:cnt_length
   I = find([N{:,1}] == cnt_str(i));
   I2 = find(strcmp(N(I,3),'Pigmeat') == 1);   
   if length(I2)>1
      for k = 1:length(meat_codes)-1
          I2 = find(strcmp(N(I,3),C{k})==1);
          N{I(I2(1)),4} = sum([N{I(I2),4}]);
          todelete(k) = I2(2);
      end
      N(I(todelete),:) = [];
    end
end

I = find(strcmp(N(:,3),'Red Meat') == 1 | strcmp(N(:,3),'Non-specified Meat') == 1);
N(I,:) = [];

%% Calculate energy, fats, and protiens for each food type
for i = 1:cnt_length
        I = find([N{:,1}] == cnt_str(i));
        for j = 1:length(I)
            I2 = find(strcmp(N(I(j),3),E(:,1))==1);
            if isempty(I2) == 0
                N{I(j),9} =  E{I2,4}.*(N{I(j),4}./100);
                N{I(j),10} =  E{I2,5}.*(N{I(j),4}./100);
                N{I(j),11} =  E{I2,6}.*(N{I(j),4}./100);
            end
        end
end

NRD = N;
%% Write out unified dataset
save('../Data/unifiedNRD.mat','NRD')