%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   a1_processFAO.m                                                            %
%                                                                              %
%   Leiden University College &                                                %
%   Institute of environmental sciences (CML), Leiden University               %
%                                                                              %    
%   Extract FAO data from raw balance sheets                                   %
%                                                                              %
%   Apply waste percentages for food groups and nations                        %
%                                                                              %
%   Paul Behrens: p.a.behrens@luc.leidenuniv.nl                                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all; clc;

% FAO country codes for EXIOBASE match, 
% see "../Data/Supplementary_Data.xlsx" for explanation of the codes
cnt_codes = [11 255 27 50 167 79 54 ...
             63 203 67 68 84 98 97 104 106 ...
             126 256 119 134 150 173 ...
             174 183 210 198 199 229 231 110 351 ... 
             33 117 21 100 138 185 10 211 223 214 ...
             162 101 202];

[num, str, raw] = xlsread('../Data/FAO_raw.xlsx','Sheet 1','C2:L52609'); % Import raw FAO data
data = [num(:,1) num(:,3) num(:,5) num(:,end)];
[food_codes] = xlsread('../Data/Supplementary_Data.xlsx','Diet_Classifications','H2:H89'); % Import food codes required

% select category classification (646 is grams) (664 is calories) (674 is protein) (684 fats)
% see "../Data/Supplementary_Data.xlsx" for explanation of the codes
data_codes = [646 664 674 684];

% extract data from sheet
for k = 1:1:length(data_codes)
    for i = 1:1:length(cnt_codes)
        for j = 1:1:length(food_codes)    
             I = find(data(:,1)==cnt_codes(i) & data(:,2) == data_codes(k) & data(:,3) == food_codes(j));
             if isempty(I)==1
                 FAO_raw(i,j,k) = 0;
                else
                 FAO_raw(i,j,k) = data(I,4);
             end
            clear I
        end
    end
end

FAO_raw = zeroinf(FAO_raw); % zero NaNs (as empty cells appear from excel in NaN format)

%% Process wastes on result
[W_table, ~, ~] = xlsread('../Data/Supplementary_Data.xlsx','Waste_Percentages','L32:S38'); % load in concordances
[FAO_food,~, ~ ] = xlsread('../Data/Supplementary_Data.xlsx','Diet_Classifications','M2:M89'); % grams per day
[FAO_cnt, ~, ~] = xlsread('../Data/Supplementary_Data.xlsx','Country_Classifications','E3:E46'); % load in countries

[I,J,K] = size(FAO_raw);

for k = 1:1:K
    for i = 1:1:I
        for j = 1:1:J
         FAO_less_waste(i,j,k) = FAO_raw(i,j,k) * (1-(W_table(FAO_cnt(i),FAO_food(j)) /100));
        end
    end
end

FAO_less_waste = zeroinf(FAO_less_waste); % zero out divide by zeros

% save out
save ../Data/FAO.mat FAO_raw FAO_less_waste