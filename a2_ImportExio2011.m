%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   a2_ImportExio2011.m                                                        %
%                                                                              %
%   Leiden University College &                                                %
%   Institute of environmental sciences (CML), Leiden University               %
%                                                                              %
%   Import EXIOBASE 3.3 for specific year                                      %
%                                                                              %
%   This script imports the EXIOBASE transaction input-output table for a      %
%   specific year and pre-computes some important matricies. Transactions are  %
%   given in basic price. The input-output table is of the form                %
%   product x product according industry technology assumptions.               %
%                                                                              %
%   Paul Behrens: p.a.behrens@luc.leidenuniv.nl                                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all; clc;

% Set Exiobase parameters
version = '3.3';
year_cnt = 2011;
year = num2str(year_cnt);

% Set paths
strpath = '~/Source/EXIOBASE/';
strpath = strcat(strpath, year, '/');

% Set variables for read in
strname = {
		'y';
		'v';
		'z';
		'r';		
		'h';
        'm';
        'n';
        'k';
	};

% Set filenames
strfile = {
	strcat('mrFinalDemand_', version, '_', year, '.txt');
    strcat('mrFactorInputs_', version, '_', year, '.txt');
    strcat('mrIot_', version, '_', year, '.txt');
	strcat('mrEmissions_', version, '_', year, '.txt');
	strcat('mrFDEmissions_', version, '_', year, '.txt');
    strcat('mrMaterials_', version, '_', year, '.txt');
	strcat('mrFDMaterials_', version, '_', year, '.txt');
    strcat('mrResources_', version, '_', year, '.txt');
	};

% Set skips for csvread
vskip = [
	2 3
	2 2
	2 3
	2 3
	2 3
    2 2
    2 2
    2 3
	];

% Read in data
for i = 1 : length(strname);
	tic
	disp(['Reading ',strname{i}]);
	tmp = dlmread([strpath, strfile{i}],'\t', vskip(i, 1), vskip(i, 2));
	toc
	disp(size(tmp));
	eval(['io.', strname{i},' = tmp;']);
	clear tmp;
end

% Compute relavant matrices
x = sum(io.z,2)+sum(io.y,2);      % total product use
io.x = x;
xinv = (1./(x+(x==0))).*(x~=0);   

io.A = io.z * diag(xinv);         % technology coefficients matrix
io.rc = io.r * diag(xinv);        % emission coefficients matrix
io.mc = io.m * diag(xinv);        % material coefficients matrix 
io.kc = io.k * diag(xinv);        % resource coefficients matrix

n = size(io.A,1);
io.L = inv(eye(n) - io.A);        % create Leontif inverse

% save out as matlab binary for space and speed
save('~/Source/EXIOBASE/2011/sys_pxp.mat','-v7.3','io');
clear all
