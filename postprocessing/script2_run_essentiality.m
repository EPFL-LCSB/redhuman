%script to run gene essentiality
clear;clc
% LOAD TFA PATH
% LOAD CPLEX
changeCobraSolver('cplex_direct', 'LP');
addpath(genpath('./utilities'))

% ... reduced models Recon 2
model2_smin = load('../models/redHUMAN_recon2_smin.mat');
model2_smin = model2_smin.redHUMAN_recon2_smin_02Sep2019_135437;
load('../data/thermoData_Recon2.mat');
load('../data/data_Leukemia.mat')

[EssGenes_leukemia2, model2_smin_leukemia] = redHUMAN_geneEss(model2_smin,thermoData_Recon2,data_Leukemia);
gem2 = model2_smin.OriginalGEM;
gem2.OriginalGEM = gem2;
[EssGenes_GEM_leukemia2, gem2_leukemia] = redHUMAN_geneEss(gem2,thermoData_Recon2,data_Leukemia);
% save('./outputs/EssGenes_Recon2.mat','EssGenes_leukemia2','EssGenes_GEM_leukemia2') 

% ... reduced models Recon 3
model3_smin = load('../models/redHUMAN_recon3_smin.mat');
model3_smin = model3_smin.redHUMAN_recon3_smin_29Aug2019_122558;
load('../data/thermoData_Recon3.mat');

[EssGenes_leukemia3, model3_smin_leukemia] = redHUMAN_geneEss(model3_smin,thermoData_Recon3,data_Leukemia);
gem3 = model3_smin.OriginalGEM;
gem3.OriginalGEM = gem3;
[EssGenes_GEM_leukemia3, gem3_leukemia] = redHUMAN_geneEss(gem3,thermoData_Recon3,data_Leukemia);
% save('./outputs/EssGenes_Recon3.mat','EssGenes_leukemia3','EssGenes_GEM_leukemia3') 

