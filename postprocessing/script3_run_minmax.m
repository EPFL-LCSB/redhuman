% script to run thermo minmax
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

minmax_model2_smin = redHUMAN_minmax(model2_smin,thermoData_Recon2,data_Leukemia);
gem2 = model2_smin.OriginalGEM;
gem2.OriginalGEM = gem2;
minmax_gem2 = redHUMAN_minmax(gem2,thermoData_Recon2,minmax_model2_smin(:,1));
% save('./outputs/MinMax_Recon2.mat','minmax_model2_smin','minmax_gem2') 


% ... reduced models Recon 3
model3_smin = load('../models/redHUMAN_recon3_smin.mat');
model3_smin = model3_smin.redHUMAN_recon3_smin_29Aug2019_122558;
load('../data/thermoData_Recon3.mat');

minmax_model3_smin = redHUMAN_minmax(model3_smin,thermoData_Recon3,data_Leukemia);
gem3 = model3_smin.OriginalGEM;
gem3.OriginalGEM = gem3;
minmax_gem3 = redHUMAN_minmax(gem3,thermoData_Recon3,minmax_model3_smin(:,1));
% save('./outputs/MinMax_Recon3.mat','minmax_model3_smin','minmax_gem3') 

                                                 
%% plot rank Recon2
minmax_model2_smin(:,4) = num2cell(cell2mat(minmax_model2_smin(:,3)) - cell2mat(minmax_model2_smin(:,2)));
minmax_gem2(:,4) = num2cell(cell2mat(minmax_gem2(:,3)) - cell2mat(minmax_gem2(:,2)));

diff_flexibility = cell2mat(minmax_gem2(:,4)) - cell2mat(minmax_model2_smin(:,4));

[B,I] = sort(abs(diff_flexibility));
t2 = [minmax_model2_smin(I,1:3) minmax_gem2(I,2:3) num2cell(diff_flexibility(I))];
perc_reduced = (cell2mat(minmax_gem2(I,4))-cell2mat(minmax_model2_smin(I,4)))./max(cell2mat(minmax_gem2(I,4)),cell2mat(minmax_model2_smin(I,4)));
% length(find(abs(perc_reduced)<0.2))/383

reactions = minmax_model2_smin(I(end-19:end));
figure;
title = 'TVA redHUMAN Recon 2';
plotMinMaxredHUMAN2(reactions,cell2mat(minmax_gem2(I(end-19:end),2:3)),...
                   cell2mat(minmax_model2_smin(I(end-19:end),2:3)),title);

[~,ii] = ismember(reactions,model2_smin.rxns);
model2_smin.subSystems(ii)

%% plot rank Recon3
minmax_model3_smin(:,4) = num2cell(cell2mat(minmax_model3_smin(:,3)) - cell2mat(minmax_model3_smin(:,2)));
minmax_gem3(:,4) = num2cell(cell2mat(minmax_gem3(:,3)) - cell2mat(minmax_gem3(:,2)));

diff_flexibility = cell2mat(minmax_gem3(:,4)) - cell2mat(minmax_model3_smin(:,4));

[B,I] = sort(abs(diff_flexibility));
t3 = [minmax_model3_smin(I,1:3) minmax_gem3(I,2:3) num2cell(diff_flexibility(I))];
perc_reduced = (cell2mat(minmax_gem3(I,4))-cell2mat(minmax_model3_smin(I,4)))./max(cell2mat(minmax_gem3(I,4)),cell2mat(minmax_model3_smin(I,4)));


reactions = minmax_model3_smin(I(end-19:end));
figure;
title = 'TVA redHUMAN Recon 3';
plotMinMaxredHUMAN2(reactions,cell2mat(minmax_gem3(I(end-19:end),2:3)),...
                   cell2mat(minmax_model3_smin(I(end-19:end),2:3)),title);

[~,ii] = ismember(reactions,model3_smin.rxns);
model3_smin.subSystems(ii)



