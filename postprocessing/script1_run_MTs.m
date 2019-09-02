%script to test the metabolic tasks for...

%!!! LOAD TFA AND CPLEX PATHS
changeCobraSolver('cplex_direct', 'LP');

addpath(genpath('./utilities'))

% ... reduced models Recon 2
model2 = load('../models/redHUMAN_recon2.mat');
model2 = model2.redHUMAN_recon2_27Aug2019_093506;
model2_smin = load('../models/redHUMAN_recon2_smin.mat');
model2_smin = model2_smin.redHUMAN_recon2_smin_02Sep2019_135437;
load('../data/thermoData_recon2');

table_MTs_model2      = list_tasks_recon2(model2,thermoData_Recon2);
table_MTs_model2_smin = list_tasks_recon2(model2_smin,thermoData_Recon2);
% save('./outputs/MTs_Recon2.mat','table_MTs_model2','table_MTs_model2_smin')

% ... reduced models Recon 3
model3 = load('../models/redHUMAN_recon3.mat');
model3 = model3.redHUMAN_recon3_28Aug2019_121053;
model3_smin = load('../models/redHUMAN_recon3_smin.mat');
model3_smin = model3_smin.redHUMAN_recon3_smin_29Aug2019_122558;
load('../data/thermoData_recon3');

table_MTs_model3      = list_tasks_recon3(model3,thermoData_Recon3);
table_MTs_model3_smin = list_tasks_recon3(model3_smin,thermoData_Recon3);
% save('./outputs/MTs_Recon3.mat','table_MTs_model3','table_MTs_model3_smin')