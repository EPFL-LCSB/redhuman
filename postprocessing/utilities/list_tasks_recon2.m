function table_MTs = list_tasks_recon2(model,thermo_data)


table_MTs = {};

%% Rephosphorylation of nucleoside triphosphates
% Task 1: Aerobic rephosphorylation of ATP from glucose
task_name = 'Aerobic rephosphorylation of ATP from glucose';
task_should_fail = 0;
task_mets_in = {'o2_e';'glc_e'};
task_mets_in_lb = [NaN;NaN];
task_mets_in_ub = [NaN;NaN];
task_mets_out = {'h2o_e';'co2_e'};
task_mets_out_lb = [NaN;NaN];
task_mets_out_ub = [NaN;NaN];
task_rxns = {'atp_c + h2o_c <=> adp_c + pi_c + h_c'};
task_rxns_lb = 1;
task1_success = test_MT(model, task_name, task_should_fail,...
                task_mets_in, task_mets_in_lb, task_mets_in_ub,... 
                task_mets_out, task_mets_out_lb,task_mets_out_ub,...
                task_rxns, task_rxns_lb, thermo_data);
table_MTs = [table_MTs;task_name num2cell(task1_success)];

%Task2: Aerobic rephosphorylation of GTP
task_name = 'Aerobic rephosphorylation of GTP from glucose';
task_should_fail = 0;
task_mets_in = {'o2_e';'glc_e'};
task_mets_in_lb = [NaN;NaN];
task_mets_in_ub = [NaN;NaN];
task_mets_out = {'h2o_e';'co2_e'};
task_mets_out_lb = [NaN;NaN];
task_mets_out_ub = [NaN;NaN];
task_rxns = {'gtp_c + h2o_c <=> gdp_c + pi_c + h_c'};
task_rxns_lb = 1;
task2_success = test_MT(model, task_name, task_should_fail,...
                task_mets_in, task_mets_in_lb, task_mets_in_ub,... 
                task_mets_out, task_mets_out_lb,task_mets_out_ub,...
                task_rxns, task_rxns_lb, thermo_data);
table_MTs = [table_MTs;task_name num2cell(task2_success)];    

%Task3: Aerobic rephosphorylation of CTP
task_name = 'Aerobic rephosphorylation of CTP from glucose';
task_should_fail = 0;
task_mets_in = {'o2_e';'glc_e'};
task_mets_in_lb = [NaN;NaN];
task_mets_in_ub = [NaN;NaN];
task_mets_out = {'h2o_e';'co2_e'};
task_mets_out_lb = [NaN;NaN];
task_mets_out_ub = [NaN;NaN];
task_rxns = {'ctp_c + h2o_c <=> cdp_c + pi_c + h_c'};
task_rxns_lb = 1;
task3_success = test_MT(model, task_name, task_should_fail,...
                task_mets_in, task_mets_in_lb, task_mets_in_ub,... 
                task_mets_out, task_mets_out_lb,task_mets_out_ub,...
                task_rxns, task_rxns_lb, thermo_data);
table_MTs = [table_MTs;task_name num2cell(task3_success)];  
            
%Task4: Aerobic rephosphorylation of UTP
task_name = 'Aerobic rephosphorylation of UTP from glucose';
task_should_fail = 0;
task_mets_in = {'o2_e';'glc_e'};
task_mets_in_lb = [NaN;NaN];
task_mets_in_ub = [NaN;NaN];
task_mets_out = {'h2o_e';'co2_e'};
task_mets_out_lb = [NaN;NaN];
task_mets_out_ub = [NaN;NaN];
task_rxns = {'utp_c + h2o_c <=> udp_c + pi_c '};
task_rxns_lb = 1;
task4_success = test_MT(model, task_name, task_should_fail,...
                task_mets_in, task_mets_in_lb, task_mets_in_ub,... 
                task_mets_out, task_mets_out_lb,task_mets_out_ub,...
                task_rxns, task_rxns_lb, thermo_data);
table_MTs = [table_MTs;task_name num2cell(task4_success)]; 

%% De novo synthesis of nucleotides
%Task5: ATP de novo synthesis
task_name = 'ATP de novo synthesis';
task_should_fail = 0;
task_mets_in = {'o2_e';'glc_e';'nh4_e';'pi_e'};
task_mets_in_lb = [NaN;NaN;NaN;NaN];
task_mets_in_ub = [NaN;NaN;NaN;NaN];
task_mets_out = {'h2o_e';'co2_e';'atp_c'};
task_mets_out_lb = [NaN;NaN;1];
task_mets_out_ub = [NaN;NaN;NaN];
task_rxns = {'atp_c + h2o_c <=> adp_c + pi_c + h_c'};
task_rxns_lb = NaN;
task5_success = test_MT(model, task_name, task_should_fail,...
                task_mets_in, task_mets_in_lb, task_mets_in_ub,... 
                task_mets_out, task_mets_out_lb,task_mets_out_ub,...
                task_rxns, task_rxns_lb, thermo_data);
table_MTs = [table_MTs;task_name num2cell(task5_success)]; 

%Task6: CTP de novo synthesis
task_name = 'CTP de novo synthesis';
task_should_fail = 0;
task_mets_in = {'o2_e';'glc_e';'nh4_e';'pi_e'};
task_mets_in_lb = [NaN;NaN;NaN;NaN];
task_mets_in_ub = [NaN;NaN;NaN;NaN];
task_mets_out = {'h2o_e';'co2_e';'ctp_c'};
task_mets_out_lb = [NaN;NaN;1];
task_mets_out_ub = [NaN;NaN;NaN];
task_rxns = {'atp_c + h2o_c <=> adp_c + pi_c + h_c'};
task_rxns_lb = NaN;
task6_success = test_MT(model, task_name, task_should_fail,...
                task_mets_in, task_mets_in_lb, task_mets_in_ub,... 
                task_mets_out, task_mets_out_lb,task_mets_out_ub,...
                task_rxns, task_rxns_lb, thermo_data);
table_MTs = [table_MTs;task_name num2cell(task6_success)]; 

%Task7: GTP de novo synthesis
task_name = 'GTP de novo synthesis';
task_should_fail = 0;
task_mets_in = {'o2_e';'glc_e';'nh4_e';'pi_e'};
task_mets_in_lb = [NaN;NaN;NaN;NaN];
task_mets_in_ub = [NaN;NaN;NaN;NaN];
task_mets_out = {'h2o_e';'co2_e';'gtp_c'};
task_mets_out_lb = [NaN;NaN;1];
task_mets_out_ub = [NaN;NaN;NaN];
task_rxns = {'atp_c + h2o_c <=> adp_c + pi_c + h_c'};
task_rxns_lb = NaN;
task7_success = test_MT(model, task_name, task_should_fail,...
                task_mets_in, task_mets_in_lb, task_mets_in_ub,... 
                task_mets_out, task_mets_out_lb,task_mets_out_ub,...
                task_rxns, task_rxns_lb, thermo_data);
table_MTs = [table_MTs;task_name num2cell(task7_success)]; 

%Task8: UTP de novo synthesis
task_name = 'UTP de novo synthesis';
task_should_fail = 0;
task_mets_in = {'o2_e';'glc_e';'nh4_e';'pi_e'};
task_mets_in_lb = [NaN;NaN;NaN;NaN];
task_mets_in_ub = [NaN;NaN;NaN;NaN];
task_mets_out = {'h2o_e';'co2_e';'utp_c'};
task_mets_out_lb = [NaN;NaN;1];
task_mets_out_ub = [NaN;NaN;NaN];
task_rxns = {'atp_c + h2o_c <=> adp_c + pi_c + h_c'};
task_rxns_lb = NaN;
task8_success = test_MT(model, task_name, task_should_fail,...
                task_mets_in, task_mets_in_lb, task_mets_in_ub,... 
                task_mets_out, task_mets_out_lb,task_mets_out_ub,...
                task_rxns, task_rxns_lb, thermo_data);
table_MTs = [table_MTs;task_name num2cell(task8_success)]; 


%Task9: dATP de novo synthesis
task_name = 'dATP de novo synthesis';
task_should_fail = 0;
task_mets_in = {'o2_e';'glc_e';'nh4_e';'pi_e'};
task_mets_in_lb = [NaN;NaN;NaN;NaN];
task_mets_in_ub = [NaN;NaN;NaN;NaN];
task_mets_out = {'h2o_e';'co2_e';'datp_c'};
task_mets_out_lb = [NaN;NaN;1];
task_mets_out_ub = [NaN;NaN;NaN];
task_rxns = {'atp_c + h2o_c <=> adp_c + pi_c + h_c'};
task_rxns_lb = NaN;
task9_success = test_MT(model, task_name, task_should_fail,...
                task_mets_in, task_mets_in_lb, task_mets_in_ub,... 
                task_mets_out, task_mets_out_lb,task_mets_out_ub,...
                task_rxns, task_rxns_lb, thermo_data);
table_MTs = [table_MTs;task_name num2cell(task9_success)]; 

%Task10: dCTP de novo synthesis
task_name = 'dCTP de novo synthesis';
task_should_fail = 0;
task_mets_in = {'o2_e';'glc_e';'nh4_e';'pi_e'};
task_mets_in_lb = [NaN;NaN;NaN;NaN];
task_mets_in_ub = [NaN;NaN;NaN;NaN];
task_mets_out = {'h2o_e';'co2_e';'dctp_c'};
task_mets_out_lb = [NaN;NaN;1];
task_mets_out_ub = [NaN;NaN;NaN];
task_rxns = {'atp_c + h2o_c <=> adp_c + pi_c'};
task_rxns_lb = NaN;
task10_success = test_MT(model, task_name, task_should_fail,...
                task_mets_in, task_mets_in_lb, task_mets_in_ub,... 
                task_mets_out, task_mets_out_lb,task_mets_out_ub,...
                task_rxns, task_rxns_lb, thermo_data);
table_MTs = [table_MTs;task_name num2cell(task10_success)]; 

%Task11: dGTP de novo synthesis
task_name = 'dGTP de novo synthesis';
task_should_fail = 0;
task_mets_in = {'o2_e';'glc_e';'nh4_e';'pi_e'};
task_mets_in_lb = [NaN;NaN;NaN;NaN];
task_mets_in_ub = [NaN;NaN;NaN;NaN];
task_mets_out = {'h2o_e';'co2_e';'dgtp_c'};
task_mets_out_lb = [NaN;NaN;1];
task_mets_out_ub = [NaN;NaN;NaN];
task_rxns = {'atp_c + h2o_c <=> adp_c + pi_c'};
task_rxns_lb = NaN;
task11_success = test_MT(model, task_name, task_should_fail,...
                task_mets_in, task_mets_in_lb, task_mets_in_ub,... 
                task_mets_out, task_mets_out_lb,task_mets_out_ub,...
                task_rxns, task_rxns_lb, thermo_data);
table_MTs = [table_MTs;task_name num2cell(task11_success)]; 

%Task12: dTTP de novo synthesis
task_name = 'dTTP de novo synthesis';
task_should_fail = 0;
task_mets_in = {'o2_e';'glc_e';'nh4_e';'pi_e'};
task_mets_in_lb = [NaN;NaN;NaN;NaN];
task_mets_in_ub = [NaN;NaN;NaN;NaN];
task_mets_out = {'h2o_e';'co2_e';'dttp_c'};
task_mets_out_lb = [NaN;NaN;1];
task_mets_out_ub = [NaN;NaN;NaN];
task_rxns = {'atp_c + h2o_c <=> adp_c + pi_c'};
task_rxns_lb = NaN;
task12_success = test_MT(model, task_name, task_should_fail,...
                task_mets_in, task_mets_in_lb, task_mets_in_ub,... 
                task_mets_out, task_mets_out_lb,task_mets_out_ub,...
                task_rxns, task_rxns_lb, thermo_data);
table_MTs = [table_MTs;task_name num2cell(task12_success)]; 

%% Uptake of essential amino acids
%Task13: Histidine uptake
task_name = 'Histidine uptake';
task_should_fail = 0;
task_mets_in = {'his_L_e'};
task_mets_in_lb = 1;
task_mets_in_ub = 1;
task_mets_out = {'his_L_c'};
task_mets_out_lb = 1;
task_mets_out_ub = NaN;
task_rxns = {};
task_rxns_lb = [];
task13_success = test_MT(model, task_name, task_should_fail,...
                task_mets_in, task_mets_in_lb, task_mets_in_ub,... 
                task_mets_out, task_mets_out_lb,task_mets_out_ub,...
                task_rxns, task_rxns_lb, thermo_data);
table_MTs = [table_MTs;task_name num2cell(task13_success)]; 

%Task14: Isoleucine uptake
task_name = 'Isoleucine uptake';
task_should_fail = 0;
task_mets_in = {'ile_L_e'};
task_mets_in_lb = 1;
task_mets_in_ub = 1;
task_mets_out = {'ile_L_c'};
task_mets_out_lb = 1;
task_mets_out_ub = NaN;
task_rxns = {};
task_rxns_lb = [];
task14_success = test_MT(model, task_name, task_should_fail,...
                task_mets_in, task_mets_in_lb, task_mets_in_ub,... 
                task_mets_out, task_mets_out_lb,task_mets_out_ub,...
                task_rxns, task_rxns_lb, thermo_data);
table_MTs = [table_MTs;task_name num2cell(task14_success)]; 

%Task15: Leucine uptake
task_name = 'Leucine uptake';
task_should_fail = 0;
task_mets_in = {'leu_L_e'};
task_mets_in_lb = 1;
task_mets_in_ub = 1;
task_mets_out = {'leu_L_c'};
task_mets_out_lb = 1;
task_mets_out_ub = NaN;
task_rxns = {};
task_rxns_lb = [];
task15_success = test_MT(model, task_name, task_should_fail,...
                task_mets_in, task_mets_in_lb, task_mets_in_ub,... 
                task_mets_out, task_mets_out_lb,task_mets_out_ub,...
                task_rxns, task_rxns_lb, thermo_data);
table_MTs = [table_MTs;task_name num2cell(task15_success)]; 

%Task16: Lysine uptake
task_name = 'Lysine uptake';
task_should_fail = 0;
task_mets_in = {'lys_L_e'};
task_mets_in_lb = 1;
task_mets_in_ub = 1;
task_mets_out = {'lys_L_c'};
task_mets_out_lb = 1;
task_mets_out_ub = NaN;
task_rxns = {};
task_rxns_lb = [];
task16_success = test_MT(model, task_name, task_should_fail,...
                task_mets_in, task_mets_in_lb, task_mets_in_ub,... 
                task_mets_out, task_mets_out_lb,task_mets_out_ub,...
                task_rxns, task_rxns_lb, thermo_data);
table_MTs = [table_MTs;task_name num2cell(task16_success)]; 

%Task17: Methionine uptake
task_name = 'Methionine uptake';
task_should_fail = 0;
task_mets_in = {'met_L_e'};
task_mets_in_lb = 1;
task_mets_in_ub = 1;
task_mets_out = {'met_L_c'};
task_mets_out_lb = 1;
task_mets_out_ub = NaN;
task_rxns = {};
task_rxns_lb = [];
task17_success = test_MT(model, task_name, task_should_fail,...
                task_mets_in, task_mets_in_lb, task_mets_in_ub,... 
                task_mets_out, task_mets_out_lb,task_mets_out_ub,...
                task_rxns, task_rxns_lb, thermo_data);
table_MTs = [table_MTs;task_name num2cell(task17_success)]; 

%Task18: Phenylalanine uptake
task_name = 'Phenylalanine uptake';
task_should_fail = 0;
task_mets_in = {'phe_L_e'};
task_mets_in_lb = 1;
task_mets_in_ub = 1;
task_mets_out = {'phe_L_c'};
task_mets_out_lb = 1;
task_mets_out_ub = NaN;
task_rxns = {};
task_rxns_lb = [];
task18_success = test_MT(model, task_name, task_should_fail,...
                task_mets_in, task_mets_in_lb, task_mets_in_ub,... 
                task_mets_out, task_mets_out_lb,task_mets_out_ub,...
                task_rxns, task_rxns_lb, thermo_data);
table_MTs = [table_MTs;task_name num2cell(task18_success)]; 


%Task19: Threonine uptake
task_name = 'Threonine uptake';
task_should_fail = 0;
task_mets_in = {'thr_L_e'};
task_mets_in_lb = 1;
task_mets_in_ub = 1;
task_mets_out = {'thr_L_c'};
task_mets_out_lb = 1;
task_mets_out_ub = NaN;
task_rxns = {};
task_rxns_lb = [];
task19_success = test_MT(model, task_name, task_should_fail,...
                task_mets_in, task_mets_in_lb, task_mets_in_ub,... 
                task_mets_out, task_mets_out_lb,task_mets_out_ub,...
                task_rxns, task_rxns_lb, thermo_data);
table_MTs = [table_MTs;task_name num2cell(task19_success)]; 

%Task20: Tryptophan uptake
task_name = 'Tryptophan uptake';
task_should_fail = 0;
task_mets_in = {'trp_L_e'};
task_mets_in_lb = 1;
task_mets_in_ub = 1;
task_mets_out = {'trp_L_c'};
task_mets_out_lb = 1;
task_mets_out_ub = NaN;
task_rxns = {};
task_rxns_lb = [];
task20_success = test_MT(model, task_name, task_should_fail,...
                task_mets_in, task_mets_in_lb, task_mets_in_ub,... 
                task_mets_out, task_mets_out_lb,task_mets_out_ub,...
                task_rxns, task_rxns_lb, thermo_data);
table_MTs = [table_MTs;task_name num2cell(task20_success)]; 

%Task21: Valine uptake
task_name = 'Valine uptake';
task_should_fail = 0;
task_mets_in = {'val_L_e'};
task_mets_in_lb = 1;
task_mets_in_ub = 1;
task_mets_out = {'val_L_c'};
task_mets_out_lb = 1;
task_mets_out_ub = NaN;
task_rxns = {};
task_rxns_lb = [];
task21_success = test_MT(model, task_name, task_should_fail,...
                task_mets_in, task_mets_in_lb, task_mets_in_ub,... 
                task_mets_out, task_mets_out_lb,task_mets_out_ub,...
                task_rxns, task_rxns_lb, thermo_data);
table_MTs = [table_MTs;task_name num2cell(task21_success)]; 


%% De novo synthesis of key intermediates
%Task22: Glycerate 3-phosphate de novo synthesis
task_name = 'Glycerate 3-phosphate de novo synthesis';
task_should_fail = 0;
task_mets_in = {'o2_e';'glc_e';'pi_e'};
task_mets_in_lb = [NaN;NaN;NaN];
task_mets_in_ub = [NaN;NaN;NaN];
task_mets_out = {'3pg_c';'co2_e';'h2o_e'};
task_mets_out_lb = [1;NaN;NaN];
task_mets_out_ub = [NaN;NaN;NaN];
task_rxns = {};
task_rxns_lb =[];
task22_success = test_MT(model, task_name, task_should_fail,...
                task_mets_in, task_mets_in_lb, task_mets_in_ub,... 
                task_mets_out, task_mets_out_lb,task_mets_out_ub,...
                task_rxns, task_rxns_lb, thermo_data);
table_MTs = [table_MTs;task_name num2cell(task22_success)]; 

%Task23: Mitochondrial acetyl-CoA de novo synthesis
task_name = 'Mitochondrial acetyl-CoA de novo synthesis';
task_should_fail = 0;
task_mets_in = {'o2_e';'glc_e'};
task_mets_in_lb = [NaN;NaN];
task_mets_in_ub = [NaN;NaN];
task_mets_out = {'co2_e';'h2o_e'};
task_mets_out_lb = [NaN;NaN];
task_mets_out_ub = [NaN;NaN];
task_rxns = {'accoa_m <=> coa_m'};
task_rxns_lb = 1;
task23_success = test_MT(model, task_name, task_should_fail,...
                task_mets_in, task_mets_in_lb, task_mets_in_ub,... 
                task_mets_out, task_mets_out_lb,task_mets_out_ub,...
                task_rxns, task_rxns_lb, thermo_data);
table_MTs = [table_MTs;task_name num2cell(task23_success)]; 

%Task24: Mitochondrial AKG de novo synthesis
task_name = 'Mitochondrial AKG de novo synthesis';
task_should_fail = 0;
task_mets_in = {'o2_e';'glc_e'};
task_mets_in_lb = [NaN;NaN];
task_mets_in_ub = [NaN;NaN];
task_mets_out = {'akg_m';'co2_e';'h2o_e'};
task_mets_out_lb = [1;NaN;NaN];
task_mets_out_ub = [NaN;NaN;NaN];
task_rxns = {};
task_rxns_lb = [];
task24_success = test_MT(model, task_name, task_should_fail,...
                task_mets_in, task_mets_in_lb, task_mets_in_ub,... 
                task_mets_out, task_mets_out_lb,task_mets_out_ub,...
                task_rxns, task_rxns_lb, thermo_data);
table_MTs = [table_MTs;task_name num2cell(task24_success)]; 

%Task25: Erythrose 4-phosphate de novo synthesis
task_name = 'Erythrose 4-phosphate de novo synthesis';
task_should_fail = 0;
task_mets_in = {'o2_e';'glc_e';'pi_e'};
task_mets_in_lb = [NaN;NaN;NaN];
task_mets_in_ub = [NaN;NaN;NaN];
task_mets_out = {'e4p_c';'co2_e';'h2o_e'};
task_mets_out_lb = [1;NaN;NaN];
task_mets_out_ub = [NaN;NaN;NaN];
task_rxns = {};
task_rxns_lb =[];
task25_success = test_MT(model, task_name, task_should_fail,...
                task_mets_in, task_mets_in_lb, task_mets_in_ub,... 
                task_mets_out, task_mets_out_lb,task_mets_out_ub,...
                task_rxns, task_rxns_lb, thermo_data);
table_MTs = [table_MTs;task_name num2cell(task25_success)]; 

%Task26: Fructose 6-phosphate de novo synthesis
task_name = 'Fructose 6-phosphate de novo synthesis';
task_should_fail = 0;
task_mets_in = {'o2_e';'glc_e';'pi_e'};
task_mets_in_lb = [NaN;NaN;NaN];
task_mets_in_ub = [NaN;NaN;NaN];
task_mets_out = {'f6p_c';'co2_e';'h2o_e'};
task_mets_out_lb = [1;NaN;NaN];
task_mets_out_ub = [NaN;NaN;NaN];
task_rxns = {};
task_rxns_lb =[];
task26_success = test_MT(model, task_name, task_should_fail,...
                task_mets_in, task_mets_in_lb, task_mets_in_ub,... 
                task_mets_out, task_mets_out_lb,task_mets_out_ub,...
                task_rxns, task_rxns_lb, thermo_data);
table_MTs = [table_MTs;task_name num2cell(task26_success)]; 

%Task27: Glyceraldehyde 3-phosphate de novo synthesis
task_name = 'Glyceraldehyde 3-phosphate de novo synthesis';
task_should_fail = 0;
task_mets_in = {'o2_e';'glc_e';'pi_e'};
task_mets_in_lb = [NaN;NaN;NaN];
task_mets_in_ub = [NaN;NaN;NaN];
task_mets_out = {'g3p_c';'co2_e';'h2o_e'};
task_mets_out_lb = [1;NaN;NaN];
task_mets_out_ub = [NaN;NaN;NaN];
task_rxns = {};
task_rxns_lb =[];
task27_success = test_MT(model, task_name, task_should_fail,...
                task_mets_in, task_mets_in_lb, task_mets_in_ub,... 
                task_mets_out, task_mets_out_lb,task_mets_out_ub,...
                task_rxns, task_rxns_lb, thermo_data);
table_MTs = [table_MTs;task_name num2cell(task27_success)]; 

%Task28: Glucose 6-phosphate de novo synthesis
task_name = 'Glucose 6-phosphate de novo synthesis';
task_should_fail = 0;
task_mets_in = {'o2_e';'glc_e';'pi_e'};
task_mets_in_lb = [NaN;NaN;NaN];
task_mets_in_ub = [NaN;NaN;NaN];
task_mets_out = {'g6p_c';'co2_e';'h2o_e'};
task_mets_out_lb = [1;NaN;NaN];
task_mets_out_ub = [NaN;NaN;NaN];
task_rxns = {};
task_rxns_lb =[];
task28_success = test_MT(model, task_name, task_should_fail,...
                task_mets_in, task_mets_in_lb, task_mets_in_ub,... 
                task_mets_out, task_mets_out_lb,task_mets_out_ub,...
                task_rxns, task_rxns_lb, thermo_data);
table_MTs = [table_MTs;task_name num2cell(task28_success)]; 

%Task29: Mitochondrial oxaloacetate de novo synthesis
task_name = 'Mitochondrial oxaloacetate de novo synthesis';
task_should_fail = 0;
task_mets_in = {'o2_e';'glc_e';'pi_e'};
task_mets_in_lb = [NaN;NaN;NaN];
task_mets_in_ub = [NaN;NaN;NaN];
task_mets_out = {'oaa_m';'co2_e';'h2o_e'};
task_mets_out_lb = [1;NaN;NaN];
task_mets_out_ub = [NaN;NaN;NaN];
task_rxns = {};
task_rxns_lb =[];
task29_success = test_MT(model, task_name, task_should_fail,...
                task_mets_in, task_mets_in_lb, task_mets_in_ub,... 
                task_mets_out, task_mets_out_lb,task_mets_out_ub,...
                task_rxns, task_rxns_lb, thermo_data);
table_MTs = [table_MTs;task_name num2cell(task29_success)]; 

%Task30: Phosphoenolpyruvate de novo synthesis
task_name = 'Phosphoenolpyruvate de novo synthesis';
task_should_fail = 0;
task_mets_in = {'o2_e';'glc_e';'pi_e'};
task_mets_in_lb = [NaN;NaN;NaN];
task_mets_in_ub = [NaN;NaN;NaN];
task_mets_out = {'pep_c';'co2_e';'h2o_e'};
task_mets_out_lb = [1;NaN;NaN];
task_mets_out_ub = [NaN;NaN;NaN];
task_rxns = {};
task_rxns_lb =[];
task30_success = test_MT(model, task_name, task_should_fail,...
                task_mets_in, task_mets_in_lb, task_mets_in_ub,... 
                task_mets_out, task_mets_out_lb,task_mets_out_ub,...
                task_rxns, task_rxns_lb, thermo_data);
table_MTs = [table_MTs;task_name num2cell(task30_success)]; 

%Task31: Pyruvate de novo synthesis
task_name = 'Pyruvate de novo synthesis';
task_should_fail = 0;
task_mets_in = {'o2_e';'glc_e';'pi_e'};
task_mets_in_lb = [NaN;NaN;NaN];
task_mets_in_ub = [NaN;NaN;NaN];
task_mets_out = {'pyr_c';'co2_e';'h2o_e'};
task_mets_out_lb = [1;NaN;NaN];
task_mets_out_ub = [NaN;NaN;NaN];
task_rxns = {};
task_rxns_lb =[];
task31_success = test_MT(model, task_name, task_should_fail,...
                task_mets_in, task_mets_in_lb, task_mets_in_ub,... 
                task_mets_out, task_mets_out_lb,task_mets_out_ub,...
                task_rxns, task_rxns_lb, thermo_data);
table_MTs = [table_MTs;task_name num2cell(task31_success)]; 

%Task32: Ribose 5-phosphate de novo synthesis
task_name = 'Ribose 5-phosphate de novo synthesis';
task_should_fail = 0;
task_mets_in = {'o2_e';'glc_e';'pi_e'};
task_mets_in_lb = [NaN;NaN;NaN];
task_mets_in_ub = [NaN;NaN;NaN];
task_mets_out = {'r5p_c';'co2_e';'h2o_e'};
task_mets_out_lb = [1;NaN;NaN];
task_mets_out_ub = [NaN;NaN;NaN];
task_rxns = {};
task_rxns_lb =[];
task32_success = test_MT(model, task_name, task_should_fail,...
                task_mets_in, task_mets_in_lb, task_mets_in_ub,... 
                task_mets_out, task_mets_out_lb,task_mets_out_ub,...
                task_rxns, task_rxns_lb, thermo_data);
table_MTs = [table_MTs;task_name num2cell(task32_success)]; 

%Task33: Mitochondrial succinyl-CoA de novo synthesis
task_name = 'Mitochondrial succinyl-CoA de novo synthesis';
task_should_fail = 0;
task_mets_in = {'o2_e';'glc_e';'pi_e'};
task_mets_in_lb = [NaN;NaN;NaN];
task_mets_in_ub = [NaN;NaN;NaN];
task_mets_out = {'co2_e';'h2o_e'};
task_mets_out_lb = [NaN;NaN];
task_mets_out_ub = [NaN;NaN];
task_rxns = {'succoa_m <=> coa_m'};
task_rxns_lb = 1;
task33_success = test_MT(model, task_name, task_should_fail,...
                task_mets_in, task_mets_in_lb, task_mets_in_ub,... 
                task_mets_out, task_mets_out_lb,task_mets_out_ub,...
                task_rxns, task_rxns_lb, thermo_data);
table_MTs = [table_MTs;task_name num2cell(task33_success)]; 

%% De novo synthesis of other compounds
%Task34: Cholesterol de novo synthesis
task_name = 'Cholesterol de novo synthesis';
task_should_fail = 0;
task_mets_in = {'o2_e';'glc_e'};
task_mets_in_lb = [1;1];
task_mets_in_ub = [NaN;NaN];
task_mets_out = {'chsterol_c';'co2_e';'h2o_e'};
task_mets_out_lb = [1;NaN;NaN];
task_mets_out_ub = [NaN;NaN;NaN];
task_rxns = {};
task_rxns_lb = [];
task34_success = test_MT(model, task_name, task_should_fail,...
                task_mets_in, task_mets_in_lb, task_mets_in_ub,... 
                task_mets_out, task_mets_out_lb,task_mets_out_ub,...
                task_rxns, task_rxns_lb, thermo_data);
table_MTs = [table_MTs;task_name num2cell(task34_success)]; 

%% Protein turnover
%Task35: Protein synthesis from AAs %%PROBLEM!!!
task_name = 'Protein synthesis from AAs';
task_should_fail = 0;
task_mets_in = {'o2_e';'glc_e';'nh4_e';'h2o_e';'arg_L_e';'his_L_e';'lys_L_e';...
                'met_L_e';'phe_L_e';'trp_L_e';'tyr_L_e';'ala_L_e';'gly_e';...
                'ser_L_e';'thr_L_e';'asp_L_e';'glu_L_e';'asn_L_e';'gln_L_e';...
                'ile_L_e';'leu_L_e';'pro_L_e';'val_L_e';'cys_L_e'};
task_mets_in_lb = [NaN;NaN;NaN;NaN;NaN;NaN;NaN;NaN;NaN;NaN;NaN;NaN;NaN;NaN;...
                   NaN;NaN;NaN;NaN;NaN;NaN;NaN;NaN;NaN];%;NaN];
task_mets_in_ub = [NaN;NaN;NaN;NaN;NaN;NaN;NaN;NaN;NaN;NaN;NaN;NaN;NaN;NaN;...
                   NaN;NaN;NaN;NaN;NaN;NaN;NaN;NaN;NaN];%;NaN];
task_mets_out = {'HC00001_c';'co2_e';'h2o_e';'urea_e'};
task_mets_out_lb = [1;NaN;NaN;NaN;NaN];
task_mets_out_ub = [NaN;NaN;NaN;NaN;NaN];
task_rxns = {};
task_rxns_lb = [];
task35_success = test_MT(model, task_name, task_should_fail,...
                task_mets_in, task_mets_in_lb, task_mets_in_ub,... 
                task_mets_out, task_mets_out_lb,task_mets_out_ub,...
                task_rxns, task_rxns_lb, thermo_data);
table_MTs = [table_MTs;task_name num2cell(task35_success)]; 

%% Electron transport chain and TCA
%Task36: Oxidative phosphorylation
task_name = 'Oxidative phosphorylation';
task_should_fail = 0;
task_mets_in = {'succ_m';'nadh_m';'h_m';'o2_e'};
task_mets_in_lb = [1;1;1;NaN];
task_mets_in_ub = [1;1;1;NaN];
task_mets_out = {'fum_m';'nad_m';'h2o_e'};
task_mets_out_lb = [1;1;NaN];
task_mets_out_ub = [1;1;NaN];
task_rxns = {'atp_m + h2o_m <=> adp_m + pi_m'};
task_rxns_lb = 1;
task36_success = test_MT(model, task_name, task_should_fail,...
                task_mets_in, task_mets_in_lb, task_mets_in_ub,... 
                task_mets_out, task_mets_out_lb,task_mets_out_ub,...
                task_rxns, task_rxns_lb, thermo_data);
table_MTs = [table_MTs;task_name num2cell(task36_success)]; 

%Task37: Oxidative decarboxylation
task_name = 'Oxidative decarboxylation';
task_should_fail = 0;
task_mets_in = {'pyr_m';'nad_m';'coa_m'};
task_mets_in_lb = [1;1;1];
task_mets_in_ub = [1;1;1];
task_mets_out = {'accoa_m';'nadh_m';'h_m';'co2_e'};
task_mets_out_lb = [1;1;1;1];
task_mets_out_ub = [1;1;1;1];
task_rxns = {};
task_rxns_lb = [];
task37_success = test_MT(model, task_name, task_should_fail,...
                task_mets_in, task_mets_in_lb, task_mets_in_ub,... 
                task_mets_out, task_mets_out_lb,task_mets_out_ub,...
                task_rxns, task_rxns_lb, thermo_data);
table_MTs = [table_MTs;task_name num2cell(task37_success)]; 

%Task38: Krebs cycle NADH
task_name = 'Krebs cycle NADH';
task_should_fail = 0;
task_mets_in = {'accoa_m';'gdp_m';'q10_m';'nad_m';'pi_m';'h2o_e'};
task_mets_in_lb = [1;1;1;3;NaN;NaN];
task_mets_in_ub = [1;1;1;3;NaN;NaN];
task_mets_out = {'coa_m';'q10h2_m';'gtp_m';'nadh_m';'co2_e';'h_c';'h_m'};
task_mets_out_lb = [1;1;1;3;NaN;NaN;NaN];
task_mets_out_ub = [1;1;1;3;NaN;NaN;NaN];
task_rxns = {};
task_rxns_lb = [];
task38_success = test_MT(model, task_name, task_should_fail,...
                task_mets_in, task_mets_in_lb, task_mets_in_ub,... 
                task_mets_out, task_mets_out_lb,task_mets_out_ub,...
                task_rxns, task_rxns_lb, thermo_data);
table_MTs = [table_MTs;task_name num2cell(task38_success)]; 

%Task39: Ubiquinol-to-proton
task_name = 'Ubiquinol-to-proton';
task_should_fail = 0;
task_mets_in = {'q10h2_m';'o2_e';'h_m'};
task_mets_in_lb = [1;NaN;NaN];
task_mets_in_ub = [1;NaN;NaN];
task_mets_out = {'q10_m';'h2o_e';'h_c'};
task_mets_out_lb = [1;NaN;6];
task_mets_out_ub = [1;NaN;NaN];
task_rxns = {};
task_rxns_lb = [];
task39_success = test_MT(model, task_name, task_should_fail,...
                task_mets_in, task_mets_in_lb, task_mets_in_ub,... 
                task_mets_out, task_mets_out_lb,task_mets_out_ub,...
                task_rxns, task_rxns_lb, thermo_data);
table_MTs = [table_MTs;task_name num2cell(task39_success)]; 

%Task40: Ubiquinol-to-ATP
task_name = 'Ubiquinol-to-ATP';
task_should_fail = 0;
task_mets_in = {'q10h2_m';'o2_e';'h_m'};
task_mets_in_lb = [1;NaN;NaN];
task_mets_in_ub = [1;NaN;NaN];
task_mets_out = {'q10_m';'h2o_e';'h_m';'h_c'};
task_mets_out_lb = [1;NaN;NaN;NaN];
task_mets_out_ub = [1;NaN;NaN;NaN];
task_rxns = {'atp_m + h2o_m <=> adp_m + pi_m'};
task_rxns_lb = 1;
task40_success = test_MT(model, task_name, task_should_fail,...
                task_mets_in, task_mets_in_lb, task_mets_in_ub,... 
                task_mets_out, task_mets_out_lb,task_mets_out_ub,...
                task_rxns, task_rxns_lb, thermo_data);
table_MTs = [table_MTs;task_name num2cell(task40_success)]; 

%% Beta oxidation of fatty acids
% Task41: Beta oxidation of saturated FA
task_name = 'Beta oxidation of saturated FA';
task_should_fail = 0;
task_mets_in = {'ocdca_e';'o2_e'};
task_mets_in_lb = [1;NaN];
task_mets_in_ub = [1;NaN];
task_mets_out = {'h2o_e';'co2_e'};
task_mets_out_lb = [NaN;NaN];
task_mets_out_ub = [NaN;NaN];
task_rxns = {'atp_c + h2o_c <=> adp_c + pi_c + h_c'};
task_rxns_lb = 1;
task41_success = test_MT(model, task_name, task_should_fail,...
                task_mets_in, task_mets_in_lb, task_mets_in_ub,... 
                task_mets_out, task_mets_out_lb,task_mets_out_ub,...
                task_rxns, task_rxns_lb, thermo_data);
table_MTs = [table_MTs;task_name num2cell(task41_success)]; 

% Task42: Beta oxidation of long-chain FA
task_name = 'Beta oxidation of long-chain FA';
task_should_fail = 0;
task_mets_in = {'M00315_e';'o2_e'};
task_mets_in_lb = [1;NaN];
task_mets_in_ub = [1;NaN];
task_mets_out = {'h2o_e';'co2_e'};
task_mets_out_lb = [NaN;NaN];
task_mets_out_ub = [NaN;NaN];
task_rxns = {'atp_c + h2o_c <=> adp_c + pi_c + h_c'};
task_rxns_lb = 1;
task42_success = test_MT(model, task_name, task_should_fail,...
                task_mets_in, task_mets_in_lb, task_mets_in_ub,... 
                task_mets_out, task_mets_out_lb,task_mets_out_ub,...
                task_rxns, task_rxns_lb, thermo_data);
table_MTs = [table_MTs;task_name num2cell(task42_success)]; 

% Task43: Beta oxidation of odd-chain FA
task_name = 'Beta oxidation of odd-chain FA';
task_should_fail = 0;
task_mets_in = {'hpdca_e';'o2_e'};
task_mets_in_lb = [1;NaN];
task_mets_in_ub = [1;NaN];
task_mets_out = {'h2o_e';'co2_e'};
task_mets_out_lb = [NaN;NaN];
task_mets_out_ub = [NaN;NaN];
task_rxns = {'atp_c + h2o_c <=> adp_c + pi_c + h_c'};
task_rxns_lb = 1;
task43_success = test_MT(model, task_name, task_should_fail,...
                task_mets_in, task_mets_in_lb, task_mets_in_ub,... 
                task_mets_out, task_mets_out_lb,task_mets_out_ub,...
                task_rxns, task_rxns_lb, thermo_data);
table_MTs = [table_MTs;task_name num2cell(task43_success)]; 

% Task44: Beta oxidation of unsaturated fatty acid (n-9)
task_name = 'Beta oxidation of unsaturated fatty acid (n-9)';
task_should_fail = 0;
task_mets_in = {'M03153_e';'o2_e'};
task_mets_in_lb = [1;NaN];
task_mets_in_ub = [1;NaN];
task_mets_out = {'h2o_e';'co2_e'};
task_mets_out_lb = [NaN;NaN];
task_mets_out_ub = [NaN;NaN];
task_rxns = {'atp_c + h2o_c <=> adp_c + pi_c + h_c'};
task_rxns_lb = 1;
task44_success = test_MT(model, task_name, task_should_fail,...
                task_mets_in, task_mets_in_lb, task_mets_in_ub,... 
                task_mets_out, task_mets_out_lb,task_mets_out_ub,...
                task_rxns, task_rxns_lb, thermo_data);
table_MTs = [table_MTs;task_name num2cell(task44_success)]; 


% Task45: Beta oxidation of unsaturated fatty acid (n-6)
task_name = 'Beta oxidation of unsaturated fatty acid (n-6)';
task_should_fail = 0;
task_mets_in = {'lnlc_e';'o2_e'};
task_mets_in_lb = [1;NaN];
task_mets_in_ub = [1;NaN];
task_mets_out = {'h2o_e';'co2_e'};
task_mets_out_lb = [NaN;NaN];
task_mets_out_ub = [NaN;NaN];
task_rxns = {'atp_c + h2o_c <=> adp_c + pi_c + h_c'};
task_rxns_lb = 1;
task45_success = test_MT(model, task_name, task_should_fail,...
                task_mets_in, task_mets_in_lb, task_mets_in_ub,... 
                task_mets_out, task_mets_out_lb,task_mets_out_ub,...
                task_rxns, task_rxns_lb, thermo_data);
table_MTs = [table_MTs;task_name num2cell(task45_success)]; 

% Task46: Uptake and beta oxidation of all NEFAs
task_name = 'Uptake and beta oxidation of all NEFAs';
task_should_fail = 0;
task_mets_in = {'lnlc_e';'lnlnca_e';'o2_e'};
task_mets_in_lb = [1;1;NaN];
task_mets_in_ub = [1;1;NaN];
task_mets_out = {'h2o_e';'co2_e'};
task_mets_out_lb = [NaN;NaN];
task_mets_out_ub = [NaN;NaN];
task_rxns = {'atp_c + h2o_c <=> adp_c + pi_c + h_c'};
task_rxns_lb = 1;
task46_success = test_MT(model, task_name, task_should_fail,...
                task_mets_in, task_mets_in_lb, task_mets_in_ub,... 
                task_mets_out, task_mets_out_lb,task_mets_out_ub,...
                task_rxns, task_rxns_lb, thermo_data);
table_MTs = [table_MTs;task_name num2cell(task46_success)]; 

%% De novo synthesis of phospholipids
% Task47: Choline uptake
task_name = 'Choline uptake';
task_should_fail = 0;
task_mets_in = {'chol_e'};
task_mets_in_lb = 1;
task_mets_in_ub = 1;
task_mets_out = {'chol_c'};
task_mets_out_lb = 1;
task_mets_out_ub = 1;
task_rxns = {};
task_rxns_lb = [];
task47_success = test_MT(model, task_name, task_should_fail,...
                task_mets_in, task_mets_in_lb, task_mets_in_ub,... 
                task_mets_out, task_mets_out_lb,task_mets_out_ub,...
                task_rxns, task_rxns_lb, thermo_data);
table_MTs = [table_MTs;task_name num2cell(task47_success)]; 

% Task48: Inositol uptake
task_name = 'Inositol uptake';
task_should_fail = 0;
task_mets_in = {'inost_e'};
task_mets_in_lb = 1;
task_mets_in_ub = 1;
task_mets_out = {'inost_c'};
task_mets_out_lb = 1;
task_mets_out_ub = 1;
task_rxns = {};
task_rxns_lb = [];
task48_success = test_MT(model, task_name, task_should_fail,...
                task_mets_in, task_mets_in_lb, task_mets_in_ub,... 
                task_mets_out, task_mets_out_lb,task_mets_out_ub,...
                task_rxns, task_rxns_lb, thermo_data);
table_MTs = [table_MTs;task_name num2cell(task48_success)]; 

% Task49: Phosphatidylcholine de novo synthesis
task_name = 'Phosphatidylcholine de novo synthesis';
task_should_fail = 0;
task_mets_in = {'chol_e';'glc_e';'o2_e';'pi_e'};%;'lnlc_e';'lnlnca_e'};
task_mets_in_lb = [NaN;NaN;NaN;NaN;NaN;NaN];
task_mets_in_ub = [NaN;NaN;NaN;NaN;NaN;NaN];
task_mets_out = {'pchol_hs_c';'h2o_e';'co2_e'};
task_mets_out_lb = [1;NaN;NaN];
task_mets_out_ub = [1;NaN;NaN];
task_rxns = {'atp_c + h2o_c <=> adp_c + pi_c'};
task_rxns_lb = NaN;
task49_success = test_MT(model, task_name, task_should_fail,...
                task_mets_in, task_mets_in_lb, task_mets_in_ub,... 
                task_mets_out, task_mets_out_lb,task_mets_out_ub,...
                task_rxns, task_rxns_lb, thermo_data);
table_MTs = [table_MTs;task_name num2cell(task49_success)]; 

% Task50: Phosphatidylethanolamine de novo synthesis
task_name = 'Phosphatidylethanolamine de novo synthesis';
task_should_fail = 0;
task_mets_in = {'ser_L_e';'glc_e';'o2_e';'pi_e'};%;'lnlc_e';'lnlnca_e';'ethanolamine};
task_mets_in_lb = [NaN;NaN;NaN;NaN;NaN;NaN];
task_mets_in_ub = [NaN;NaN;NaN;NaN;NaN;NaN];
task_mets_out = {'pe_hs_c';'h2o_e';'co2_e'};
task_mets_out_lb = [1;NaN;NaN];
task_mets_out_ub = [1;NaN;NaN];
task_rxns = {'atp_c + h2o_c <=> adp_c + pi_c'};
task_rxns_lb = NaN;
task50_success = test_MT(model, task_name, task_should_fail,...
                task_mets_in, task_mets_in_lb, task_mets_in_ub,... 
                task_mets_out, task_mets_out_lb,task_mets_out_ub,...
                task_rxns, task_rxns_lb, thermo_data);
table_MTs = [table_MTs;task_name num2cell(task50_success)]; 

% Task51: Phosphatidylserine de novo synthesis
task_name = 'Phosphatidylserine de novo synthesis';
task_should_fail = 0;
task_mets_in = {'ser_L_c';'glc_e';'o2_e';'pi_e'};%;'lnlc_e';'lnlnca_e'};
task_mets_in_lb = [NaN;NaN;NaN;NaN];
task_mets_in_ub = [NaN;NaN;NaN;NaN];
task_mets_out = {'ps_hs_c';'h2o_e';'co2_e'};
task_mets_out_lb = [1;NaN;NaN];
task_mets_out_ub = [1;NaN;NaN];
task_rxns = {'atp_c + h2o_c <=> adp_c + pi_c'};
task_rxns_lb = NaN;
task51_success = test_MT(model, task_name, task_should_fail,...
                task_mets_in, task_mets_in_lb, task_mets_in_ub,... 
                task_mets_out, task_mets_out_lb,task_mets_out_ub,...
                task_rxns, task_rxns_lb, thermo_data);
table_MTs = [table_MTs;task_name num2cell(task51_success)]; 

% Task52: Phosphatidylinositol de novo synthesis
task_name = 'Phosphatidylinositol de novo synthesis';
task_should_fail = 0;
task_mets_in = {'glc_e';'o2_e';'pi_e'};%;'lnlc_e';'lnlnca_e';'inost_e'};
task_mets_in_lb = [NaN;NaN;NaN;NaN];
task_mets_in_ub = [NaN;NaN;NaN;NaN];
task_mets_out = {'pail_hs_c';'h2o_e';'co2_e'};
task_mets_out_lb = [NaN;NaN;NaN];
task_mets_out_ub = [NaN;NaN;NaN];
task_rxns = {'atp_c + h2o_c <=> adp_c + pi_c'};
task_rxns_lb = NaN;
task52_success = test_MT(model, task_name, task_should_fail,...
                task_mets_in, task_mets_in_lb, task_mets_in_ub,... 
                task_mets_out, task_mets_out_lb,task_mets_out_ub,...
                task_rxns, task_rxns_lb, thermo_data);
table_MTs = [table_MTs;task_name num2cell(task52_success)]; 


%% Vitamins and co-factors
% Task53: Thiamin phosphorylation to TPP
task_name = 'Thiamin phosphorylation to TPP';
task_should_fail = 0;
task_mets_in = {'thm_e';'pi_e';'o2_e';'h2o_e'};
task_mets_in_lb = [NaN;NaN;NaN;NaN];
task_mets_in_ub = [NaN;NaN;NaN;NaN];
task_mets_out = {'thm-PP_c';'h2o_e'};
task_mets_out_lb = [1;NaN];
task_mets_out_ub = [NaN;NaN];
task_rxns = {'atp_c + h2o_c <=> adp_c + pi_c'};
task_rxns_lb = NaN;
task53_success = test_MT(model, task_name, task_should_fail,...
                task_mets_in, task_mets_in_lb, task_mets_in_ub,... 
                task_mets_out, task_mets_out_lb,task_mets_out_ub,...
                task_rxns, task_rxns_lb, thermo_data);
table_MTs = [table_MTs;task_name num2cell(task53_success)]; 

% Task54: Coenzyme A synthesis from pantothenate
task_name = 'Coenzyme A synthesis from pantothenate';
task_should_fail = 0;
task_mets_in = {'pnto_R_e';'cys_L_e';'glc_e';'o2_e';'h2o_e';'pi_e';'nh4_e'};
task_mets_in_lb = [NaN;NaN;NaN;NaN;NaN;NaN;NaN];
task_mets_in_ub = [NaN;NaN;NaN;NaN;NaN;NaN;NaN];
task_mets_out = {'coa_c';'h2o_e';'co2_e'};
task_mets_out_lb = [1;NaN;NaN];
task_mets_out_ub = [1;NaN;NaN];
task_rxns = {'atp_c + h2o_c <=> adp_c + pi_c'};
task_rxns_lb = NaN;
task54_success = test_MT(model, task_name, task_should_fail,...
                task_mets_in, task_mets_in_lb, task_mets_in_ub,... 
                task_mets_out, task_mets_out_lb,task_mets_out_ub,...
                task_rxns, task_rxns_lb, thermo_data);
table_MTs = [table_MTs;task_name num2cell(task54_success)]; 

% Task55: FAD synthesis from riboflavin
task_name = 'FAD synthesis from riboflavin';
task_should_fail = 0;
task_mets_in = {'o2_e';'glc_D_e';'nh4_e';'pi_e';'ribflv_e'};
task_mets_in_lb = [NaN;NaN;NaN;NaN;NaN];
task_mets_in_ub = [NaN;NaN;NaN;NaN;NaN];
task_mets_out = {'fad_c';'h2o_e';'co2_e'};
task_mets_out_lb = [1;NaN;NaN];
task_mets_out_ub = [1;NaN;NaN];
task_rxns = {};
task_rxns_lb = [];
task55_success = test_MT(model, task_name, task_should_fail,...
                task_mets_in, task_mets_in_lb, task_mets_in_ub,... 
                task_mets_out, task_mets_out_lb,task_mets_out_ub,...
                task_rxns, task_rxns_lb, thermo_data);
table_MTs = [table_MTs;task_name num2cell(task55_success)]; 

% Task56: Heme biosynthesis
task_name = 'Heme biosynthesis';
task_should_fail = 0;
task_mets_in = {'o2_e';'nh4_e';'glc_e';'fe2_e'};
task_mets_in_lb = [NaN;NaN;NaN;NaN];
task_mets_in_ub = [NaN;NaN;NaN;NaN];
task_mets_out = {'pheme_c';'h2o_e';'co2_e'};
task_mets_out_lb = [1;NaN;NaN];
task_mets_out_ub = [1;NaN;NaN];
task_rxns = {'atp_c + h2o_c <=> adp_c + pi_c'};
task_rxns_lb = NaN;
task56_success = test_MT(model, task_name, task_should_fail,...
                task_mets_in, task_mets_in_lb, task_mets_in_ub,... 
                task_mets_out, task_mets_out_lb,task_mets_out_ub,...
                task_rxns, task_rxns_lb, thermo_data);
table_MTs = [table_MTs;task_name num2cell(task56_success)]; 


%% Growth
% Task 57: Growth in Ham's media
% task_name = 'Growth in Ham''s media';
% task_should_fail = 0;
% task_mets_in = {'arg_L_e';'his_L_e';'lys_L_e';'met_L_e';'phe_L_e';'trp_L_e';...
%                 'tyr_L_e';'ala_L_e';'gly_e';'ser_L_e';'thr_L_e';'asp_L_e';...
%                 'glu_L_e';'asn_L_e';'gln_L_e';'ile_L_e';'leu_L_e';'pro_L_e';...
%                 'val_L_e';'cys_L_e';...%'thm_e';'hxan_e';'fol_e';'btn_e';'pnto_R_e';...
%                 %'chol_e';'inost_e';'ncam_e';'pydxn_e';'ribflv_e';'thymd_e';...
%                 %'aqcobal_e';'lipoate_e';...
%                 'glc_e';'so4_e';'lnlc_e';'lnlnca_e';'lnlc_e';'lnlnca_e';...
%                 'o2_e';'h2o_e';'retn_e';'fe2_e';'pi_e'};%;...
%                 %'apocytochrome-C_c';'apoA1_c'};
% task_mets_in_lb = [NaN;NaN;NaN;NaN;NaN;NaN;NaN;NaN;NaN;NaN;NaN;NaN;NaN;NaN;NaN;NaN;NaN;NaN;NaN;NaN];
% task_mets_in_ub = [NaN;NaN;NaN;NaN;NaN;NaN;NaN;NaN;NaN;NaN;NaN;NaN;NaN;NaN;NaN;NaN;NaN;NaN;NaN;NaN];
% task_mets_out = {'h2o_e';'co2_e';'urea_e';'H2S_e'};
% task_mets_out_lb = [NaN;NaN;NaN;NaN];
% task_mets_out_ub = [NaN;NaN;NaN;NaN];
% task_rxns = {'atp_c + h2o_c <=> adp_c + pi_c'};
% task_rxns_lb = NaN;
% task56_success = test_MT(model, task_name, task_should_fail,...
%                 task_mets_in, task_mets_in_lb, task_mets_in_ub,... 
%                 task_mets_out, task_mets_out_lb,task_mets_out_ub,...
%                 task_rxns, task_rxns_lb, DB_AlbertyUpdate);
% table_MTs = [table_MTs;task_name num2cell(task56_success)]; 
