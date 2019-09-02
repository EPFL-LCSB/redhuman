function task_success = test_MT(model, task_name, task_should_fail,...
                                task_mets_in, task_mets_in_lb, task_mets_in_ub,... 
                                task_mets_out, task_mets_out_lb,task_mets_out_ub,...
                                task_rxns, task_rxns_lb, thermo_data)
                            
% function to test metabolic task
% INPUTS: 
% - model: TFA model
% - task_name: name for the task
% - task_met: metabolite that we want to produce in the specific task
% - task_should_fail: 1 if task should fail, 0 otherwise
% - task_mets_in: set of metabolites allowed to be uptaken from the medium.
%                 If empty, medium left as in the model
% - task_mets_in_lb: lower bound for task_mets_in
% - task_mets_in_ub: upper bound for task_mets_in
% - task_mets_out: set of metabolites allowed to be secreted from the medium.
%                 If empty, medium left as in the model
% - task_mets_out_lb: lower bound for task_mets_out
% - task_mets_out_lb: upper bound for task_mets_out
% - task_rxns: reactions that are forced to carry flux during the specific
%              task
% - task_rxns_lb: lower bound for task_rxns
% - task_rxns_bounds: vector with bounds for the reactions. If empty, bounds
%                     left as in the model
% 
% OUTPUTS
% - task_success: 1 if success, 0 otherwise


% make sure there are no reactions being forced
model.var_lb(1:2*length(model.rxns)) = 0;

drains = extract_drains(model);
drains_mets_in  = strcat('EX_', task_mets_in);
drains_mets_out = strcat('EX_', task_mets_out);


task_mets_not_drain_model = task_mets_in(find(ismember(drains_mets_in,drains)==0));
if ~isempty(task_mets_not_drain_model)
    for i = 1:length(task_mets_not_drain_model)
        comps{i,1} = task_mets_not_drain_model{i}(end);
    end
    if isempty(setdiff(comps,{'e'}))
        disp('The uptake reactions for some metabolites are not part of the model')
        task_success = NaN;
        return
    else
        new_drains = task_mets_in(find(ismember(drains_mets_in,drains)==0));
        ind_met_out = find(ismember(drains_mets_in,drains)==0);
        for i = 1:length(new_drains)
            model = addDemandReaction(model,task_mets_in(ind_met_out(i)),1,-100,0);
            model.rxns(end) = strcat('EX_',task_mets_in(ind_met_out(i)));
        end
    end
end

task_mets_not_drain_model = task_mets_out(find(ismember(drains_mets_out,drains)==0));
if  ~isempty(task_mets_not_drain_model)
    for i = 1:length(task_mets_not_drain_model)
        comps{i,1} = task_mets_not_drain_model{i}(end);
    end
    if isempty(setdiff(comps,{'e'}))
    disp('The secretion reactions for some metabolites are not part of the model')
    task_success = NaN;
    return
    else
        new_drains = task_mets_out(find(ismember(drains_mets_out,drains)==0));
        ind_met_out = find(ismember(drains_mets_out,drains)==0);
        for i = 1:length(new_drains)
            if ~ismember(strcat('EX_',task_mets_out(ind_met_out(i))),model.rxns)
                model = addDemandReaction(model,task_mets_out(ind_met_out(i)),1,0,100);
                model.rxns(end) = strcat('EX_',task_mets_out(ind_met_out(i)));
            end
        end
    end
end

for i = 1:length(task_rxns)
    [model,rxnIDexists] = addReaction(model,strcat('task_rxn',num2str(i)),task_rxns{i});
    if ~isempty(rxnIDexists)
        model.rxns(rxnIDexists) = {strcat('task_rxn',num2str(i))};
    end
    if ~isnan(task_rxns_lb(i))
        model.lb(end) = task_rxns_lb(i);
    end
end

%fix to be able to convert to thermo.
if length(model.mets)>length(model.metSEEDID)
    num_new_mets = length(model.mets)-length(model.metSEEDID);
    ii = length(model.metSEEDID);
    new_mets = model.mets(ii+1:ii+num_new_mets);
    for i = 1:length(new_mets)
        model.metCompSymbol{ii+i} = new_mets{i}(end);
        if ismember({new_mets{i}(1:end-2)},thermo_data.compound.modelMets)
            model.metSEEDID(ii+i) = model.metSEEDID(find(ismember({new_mets{i}(1:end-2)},thermo_data.compound.modelMets)));
        else
            model.metSEEDID(ii+i) = {'NA'};
        end
    end
    
end

% convert to thermo
model_orig = model;
model = TFAprepconv(model,thermo_data);
% recover original bounds
[~,ii] = ismember(model_orig.varNames,model.varNames);
model.var_lb(ii(find(ii))) = model_orig.var_lb(find(ii));
model.var_ub(ii(find(ii))) = model_orig.var_ub(find(ii));

uptk_to_block = setdiff(strcat('R_',drains),strcat('R_EX_',task_mets_in));
model.var_ub(find(ismember(model.varNames,uptk_to_block))) = 0;
model.var_ub(find(ismember(model.varNames,strcat('R_EX_',task_mets_in)))) = 100;
scrtns_to_block = setdiff(strcat('F_',drains),strcat('F_EX_',task_mets_out));
model.var_ub(find(ismember(model.varNames,scrtns_to_block))) = 0;
model.var_ub(find(ismember(model.varNames,strcat('F_EX_',task_mets_out)))) = 100;
for i = 1:length(task_mets_in)
    if ~isnan(task_mets_in_lb(i))
        model.var_lb(find(ismember(model.varNames,strcat('R_EX_',task_mets_in(i))))) = task_mets_in_lb(i);
    end
    
    if ~isnan(task_mets_in_ub(i))
        model.var_ub(find(ismember(model.varNames,strcat('R_EX_',task_mets_in(i))))) = task_mets_in_ub(i);
    end 
end

for i = 1:length(task_mets_out)
    if ~isnan(task_mets_out_lb(i))
        model.var_lb(find(ismember(model.varNames,strcat('F_EX_',task_mets_out(i))))) = task_mets_out_lb(i);
    end
    
    if ~isnan(task_mets_out_ub(i))
        model.var_ub(find(ismember(model.varNames,strcat('F_EX_',task_mets_out(i))))) = task_mets_out_ub(i);
    end
end


for i = 1:length(task_rxns)
    if ~isnan(task_rxns_lb(i))      
        model.rhs(find(ismember(model.varNames,strcat('NF_task_rxn',num2str(i))))) = task_rxns_lb(i);
    end
end

model.f(find(model.f)) = 0;% we don't need to maximize an objective, we only need to know if it is feasible

sol = solveTFAmodelCplex(model,30);


if (~isempty(sol.x) && ~task_should_fail) || (isempty(sol.x) && task_should_fail)
    task_success = 1;
else
    task_success = 0;
end

[strcat(task_name,':') num2cell(task_success)]

end