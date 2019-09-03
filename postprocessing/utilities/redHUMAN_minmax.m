function minmax = redHUMAN_minmax(model,thermo_data,data_Leukemia,rxns)

    if nargin<4 %then do it for intracellular rxns
        tpts_rxns = model.rxns(find(model.isTrans));
        drains = extract_drains(model);
        rxns = setdiff(model.rxns,[tpts_rxns;drains]);
    end

    % Build the specific model for the chosen cancer types
    if contains(model.description,'recon3') || contains(model.description,'Recon3')
        data_Leukemia(:,1) = strrep(data_Leukemia(:,1),'glc_e','glc_D_e');
    end
    [~,ii] = ismember(data_Leukemia(1:27,1),model.rxns);
    model.lb(ii) = cell2mat(data_Leukemia(1:27,2));
    model.ub(ii) = cell2mat(data_Leukemia(1:27,3));
    
    model = TFAprepconv(model,thermo_data);
    
    [~,jj] = ismember(data_Leukemia(28:end,1),model.varNames);
    model.var_lb(jj(find(jj))) = cell2mat(data_Leukemia(27+find(jj),2));
    model.var_ub(jj(find(jj))) = cell2mat(data_Leukemia(27+find(jj),3));
    
    model.data = [data_Leukemia(1:27,:);data_Leukemia(27+find(jj),:)];
    
    % Compute minmax
    sol = solveTFAmodelCplex(model,60);
    model.var_lb(find(model.f)) = 0.99*sol.val;
    minmax = runTMinMax(model,strcat('NF_',rxns),60);
    
    minmax(minmax > 100)  = 100;
    minmax(minmax < -100) = -100;
    
    minmax = [rxns num2cell(minmax)];
end