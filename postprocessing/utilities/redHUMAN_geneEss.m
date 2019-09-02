function [EssGenes_leukemia, model] = redHUMAN_geneEss(model,thermo_data,data_Leukemia)


    if isfield(model,'info_LMDP')
        model = removeLMPDrxns(model);
    end

    model = genes_redGEM(model);
    
    %% 2. Build the specific model for the chosen cancer types
    if contains(model.description,'recon3')||contains(model.description,'Recon3')
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
    
    
    % gene essentiality
    essThr = 0;

    indF = getAllVar(model,{'F'});
    indR = getAllVar(model,{'R'});
    model.var_lb(indR) = 0;
    model.var_lb(indF) = 0;
    EssGenes_leukemia = runThermoGeneEssentiality(model,essThr);

end