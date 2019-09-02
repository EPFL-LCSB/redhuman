function model = removeLMPDrxns(model)

    gem = model.OriginalGEM;
    lmpd_rxns = model.rxns(contains(model.rxns,'LMPD_'));
    model = removeRxns(model,lmpd_rxns);
    
    for i = 1:size(model.info_LMDP,1)
        [~,ind_rxns] = ismember(model.info_LMDP{i,4},gem.rxns);
        for j = 1:length(ind_rxns)
            rxnName = model.info_LMDP{i,4}{j};
            metaboliteList = gem.mets(find(gem.S(:,ind_rxns(j))));
            stoichCoeffList = gem.S(find(gem.S(:,ind_rxns(j))),ind_rxns(j));
            revFlag = true;
            lowerBound = gem.lb(find(ismember(gem.rxns,rxnName)));
            upperBound = gem.ub(find(ismember(gem.rxns,rxnName)));
            objCoeff = 0;
            subSystem = gem.subSystems{find(ismember(gem.rxns,rxnName))};
            grRule = gem.grRules(find(ismember(gem.rxns,rxnName)));
            
            [model,~] = addReactionFromLMPD(model,rxnName,metaboliteList,stoichCoeffList,revFlag,lowerBound,upperBound,objCoeff,subSystem,grRule);
        end
    end
    
    % fix metabolite fields
    [~,ind_mets] = ismember(model.mets,gem.mets);
    model.metSEEDID = gem.metSEEDID(ind_mets);
    model.metCompSymbol = gem.metCompSymbol(ind_mets);
    model.metNames = gem.metNames(ind_mets);
    model.metFormulas = gem.metFormulas(ind_mets);
    
    % fix gene information
    model = genes_redGEM(model);


end


