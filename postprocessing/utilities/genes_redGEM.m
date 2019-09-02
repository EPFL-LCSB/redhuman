% OBTAIN GENES FROM ORIGINAL MODEL
function model = genes_redGEM(model)
    gem = model.OriginalGEM;

    model.lb(find(model.c)) = 0;
    model.var_lb(find(model.f)+1) = 0;

    [~,ind_rxns] = ismember(model.rxns, gem.rxns);
    ind_rxns = ind_rxns(find(ind_rxns)); % because of the lumps
    rxnGeneMatrix = gem.rxnGeneMat(ind_rxns,:);

    ind_genes = find(sum(rxnGeneMatrix,1));

    model.rxnGeneMat = rxnGeneMatrix(:,ind_genes);
    model.genes = gem.genes(ind_genes);
    model.rules = gem.rules(ind_rxns);
    model.grRules = gem.grRules(ind_rxns);

    rules = model.rules;

    for i = 1:length(ind_genes)
        rules = strrep(rules,['(',num2str(ind_genes(i)),')'],['(_',num2str(ind_genes(i)),'_)']);
    end

    for i = 1:length(ind_genes)
       rules = strrep(rules,['(_',num2str(ind_genes(i)),'_)'],['(',num2str(i),')']);
    end

    model.rules = rules;

    % include the lumps (zeros everywhere)
    ind_lumps = find(contains(model.rxns,'LMPD_'));
    model.rules(ind_lumps) = {''};
    model.grRules(ind_lumps) = {''};
    model.rxnGeneMat(ind_lumps,:) = zeros(length(ind_lumps),length(model.genes));

    %singleGeneDeletion(model);
end