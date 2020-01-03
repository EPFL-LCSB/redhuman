function table = table_gene_rxns_ess(model,ess_genes,ess_genes_GEM)

table=[];
table2=[];
% table_rxns=[];
for i = 1:length(ess_genes)
    [~,ii] = ismember(ess_genes(i),model.genes);
    rxns_this_gene = model.rxns(find(model.rxnGeneMat(:,ii)));
    GPR_this_gene = model.grRules(find(model.rxnGeneMat(:,ii)));
%     table_rxns = [table_rxns;rxns_this_gene];
    table = [table; ess_genes(i) {''} {''} num2cell(ismember(ess_genes(i),ess_genes_GEM));...
        repmat({''},length(rxns_this_gene),1) rxns_this_gene GPR_this_gene repmat({''},length(rxns_this_gene),1)];
    
    concat_rxns_this_gene=rxns_this_gene{1};
    for j=2:length(rxns_this_gene)
        concat_rxns_this_gene = strcat(concat_rxns_this_gene,{' | '},rxns_this_gene{j});
    end
    table2 = [table2; ess_genes(i) concat_rxns_this_gene num2cell(ismember(ess_genes(i),ess_genes_GEM))];
end
