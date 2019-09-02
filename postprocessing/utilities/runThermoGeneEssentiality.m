function [essential_genes_tfa] = runThermoGeneEssentiality(tmodel,essThr)

indNF = getAllVar(tmodel,{'NF'});
tmodel.var_lb(find(tmodel.f))=1e-6;
solt = solveTFAmodelCplex(tmodel,30);

[grRatio_genetfa,grRateKO_genetfa] = thermoSingleGeneDeletion(tmodel, 'TFA', tmodel.genes, 0, 0, 0, essThr, indNF);

if any(isnan(grRatio_genetfa)) %keep track of the NaN KO by comparing essential_genes_tfaNaN with essential_genes_tfa
    grRateKO_genetfaNaN = grRateKO_genetfa;
    essential_genes_tfaNaN = tmodel.genes(grRateKO_genetfaNaN(:,1)<=essThr*solt.val);
    yesgenetfaNaN = ismember(tmodel.genes,essential_genes_tfaNaN);
end

grRateKO_genetfa(isnan(grRateKO_genetfa)) = 0; %by default NaN is considered an essential KO
essential_genes_tfa = tmodel.genes(grRateKO_genetfa(:,1)<=essThr*solt.val);
yesgenetfa = ismember(tmodel.genes,essential_genes_tfa);

end
