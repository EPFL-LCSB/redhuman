function [drains,drain_mets]=extract_drains(model,drains_medium)

if nargin<2
    drains_medium={};
end

drain_mets={};
drains={};
for i=1:length(model.rxns)
    no_of_mets=length(find(model.S(:, i)));
    if no_of_mets==1
        drains=[drains; model.rxns(i)];
        sto=model.S(:, i);
        drain_mets=[drain_mets; model.mets(find(sto))];
        sto=sto(find(sto));
        if sto==1
            model.S(:, i)=-model.S(:, i);
            ub=model.ub(i);
            lb=model.lb(i);
            model.ub(i)=-lb;
            model.lb(i)=-ub;
        end
    end
end

if ~isempty(drains_medium)
        drains_aux=drains;
        drains={};
        for i=1:length(drains_aux)
            if ismember(drains_aux(i),drains_medium)
                drains=[drains;drains_aux(i)];
            end
        end
end
end