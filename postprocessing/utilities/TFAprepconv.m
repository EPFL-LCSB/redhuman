function [tmodel,solTFA] = TFAprepconv(model,thermo_data,minmaxTightening,addNetFluxes)

    if nargin<3
        addNetFluxes = 1;
    end

    if nargin<3
        minmaxTightening = 0;
    end

    if minmaxTightening==1
        minmax = runMinMax(model);
        model.lb = minmax(:,1);
        model.ub = minmax(:,2);
    end

    %convert to TFA model
    tmodel = prepModelforTFA(model,thermo_data,model.CompartmentData);
    tmodel = convToTFA(tmodel,thermo_data);

    if addNetFluxes
        tmodel = addNetFluxVariables(tmodel);
    end

end