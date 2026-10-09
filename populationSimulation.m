%% populationSimulation.m
% Input:
%   theta_template - parameter struct
%   fitParams      - cell array of parameters being calibrated
%   popParams      - cell array of parameters varied across patients
%   p              - current optimizer values for fitParams
%   NP             - number of virtual patients
% Output:
%   stats - struct with median/variance of features

function [stats, Tonset_all, Tres_all] = populationSimulation(theta_template,fitParams,popParams,p,NP)

theta = theta_template;

% insert calibration parameters
for j = 1:length(fitParams)
    theta.(fitParams{j}).val = p(j);
end

Tonset = zeros(NP,1);
Tres   = zeros(NP,1);

parfor i = 1:NP
    theta_i = theta;
    % sample heterogeneity
    for k = 1:length(popParams)
        low  = theta.(popParams{k}).min;
        high = theta.(popParams{k}).max;
        theta_i.(popParams{k}).val = low + rand*(high-low);
    end
    
    params = IL6ParamUnpack(theta_i);

    try
        [t,x] = IL6RunSim(params);
        feats = toxicityFeatures(t,x);
        Tonset(i) = feats.Tonset;
        Tres(i)   = feats.Tres;
        
        % Optional: replace NaN with end of simulation
        if isnan(Tonset(i)); Tonset(i) = params.tEnd; end
        if isnan(Tres(i));   Tres(i)   = params.tEnd; end

    catch
        % if ODE fails, assign large penalty
        Tonset(i) = params.tEnd;
        Tres(i)   = params.tEnd;
    end
end

% population statistics
stats.medTonset = median(Tonset);
stats.medTres   = median(Tres);

stats.varTonset = var(Tonset);
stats.varTres   = var(Tres);

Tonset_all = Tonset;
Tres_all   = Tres;

end