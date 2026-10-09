%% populationLoss.m
% Input:
%   p         - current guess for fitParams
%   theta     - parameter template
%   fitParams - cell array of parameters being calibrated
%   popParams - cell array of parameters varied in population
% Output:
%   J - scalar loss

function J = populationLoss(p,theta,fitParams,popParams)

NP = 100; % number of virtual patients

[stats, Tonset_all, Tres_all] = populationSimulation(theta,fitParams,popParams,p,NP);

% Check for failed simulations
if ~isfinite(stats.medTonset) || ~isfinite(stats.medTres)
    J = 1e6;
    return
end

% literature targets
Tonset_target = 60*24;  % 60 days
Tres_target   = 15*7*24; % 8 weeks

% compute loss (relative)
J = ((stats.medTonset - Tonset_target)/Tonset_target)^2 ...
  + ((stats.medTres   - Tres_target)/Tres_target)^2 ...
  + 0.1*stats.varTonset ...
  + 0.1*stats.varTres;

end