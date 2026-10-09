%% runPopulationCalibration.m
% Population-based calibration for IL-6 / toxicity ODE model

clc; clear; close all;

%% ------------------------------
% Load your parameter template

theta = IL6DosingParameters();  % your existing function

theta.useIDosing.val = true;       % ICI active for all patients
theta.nDoses_I.val = 4;            % total doses
theta.doseInterval_I.val = 24*7*3; % 3 weeks
theta.doseI.val = 1;               % dose size

% Optional: turn off Anti-IL6 for calibration
theta.useADosing.val = false;

%% ------------------------------
% Define which parameters to calibrate (fit globally)

%fitParams = {'rx','rxres','etox'};  % global calibration

% Parameters to vary across virtual population
popParams = {'E_initial','L_initial', 'C_initial'};  % heterogeneity

% Initial guess for optimizer (fitParams)
p0 = [0.1 9e-4];

% Lower / upper bounds
%lb = [4e-6 5e-4];
%ub = [1 1.4e-3];

% Define parameters to fit
fitParams = {'rx','rxres','etox', 'ltox'};

nVars = length(fitParams);

% Lower / upper bounds
lb = zeros(nVars,1);
ub = zeros(nVars,1);

for i = 1:nVars
    lb(i) = theta.(fitParams{i}).min;
    ub(i) = theta.(fitParams{i}).max;
end

%% ------------------------------
% Run the optimizer

options = optimoptions('particleswarm','Display','iter','UseParallel',true);

p_opt = particleswarm(@(p) populationLoss(p,theta,fitParams,popParams), ...
                      length(fitParams), lb, ub, options);

disp('Optimized parameters:');
for j = 1:length(fitParams)
    fprintf('%s = %.5e\n', fitParams{j}, p_opt(j));
    theta.(fitParams{j}).val = p_opt(j);
end

%% ------------------------------
% Analyze optimized population statistics

NP_analysis = 100;

[stats, Tonset_all, Tres_all] = populationSimulation(theta,fitParams,popParams,p_opt,NP_analysis);

% Convert to days
Tonset_days = Tonset_all / 24;
Tres_days   = Tres_all / 24;

% Remove invalid values
Tonset_days = Tonset_days(isfinite(Tonset_days));
Tres_days   = Tres_days(isfinite(Tres_days));

% Display summary
fprintf('\n--- Population Statistics (Optimized) ---\n');
fprintf('Median Tonset (days): %.2f\n', median(Tonset_days));
fprintf('Mean   Tonset (days): %.2f\n', mean(Tonset_days));
fprintf('Median Tres   (days): %.2f\n', median(Tres_days));
fprintf('Mean   Tres   (days): %.2f\n', mean(Tres_days));

figure;

subplot(1,2,1)
histogram(Tonset_days,20)
xlabel('Tonset (days)')
ylabel('Count')
title('Distribution of Toxicity Onset')

subplot(1,2,2)
histogram(Tres_days,20)
xlabel('Tres (days)')
ylabel('Count')
title('Distribution of Toxicity Resolution')
%% ------------------------------
% Validate population and plot normalized X(t)

NPval = 100;  % small population for visualization

figure; hold on
for i = 1:NPval
    theta_i = theta;

    % Sample heterogeneous parameters for this virtual patient
    for k = 1:length(popParams)
        low  = theta_i.(popParams{k}).min;
        high = theta_i.(popParams{k}).max;
        theta_i.(popParams{k}).val = low + rand*(high-low);
    end

    % Run simulation
    [t,x] = IL6RunSim(IL6ParamUnpack(theta_i));

    % Extract toxicity X(t) and normalize
    Xnorm = x(:,6) / max(x(:,6));

    % Plot normalized toxicity
    plot(t/24, Xnorm, 'LineWidth', 1.5);
end

xlabel('Time (days)');
ylabel('Normalized toxicity X(t)');
title('Population-based toxicity trajectories (normalized)');
grid on;
hold off;

%% ------------------------------
% Save calibrated parameters

save('theta_calibrated.mat','theta','fitParams','popParams');