%% Supplementary Figure X 
% ICI_tumor_control_heatmaps
% Parameter grids
C0_vals = logspace(-2, 0, 25);
E0_vals = logspace(-2, 0, 25);

nC = numel(C0_vals);
nE = numel(E0_vals);

% Storage
C_end_OFF = zeros(nE, nC);
C_end_ON  = zeros(nE, nC);

% Sweep loop
for iC = 1:nC
    for iE = 1:nE

        % Base parameters
        theta = IL6DosingParameters();
        theta.C_initial.val = C0_vals(iC);
        theta.E_initial.val = E0_vals(iE);

        % ICI OFF
        theta_OFF = theta;
        theta_OFF.useIDosing.val = false;

        params_OFF = IL6ParamUnpack(theta_OFF);
        [t_OFF, x_OFF] = IL6RunSim(params_OFF);

        C_end_OFF(iE, iC) = x_OFF(end,1);

        % ICI ON
        theta_ON = theta;
        theta_ON.useIDosing.val = true;

        params_ON = IL6ParamUnpack(theta_ON);
        [t_ON, x_ON] = IL6RunSim(params_ON);

        C_end_ON(iE, iC) = x_ON(end,1);

    end
end

% Convert tumor burden to log10 scale
C_end_OFF_plot = log10(C_end_OFF);
C_end_ON_plot  = log10(C_end_ON);

% Shared color limits
cmin = min([C_end_OFF_plot(:); C_end_ON_plot(:)]);
cmax = max([C_end_OFF_plot(:); C_end_ON_plot(:)]);

% Figure
fig = figure('Position',[100 100 1000 550], ...
             'Color','w');

t = tiledlayout(1,2, ...
    'TileSpacing','compact', ...
    'Padding','compact');

% Common formatting
fontSize = 12;
labelSize = 14;
titleSize = 15;

% ICI OFF
ax1 = nexttile;

imagesc(C0_vals, E0_vals, C_end_OFF_plot);
set(gca,'YDir','normal');

caxis([cmin cmax]);
colormap(ax1,turbo);

set(gca, ...
    'XScale','linear', ...
    'YScale','linear', ...
    'FontSize',fontSize, ...
    'FontName','Times New Roman', ...
    'LineWidth',1, ...
    'TickDir','out', ...
    'Box','off', ...
    'Layer','top');

xlabel('$\mathrm{Initial\ Cancer\ Density},\ C_{\mathrm{initial}}$', ...
    'Interpreter','latex', ...
    'FontName','Times New Roman', ...
    'FontSize',labelSize);

ylabel('$\mathrm{Initial\ Effector\ Density},\ E_{\mathrm{initial}}$', ...
    'Interpreter','latex', ...
    'FontName','Times New Roman', ...
    'FontSize',labelSize);

title('Without ICI', ...
    'FontSize',titleSize, ...
    'FontWeight','normal');

% ICI ON
ax2 = nexttile;

imagesc(C0_vals, E0_vals, C_end_ON_plot);
set(gca,'YDir','normal');

caxis([cmin cmax]);
colormap(ax2,turbo);

set(gca, ...
    'XScale','linear', ...
    'YScale','linear', ...
    'FontSize',fontSize, ...
    'FontName','Times New Roman', ...
    'LineWidth',1, ...
    'TickDir','out', ...
    'Box','off', ...
    'Layer','top');

xlabel('$\mathrm{Initial\ Cancer\ Density},\ C_{\mathrm{initial}}$', ...
    'Interpreter','latex', ...
    'FontName','Times New Roman', ...
    'FontSize',labelSize);

ylabel('$\mathrm{Initial\ Effector\ Density},\ E_{\mathrm{initial}}$', ...
    'Interpreter','latex', ...
    'FontName','Times New Roman', ...
    'FontSize',labelSize);

title('With ICI', ...
    'FontSize',titleSize, ...
    'FontWeight','normal');

% Shared colorbar underneath both heatmaps
cb = colorbar(ax2,'southoutside');

cb.Layout.Tile = 'south';

cb.Label.String = '$\log_{10}$(Tumor Density)';
cb.Label.Interpreter = 'latex';
cb.Label.FontSize = labelSize;

% Numerical colorbar labels
cb.Ticks = linspace(cmin,cmax,5);
cb.TickLabels = arrayfun(@(x)sprintf('%.2f',x), ...
                          cb.Ticks, ...
                          'UniformOutput',false);

% Ensure identical axis limits
xlim(ax1,[min(C0_vals) max(C0_vals)]);
xlim(ax2,[min(C0_vals) max(C0_vals)]);

ylim(ax1,[min(E0_vals) max(E0_vals)]);
ylim(ax2,[min(E0_vals) max(E0_vals)]);

exportgraphics(fig,'ICI_tumor_control_heatmaps.png','Resolution',600);
exportgraphics(fig,'ICI_tumor_control_heatmaps.pdf','ContentType','vector');

%% Figure 6
% Tumor Outcome vs Anti-IL-6 timing (2x3)

% Tumor outcome and toxicity with/without anti-IL-6 administration
% Across three initial tumor burdens

clear;
clc;

% ==============================
% Figure parameters
% ==============================

fontSize   = 11;
labelSize  = 12;
titleSize  = 13;
headerSize = 13;

lineWidth  = 2.0;
markerSize = 5;

% ==============================
% Base parameters
% ==============================

paramsStruct = IL6DosingParameters();

% Convert struct-of-structs -> flat struct of values
params = structfun(@(s) s.val, paramsStruct, ...
                   'UniformOutput', false);

% Enable anti-CTLA-4 dosing
params.useIDosing = true;

% Sweep: time of anti-IL-6 administration
doseTimes = linspace(2,3999,30);
nDose = length(doseTimes);

% Population sweep size
nPop = 500;

% Initial tumor conditions
C_initial_vals = [0.05, 0.20, 0.50];

nRuns = length(C_initial_vals);

% ==============================
% Storage
% ==============================

tumor_with     = cell(nRuns,1);
tumor_without  = cell(nRuns,1);

tox_with       = cell(nRuns,1);
tox_without    = cell(nRuns,1);

% ==============================
% Main simulation sweep
% ==============================

for k = 1:nRuns
    
    fprintf('\nRunning C_initial = %.2f\n', ...
            C_initial_vals(k));
    
    % Set initial tumor size
    params.C_initial = C_initial_vals(k);
    
    % Storage for this run
    tumor_with{k}    = zeros(nDose,nPop);
    tumor_without{k} = zeros(nDose,nPop);
    
    tox_with{k}      = zeros(nDose,nPop);
    tox_without{k}   = zeros(nDose,nPop);
    
    for i = 1:nDose
        
        fprintf('  Dose time %.1f / %.1f hours\n', ...
                doseTimes(i),doseTimes(end));
        
        params.doseStart_A = doseTimes(i);
        
        for j = 1:nPop
            
            % ------------------------------
            % Sample biological variability
            % ------------------------------
            
            params.E_initial = ...
                paramsStruct.E_initial.min + ...
                (paramsStruct.E_initial.max - ...
                 paramsStruct.E_initial.min) * rand();
            
            params.L_initial = ...
                paramsStruct.L_initial.min + ...
                (paramsStruct.L_initial.max - ...
                 paramsStruct.L_initial.min) * rand();
            
            % ==============================
            % WITH anti-IL-6
            % ==============================
            
            params.useADosing = true;
            
            [t,x] = IL6RunSim(params);
            
            tumor_with{k}(i,j) = x(end,1);
            tox_with{k}(i,j)   = max(x(:,6));
            
            % ==============================
            % WITHOUT anti-IL-6
            % ==============================
            
            params.useADosing = false;
            
            [t,x] = IL6RunSim(params);
            
            tumor_without{k}(i,j) = x(end,1);
            tox_without{k}(i,j)   = max(x(:,6));
            
        end
    end
end

% ==============================
% Normalize toxicity
% ==============================

globalMax = max([ ...
    tox_with{1}(:); ...
    tox_with{2}(:); ...
    tox_with{3}(:); ...
    tox_without{1}(:); ...
    tox_without{2}(:); ...
    tox_without{3}(:)]);

for k = 1:nRuns
    
    tox_with{k}    = tox_with{k} / globalMax;
    tox_without{k} = tox_without{k} / globalMax;
    
end

% ==============================
% Summary statistics
% ==============================

med_tw  = cell(nRuns,1);
min_tw  = cell(nRuns,1);
max_tw  = cell(nRuns,1);

med_two = cell(nRuns,1);
min_two = cell(nRuns,1);
max_two = cell(nRuns,1);

med_xw  = cell(nRuns,1);
min_xw  = cell(nRuns,1);
max_xw  = cell(nRuns,1);

med_xwo = cell(nRuns,1);
min_xwo = cell(nRuns,1);
max_xwo = cell(nRuns,1);

for k = 1:nRuns
    
    % Tumor
    med_tw{k} = median(tumor_with{k},2);
    min_tw{k} = min(tumor_with{k},[],2);
    max_tw{k} = max(tumor_with{k},[],2);
    
    med_two{k} = median(tumor_without{k},2);
    min_two{k} = min(tumor_without{k},[],2);
    max_two{k} = max(tumor_without{k},[],2);
    
    % Toxicity
    med_xw{k} = median(tox_with{k},2);
    min_xw{k} = min(tox_with{k},[],2);
    max_xw{k} = max(tox_with{k},[],2);
    
    med_xwo{k} = median(tox_without{k},2);
    min_xwo{k} = min(tox_without{k},[],2);
    max_xwo{k} = max(tox_without{k},[],2);
    
end

% ==============================
% Colors
% ==============================

% With anti-IL-6
colorWith = [0.85 0.20 0.20];

% Without anti-IL-6
colorWithout = [0.20 0.20 0.85];

% ==============================
% Create figure
% ==============================

fig = figure( ...
    'Name','Anti-IL-6 Dosing Timing', ...
    'NumberTitle','off', ...
    'Color','w', ...
    'Position',[100 100 1250 700]);



% ==============================
% Plot all six panels
% ==============================

ax = gobjects(2,3);

for k = 1:nRuns
    
    % ==========================================================
    % TOP ROW — TUMOR OUTCOME
    % ==========================================================
    
    ax(1,k) = subplot(2,3,k);
    hold on;
    
    % With anti-IL-6 range
    fill([doseTimes fliplr(doseTimes)], ...
         [min_tw{k}' fliplr(max_tw{k}')], ...
         colorWith, ...
         'FaceAlpha',0.15, ...
         'EdgeColor','none');
    
    % Without anti-IL-6 range
    fill([doseTimes fliplr(doseTimes)], ...
         [min_two{k}' fliplr(max_two{k}')], ...
         colorWithout, ...
         'FaceAlpha',0.15, ...
         'EdgeColor','none');
    
    % Median lines
    plot(doseTimes,med_tw{k}, ...
         'Color',colorWith, ...
         'LineWidth',lineWidth);
    
    plot(doseTimes,med_two{k}, ...
         'Color',colorWithout, ...
         'LineWidth',lineWidth);
    
    % Formatting
    xlabel('Anti-IL-6 administration time (hours)', ...
        'FontName','Times New Roman', ...
        'FontSize',labelSize);
    
    ylabel('Final tumor burden', ...
        'FontName','Times New Roman', ...
        'FontSize',labelSize);
    
    title('Tumor Burden', ...
        'FontName','Times New Roman', ...
        'FontSize',titleSize, ...
        'FontWeight','normal');
    
    set(gca, ...
        'FontName','Times New Roman', ...
        'FontSize',fontSize, ...
        'LineWidth',1, ...
        'TickDir','out', ...
        'Box','off', ...
        'Layer','top');
    
    grid on;
    ax(1,k).GridAlpha = 0.15;
    ax(1,k).MinorGridAlpha = 0.08;
    
    xlim([min(doseTimes) max(doseTimes)]);
    
    
    % ==========================================================
    % BOTTOM ROW — TOXICITY
    % ==========================================================
    
    ax(2,k) = subplot(2,3,k+3);
    hold on;
    
    % With anti-IL-6 range
    fill([doseTimes fliplr(doseTimes)], ...
         [min_xw{k}' fliplr(max_xw{k}')], ...
         colorWith, ...
         'FaceAlpha',0.15, ...
         'EdgeColor','none');
    
    % Without anti-IL-6 range
    fill([doseTimes fliplr(doseTimes)], ...
         [min_xwo{k}' fliplr(max_xwo{k}')], ...
         colorWithout, ...
         'FaceAlpha',0.15, ...
         'EdgeColor','none');
    
    % Median lines
    plot(doseTimes,med_xw{k}, ...
         'Color',colorWith, ...
         'LineWidth',lineWidth);
    
    plot(doseTimes,med_xwo{k}, ...
         'Color',colorWithout, ...
         'LineWidth',lineWidth);
    
    % Formatting
    xlabel('Anti-IL-6 administration time (hours)', ...
        'FontName','Times New Roman', ...
        'FontSize',labelSize);
    
    ylabel('Normalized peak toxicity', ...
        'FontName','Times New Roman', ...
        'FontSize',labelSize);
    
    title('Normalized Toxicity', ...
        'FontName','Times New Roman', ...
        'FontSize',titleSize, ...
        'FontWeight','normal');
    
    set(gca, ...
        'FontName','Times New Roman', ...
        'FontSize',fontSize, ...
        'LineWidth',1, ...
        'TickDir','out', ...
        'Box','off', ...
        'Layer','top');
    
    grid on;
    ax(2,k).GridAlpha = 0.15;
    ax(2,k).MinorGridAlpha = 0.08;
    
    xlim([min(doseTimes) max(doseTimes)]);
    ylim([0 1]);
    
end

% ==============================
% Column headers
% ==============================

annotation('textbox',[0.145 0.95 0.20 0.035], ...
    'String','(A)  $C_{\mathrm{initial}} = 0.05$', ...
    'Interpreter','latex', ...
    'FontName','Times New Roman', ...
    'FontSize',headerSize, ...
    'FontWeight','bold', ...
    'HorizontalAlignment','center', ...
    'VerticalAlignment','middle', ...
    'EdgeColor','none');

annotation('textbox',[0.400 0.95 0.20 0.035], ...
    'String','(B)  $C_{\mathrm{initial}} = 0.20$', ...
    'Interpreter','latex', ...
    'FontName','Times New Roman', ...
    'FontSize',headerSize, ...
    'FontWeight','bold', ...
    'HorizontalAlignment','center', ...
    'VerticalAlignment','middle', ...
    'EdgeColor','none');

annotation('textbox',[0.655 0.95 0.20 0.035], ...
    'String','(C)  $C_{\mathrm{initial}} = 0.50$', ...
    'Interpreter','latex', ...
    'FontName','Times New Roman', ...
    'FontSize',headerSize, ...
    'FontWeight','bold', ...
    'HorizontalAlignment','center', ...
    'VerticalAlignment','middle', ...
    'EdgeColor','none');

% ==============================
% Custom legend — colored boxes
% ==============================

% Remove standard legend if present
legend(ax(1,1),'off');

% ---- With anti-IL-6 ----
annotation('rectangle',[0.54 0.02 0.018 0.018], ...
    'FaceColor',colorWith, ...
    'EdgeColor','none');

annotation('textbox',[0.56 0.02 0.18 0.018], ...
    'String','With anti-IL-6', ...
    'FontName','Times New Roman', ...
    'FontSize',11, ...
    'HorizontalAlignment','left', ...
    'VerticalAlignment','middle', ...
    'EdgeColor','none');

% ---- Without anti-IL-6 ----
annotation('rectangle',[0.4 0.02 0.018 0.018], ...
    'FaceColor',colorWithout, ...
    'EdgeColor','none');

annotation('textbox',[0.42 0.02 0.180 0.018], ...
    'String','Without anti-IL-6', ...
    'FontName','Times New Roman', ...
    'FontSize',11, ...
    'HorizontalAlignment','left', ...
    'VerticalAlignment','middle', ...
    'EdgeColor','none');
% ==============================
% Export
% ==============================

exportgraphics(fig, ...
    'AntiIL6_DosingTiming_Comparison.pdf', ...
    'ContentType','vector');

exportgraphics(fig, ...
    'AntiIL6_DosingTiming_Comparison.png', ...
    'Resolution',600);

%% Figure 2
% 2x3 Heatmap Figure: Tumor Burden & Toxicity vs E_initial & L_initial
% For three initial tumor sizes
% ==============================

% Load parameters
paramsStruct = IL6DosingParameters();
params = structfun(@(s) s.val, paramsStruct, 'UniformOutput', false);

% Enable therapies
params.useIDosing = true;
params.useADosing = false;

% Sweep ranges
nE = 25;
nL = 25;

E_vals = linspace(paramsStruct.E_initial.min, paramsStruct.E_initial.max, nE);
L_vals = linspace(paramsStruct.L_initial.min, paramsStruct.L_initial.max, nL);

% Define initial tumor sizes
C_init_vals = [0.05, 0.2, 0.5];

% Preallocate heatmaps
tumorHeat = zeros(nE, nL);
toxHeat   = zeros(nE, nL);

% ==============================
% Figure parameters
% ==============================

fontSize  = 11;
labelSize = 12;
titleSize = 12;

fig = figure('Name','Tumor & Toxicity Heatmaps', ...
    'NumberTitle','off', ...
    'Color','w', ...
    'Position',[200 200 1100 650]);

% ==============================
% Generate heatmaps
% ==============================

for c_idx = 1:length(C_init_vals)

    % Set initial tumor size
    params.C_initial = C_init_vals(c_idx);

    % Sweep over E_initial and L_initial
    for i = 1:nE
        for j = 1:nL

            params.E_initial = E_vals(i);
            params.L_initial = L_vals(j);

            % Run simulation
            [t,x] = IL6RunSim(params);

            % Extract final tumor and max toxicity
            tumorHeat(i,j) = x(end,1);
            toxHeat(i,j)   = max(x(:,6));
        end
    end

    % Normalize toxicity for visualization
    toxHeatNorm = toxHeat / max(toxHeat(:));

    % ==============================
    % TOP ROW — Tumor burden
    % ==============================

    ax1 = subplot(2,3,c_idx);

    imagesc(L_vals,E_vals,tumorHeat);

    set(ax1, ...
        'YDir','normal', ...
        'FontName','Times New Roman', ...
        'FontSize',fontSize, ...
        'LineWidth',1, ...
        'TickDir','out', ...
        'Box','off', ...
        'XColor','k', ...
        'YColor','k', ...
        'Layer','top');

    xlabel('$L_{\mathrm{initial}}$ (IL-6)', ...
        'Interpreter','latex', ...
        'FontName','Times New Roman', ...
        'FontSize',labelSize);

    ylabel('$E_{\mathrm{initial}}$ (Effector cells)', ...
        'Interpreter','latex', ...
        'FontName','Times New Roman', ...
        'FontSize',labelSize);

    title('Tumor Burden', ...
        'FontName','Times New Roman', ...
        'FontSize',titleSize, ...
        'Color','k', ...
        'FontWeight','normal');

    cb1 = colorbar;
    cb1.FontName = 'Times New Roman';
    cb1.FontSize = fontSize;
    cb1.Color = 'k';
    cb1.Label.Color = 'k';

    colormap(ax1,parula);


    % ==============================
    % BOTTOM ROW — Normalized toxicity
    % ==============================

    ax2 = subplot(2,3,c_idx+3);

    imagesc(L_vals,E_vals,toxHeatNorm);

    set(ax2, ...
        'YDir','normal', ...
        'FontName','Times New Roman', ...
        'FontSize',fontSize, ...
        'LineWidth',1, ...
        'TickDir','out', ...
        'Box','off', ...
        'XColor','k', ...
        'YColor','k', ...
        'Layer','top');

    xlabel('$L_{\mathrm{initial}}$ (IL-6)', ...
        'Interpreter','latex', ...
        'FontName','Times New Roman', ...
        'FontSize',labelSize);

    ylabel('$E_{\mathrm{initial}}$ (Effector cells)', ...
        'Interpreter','latex', ...
        'FontName','Times New Roman', ...
        'FontSize',labelSize);

    title('Normalized Toxicity', ...
        'FontName','Times New Roman', ...
        'FontSize',titleSize, ...
        'Color','k', ...
        'FontWeight','normal');

    cb2 = colorbar;
    cb2.FontName = 'Times New Roman';
    cb2.FontSize = fontSize;
    cb2.Color = 'k';
    cb2.Label.Color = 'k';

    colormap(ax2,hot);

end


% ============================================================
% COLUMN HEADERS — Initial Tumor Density
% ============================================================

% Retrieve axis positions
axTop1 = subplot(2,3,1);
axTop2 = subplot(2,3,2);
axTop3 = subplot(2,3,3);

axBot1 = subplot(2,3,4);
axBot2 = subplot(2,3,5);
axBot3 = subplot(2,3,6);

posTop1 = axTop1.Position;
posTop2 = axTop2.Position;
posTop3 = axTop3.Position;

posBot1 = axBot1.Position;
posBot2 = axBot2.Position;
posBot3 = axBot3.Position;


% ==============================
% Column header 1
% ==============================

annotation(fig,'textbox', ...
    [posTop1(1), posTop1(2)+posTop1(4)+0.035, posTop1(3), 0.035], ...
    'String','$C_{\mathrm{initial}} = 0.05$', ...
    'Interpreter','latex', ...
    'FontName','Times New Roman', ...
    'FontSize',13, ...
    'FontWeight','normal', ...
    'HorizontalAlignment','center', ...
    'VerticalAlignment','middle', ...
    'Color','k', ...
    'EdgeColor','none');


% ==============================
% Column header 2
% ==============================

annotation(fig,'textbox', ...
    [posTop2(1), posTop2(2)+posTop2(4)+0.035, posTop2(3), 0.035], ...
    'String','$C_{\mathrm{initial}} = 0.20$', ...
    'Interpreter','latex', ...
    'FontName','Times New Roman', ...
    'FontSize',13, ...
    'FontWeight','normal', ...
    'HorizontalAlignment','center', ...
    'VerticalAlignment','middle', ...
    'Color','k', ...
    'EdgeColor','none');


% ==============================
% Column header 3
% ==============================

annotation(fig,'textbox', ...
    [posTop3(1), posTop3(2)+posTop3(4)+0.035, posTop3(3), 0.035], ...
    'String','$C_{\mathrm{initial}} = 0.50$', ...
    'Interpreter','latex', ...
    'FontName','Times New Roman', ...
    'FontSize',13, ...
    'FontWeight','normal', ...
    'HorizontalAlignment','center', ...
    'VerticalAlignment','middle', ...
    'Color','k', ...
    'EdgeColor','none');


% ============================================================
% VISUAL GROUPING — Couple Top & Bottom Panels
% ============================================================

bracketWidth  = 0.012;
bracketOffset = 0.008;
bracketColor = [0.5 0.5 0.5];
bracketLineWidth = 1;


% ==============================
% Column 1 bracket
% ==============================

x1 = posTop1(1) + posTop1(3) + bracketOffset;
y1 = posBot1(2);
y2 = posTop1(2) + posTop1(4);

annotation(fig,'line',[x1 x1],[y1 y2], ...
    'Color',bracketColor, ...
    'LineWidth',bracketLineWidth);

annotation(fig,'line',[x1-bracketWidth x1],[y1 y1], ...
    'Color',bracketColor, ...
    'LineWidth',bracketLineWidth);

annotation(fig,'line',[x1-bracketWidth x1],[y2 y2], ...
    'Color',bracketColor, ...
    'LineWidth',bracketLineWidth);


% ==============================
% Column 2 bracket
% ==============================

x2 = posTop2(1) + posTop2(3) + bracketOffset;
y1 = posBot2(2);
y2 = posTop2(2) + posTop2(4);

annotation(fig,'line',[x2 x2],[y1 y2], ...
    'Color',bracketColor, ...
    'LineWidth',bracketLineWidth);

annotation(fig,'line',[x2-bracketWidth x2],[y1 y1], ...
    'Color',bracketColor, ...
    'LineWidth',bracketLineWidth);

annotation(fig,'line',[x2-bracketWidth x2],[y2 y2], ...
    'Color',bracketColor, ...
    'LineWidth',bracketLineWidth);


% ==============================
% Column 3 bracket
% ==============================

x3 = posTop3(1) + posTop3(3) + bracketOffset;
y1 = posBot3(2);
y2 = posTop3(2) + posTop3(4);

annotation(fig,'line',[x3 x3],[y1 y2], ...
    'Color',bracketColor, ...
    'LineWidth',bracketLineWidth);

annotation(fig,'line',[x3-bracketWidth x3],[y1 y1], ...
    'Color',bracketColor, ...
    'LineWidth',bracketLineWidth);

annotation(fig,'line',[x3-bracketWidth x3],[y2 y2], ...
    'Color',bracketColor, ...
    'LineWidth',bracketLineWidth);


% ============================================================
% EXPORT
% ============================================================

exportgraphics(fig, ...
    'Tumor_Toxicity_Heatmaps.pdf', ...
    'ContentType','vector');

exportgraphics(fig, ...
    'Tumor_Toxicity_Heatmaps.png', ...
    'Resolution',600);


%% Figure 4
% Anti-IL-6 Benefit Region
% Rows = Anti-IL-6 administration timing
% Columns = Initial tumor density
% ==============================

% Load parameters
paramsStruct = IL6DosingParameters();
params = structfun(@(s) s.val, paramsStruct, 'UniformOutput', false);

% Sweep ranges
nE = 30;
nL = 30;

E_vals = linspace(paramsStruct.E_initial.min, ...
                  paramsStruct.E_initial.max, nE);

L_vals = linspace(paramsStruct.L_initial.min, ...
                  paramsStruct.L_initial.max, nL);

% Initial tumor sizes (columns)
C_init_vals = [0.02, 0.1, 0.2];

% Anti-IL-6 administration times (rows)
doseStart_vals = [2, 500, 1000];

% ==============================
% Figure formatting
% ==============================

fontSize   = 11;
labelSize  = 12;
headerSize = 13;

fig = figure( ...
    'Name','Anti-IL-6 Benefit Regions', ...
    'NumberTitle','off', ...
    'Color','w', ...
    'Position',[100 100 1100 850]);

% Binary colormap
% 0 = no benefit
% 1 = tumor reduction

cmap = [0.82 0.82 0.82;
        0.00 0.65 0.00];

% ==============================
% Main grid dimensions
% ==============================

left   = 0.14;
bottom = 0.17;

gridWidth  = 0.75;
gridHeight = 0.75;

% Increase spacing between panels
gapX = 0.045;
gapY = 0.075;

% Individual panel dimensions
axWidth  = (gridWidth - 2*gapX)/3;
axHeight = (gridHeight - 2*gapY)/3;

% ==============================
% Generate heatmaps
% ==============================

for r = 1:3

    params.doseStart_A = doseStart_vals(r);

    for c = 1:3

        params.C_initial = C_init_vals(c);

        mask = zeros(nE,nL);

        % ------------------------------
        % Sweep E_initial and L_initial
        % ------------------------------

        for i = 1:nE
            for j = 1:nL

                params.E_initial = E_vals(i);
                params.L_initial = L_vals(j);

                % No anti-IL-6
                params.useADosing = false;
                [~,x_noA] = IL6RunSim(params);

                tumor_noA = x_noA(end,1);

                % With anti-IL-6
                params.useADosing = true;
                [~,x_withA] = IL6RunSim(params);

                tumor_withA = x_withA(end,1);

                % Benefit
                if tumor_withA < tumor_noA
                    mask(i,j) = 1;
                end

            end
        end

        % ------------------------------
        % Axis position
        % ------------------------------

        xPos = left + (c-1)*(axWidth + gapX);
        yPos = bottom + (3-r)*(axHeight + gapY);

        ax = axes( ...
            'Parent',fig, ...
            'Position',[xPos yPos axWidth axHeight]);

        % Heatmap
        imagesc(ax,L_vals,E_vals,mask);

        set(ax, ...
            'YDir','normal', ...
            'FontName','Times New Roman', ...
            'FontSize',fontSize, ...
            'LineWidth',1, ...
            'TickDir','out', ...
            'Box','off', ...
            'XColor','k', ...
            'YColor','k', ...
            'Layer','top');

        colormap(ax,cmap);
        caxis(ax,[0 1]);

        % ------------------------------
        % ONLY label shared axes
        % ------------------------------

        % X-axis labels only on bottom row
        if r == 3

            xlabel(ax, ...
                '$L_{\mathrm{initial}}$', ...
                'Interpreter','latex', ...
                'FontName','Times New Roman', ...
                'FontSize',labelSize);

        else

            ax.XTickLabel = [];

        end

        % Y-axis labels only on left column
        if c == 1

            ylabel(ax, ...
                '$E_{\mathrm{initial}}$', ...
                'Interpreter','latex', ...
                'FontName','Times New Roman', ...
                'FontSize',labelSize);

        else

            ax.YTickLabel = [];

        end

    end
end

% ==============================
% Column headers
% ==============================

columnLabels = { ...
    '$C_{\mathrm{initial}} = 0.02$', ...
    '$C_{\mathrm{initial}} = 0.10$', ...
    '$C_{\mathrm{initial}} = 0.20$'};

for c = 1:3

    xCenter = left + ...
        (c-1)*(axWidth + gapX) + ...
        axWidth/2;

    annotation(fig,'textbox', ...
        [xCenter-0.09, 0.935, 0.18, 0.035], ...
        'String',columnLabels{c}, ...
        'Interpreter','latex', ...
        'FontName','Times New Roman', ...
        'FontSize',headerSize, ...
        'FontWeight','normal', ...
        'HorizontalAlignment','center', ...
        'VerticalAlignment','middle', ...
        'Color','k', ...
        'EdgeColor','none');

end

% ==============================
% Row labels
% ==============================

rowLabels = { ...
    '2 Hours', ...
    '500 Hours', ...
    '1000 Hours'};

for r = 1:3

    yCenter = bottom + ...
        (3-r)*(axHeight + gapY) + ...
        axHeight/2;

    annotation(fig,'textbox', ...
        [0.045, yCenter-0.025, 0.04, 0.02], ...
        'String',rowLabels{r}, ...
        'FontName','Times New Roman', ...
        'FontSize',headerSize, ...
        'FontWeight','normal', ...
        'HorizontalAlignment','center', ...
        'VerticalAlignment','middle', ...
        'Color','k', ...
        'EdgeColor','none');

end

% ==============================
% Row descriptor
% ==============================

annotation(fig,'textbox', ...
    [0.015, bottom+0.17, 0.035, 0.28], ...
    'String',{'Anti-IL-6','administration','time'}, ...
    'FontName','Times New Roman', ...
    'FontSize',11, ...
    'HorizontalAlignment','center', ...
    'VerticalAlignment','middle', ...
    'Rotation',90, ...
    'Color','k', ...
    'EdgeColor','none');


% ==============================
% Legend
% ==============================

annotation(fig,'rectangle', ...
    [0.35 0.055 0.025 0.025], ...
    'FaceColor',cmap(1,:), ...
    'EdgeColor','k', ...
    'LineWidth',0.75);

annotation(fig,'textbox', ...
    [0.385 0.048 0.12 0.04], ...
    'String','No benefit', ...
    'FontName','Times New Roman', ...
    'FontSize',11, ...
    'HorizontalAlignment','left', ...
    'VerticalAlignment','middle', ...
    'Color','k', ...
    'EdgeColor','none');

annotation(fig,'rectangle', ...
    [0.53 0.055 0.025 0.025], ...
    'FaceColor',cmap(2,:), ...
    'EdgeColor','k', ...
    'LineWidth',0.75);

annotation(fig,'textbox', ...
    [0.565 0.048 0.16 0.04], ...
    'String','Tumor reduction', ...
    'FontName','Times New Roman', ...
    'FontSize',11, ...
    'HorizontalAlignment','left', ...
    'VerticalAlignment','middle', ...
    'Color','k', ...
    'EdgeColor','none');

% Export

% ==============================
drawnow;
set(fig,'PaperPositionMode','auto');

exportgraphics(fig, ...
    'Anti_IL6_Benefit_Regions.pdf', ...
    'ContentType','vector', ...
    'BackgroundColor','white');
%% Figure 5
% 3x3 Heatmaps ==============================
% 3x3 Heatmaps: Optimal Timing + Outcomes
% Optimal Anti-IL-6 Administration
% ==============================

% Load parameters
paramsStruct = IL6DosingParameters();
params = structfun(@(s) s.val, paramsStruct, 'UniformOutput', false);

% Enable dosing
params.useIDosing = true;
params.useADosing = true;

% Dose sweep
doseTimes = linspace(2,3999,40);

% Initial tumor sizes
C_init_vals = [0.02,0.1,0.2];

% Initial condition grids
E_vals = linspace(0.01,0.5,20);
L_vals = linspace(0.01,5,20);

% ==============================
% Figure formatting
% ==============================

fontSize   = 11;
labelSize  = 12;
titleSize  = 12;
headerSize = 13;

fig = figure( ...
    'Name','Optimal Anti-IL-6 Timing + Outcomes', ...
    'NumberTitle','off', ...
    'Color','w', ...
    'Units','inches', ...
    'Position',[1 1 12 9]);

% ==============================
% Colormaps
% ==============================

cmapTiming = parula(256);
cmapTumor  = parula(256);
cmapTox    = hot(256);

% ==============================
% Layout
% ==============================

left   = 0.12;
bottom = 0.13;

gridWidth  = 0.78;
gridHeight = 0.78;

gapX = 0.055;
gapY = 0.075;

axWidth  = (gridWidth - 2*gapX)/3;
axHeight = (gridHeight - 2*gapY)/3;

% ==============================
% Storage for color limits
% ==============================

allTiming = [];
allTumor  = [];
allTox    = [];

% ==============================
% Generate heatmaps
% ==============================

for c_idx = 1:length(C_init_vals)

    optimalDoseMat = nan(length(E_vals),length(L_vals));
    tumorMat       = zeros(length(E_vals),length(L_vals));
    toxMat         = zeros(length(E_vals),length(L_vals));

    for e_idx = 1:length(E_vals)

        for l_idx = 1:length(L_vals)

            % Initial conditions
            params.E_initial = E_vals(e_idx);
            params.L_initial = L_vals(l_idx);
            params.C_initial = C_init_vals(c_idx);

            % =========================
            % BASELINE: No anti-IL-6
            % =========================

            params.useADosing = false;

            [~,x_base] = IL6RunSim(params);

            tumor_base = x_base(end,1);
            tox_base   = max(x_base(:,6));

            % =========================
            % WITH anti-IL-6
            % =========================

            params.useADosing = true;

            finalTumors = zeros(length(doseTimes),1);
            peakTox     = zeros(length(doseTimes),1);

            for d_idx = 1:length(doseTimes)

                params.doseStart_A = doseTimes(d_idx);

                [~,x] = IL6RunSim(params);

                finalTumors(d_idx) = x(end,1);
                peakTox(d_idx)     = max(x(:,6));

            end

            % =========================
            % Find minimum tumor
            % =========================

            [minTumor,optIdx] = min(finalTumors);

            optDose = doseTimes(optIdx);

            % =========================
            % Validity conditions
            % =========================

            isBenefit = minTumor < tumor_base;

            % =========================
            % Store results
            % =========================

            if isBenefit

                optimalDoseMat(e_idx,l_idx) = optDose;
                tumorMat(e_idx,l_idx)       = minTumor;
                toxMat(e_idx,l_idx)         = peakTox(optIdx);
            
            else
            
                tumorMat(e_idx,l_idx) = tumor_base;
                toxMat(e_idx,l_idx)   = tox_base;
            
            end

        end
    end

    % ==============================
    % Normalize toxicity
    % ==============================

    toxMat = toxMat ./ max(toxMat(:));

    % Store valid values for shared limits
    allTiming = [allTiming; optimalDoseMat(~isnan(optimalDoseMat))];
    allTumor  = [allTumor; tumorMat(:)];
    allTox    = [allTox; toxMat(:)];

    % ==============================
    % Plot positions
    % ==============================

    % Column = initial tumor size
    % Row 1 = optimal timing
    % Row 2 = final tumor
    % Row 3 = toxicity

    xPos = left + ...
        (c_idx-1)*(axWidth + gapX);

    % ==============================
    % Row 1 — Optimal timing
    % ==============================

    yPos = bottom + 2*(axHeight + gapY);

    ax1 = axes( ...
        'Parent',fig, ...
        'Position',[xPos yPos axWidth axHeight]);

    imagesc(ax1,L_vals,E_vals,optimalDoseMat);

    set(ax1, ...
        'YDir','normal', ...
        'Color',[0.70 0.70 0.70], ...
        'FontName','Times New Roman', ...
        'FontSize',fontSize, ...
        'LineWidth',1, ...
        'TickDir','out', ...
        'Box','off', ...
        'XColor','k', ...
        'YColor','k', ...
        'Layer','top');

    h = get(ax1,'Children');
    set(h,'AlphaData',~isnan(optimalDoseMat));

    colormap(ax1,cmapTiming);
    caxis(ax1,[min(doseTimes) max(doseTimes)]);

  

    if c_idx == 1

        ylabel(ax1, ...
            '$E_{\mathrm{initial}}$', ...
            'Interpreter','latex', ...
            'FontName','Times New Roman', ...
            'Color','k', ...
            'FontSize',labelSize);

    else

        ax1.YTickLabel = [];

    end

    % ==============================
    % Row 2 — Final tumor
    % ==============================

    yPos = bottom + axHeight + gapY;

    ax2 = axes( ...
        'Parent',fig, ...
        'Position',[xPos yPos axWidth axHeight]);

    imagesc(ax2,L_vals,E_vals,tumorMat);

    set(ax2, ...
        'YDir','normal', ...
        'FontName','Times New Roman', ...
        'FontSize',fontSize, ...
        'LineWidth',1, ...
        'TickDir','out', ...
        'Box','off', ...
        'XColor','k', ...
        'YColor','k', ...
        'Layer','top');

    colormap(ax2,cmapTumor);

   

    if c_idx == 1

        ylabel(ax2, ...
            '$E_{\mathrm{initial}}$', ...
            'Interpreter','latex', ...
            'FontName','Times New Roman', ...
            'Color','k', ...
            'FontSize',labelSize);

    else

        ax2.YTickLabel = [];

    end

    % ==============================
    % Row 3 — Relative toxicity
    % ==============================

    yPos = bottom;

    ax3 = axes( ...
        'Parent',fig, ...
        'Position',[xPos yPos axWidth axHeight]);

    imagesc(ax3,L_vals,E_vals,toxMat);

    set(ax3, ...
        'YDir','normal', ...
        'FontName','Times New Roman', ...
        'FontSize',fontSize, ...
        'LineWidth',1, ...
        'TickDir','out', ...
        'Box','off', ...
        'XColor','k', ...
        'YColor','k', ...
        'Layer','top');

    colormap(ax3,cmapTox);

    

    if c_idx == 1

        ylabel(ax3, ...
            '$E_{\mathrm{initial}}$', ...
            'Interpreter','latex', ...
            'FontName','Times New Roman', ...
            'Color','k', ...
            'FontSize',labelSize);

    else

        ax3.YTickLabel = [];

    end

    xlabel(ax3, ...
        '$L_{\mathrm{initial}}$', ...
        'Interpreter','latex', ...
        'FontName','Times New Roman', ...
        'Color','k', ...
        'FontSize',labelSize);

end

% ==============================
% Column headers
% ==============================

columnLabels = { ...
    '$C_{\mathrm{initial}} = 0.02$', ...
    '$C_{\mathrm{initial}} = 0.10$', ...
    '$C_{\mathrm{initial}} = 0.20$'};

for c = 1:3

    xCenter = left + ...
        (c-1)*(axWidth + gapX) + ...
        axWidth/2;

    annotation(fig,'textbox', ...
        [xCenter-0.09,0.925,0.18,0.035], ...
        'String',columnLabels{c}, ...
        'Interpreter','latex', ...
        'FontName','Times New Roman', ...
        'FontSize',headerSize, ...
        'FontWeight','normal', ...
        'Color','k', ...
        'HorizontalAlignment','center', ...
        'VerticalAlignment','middle', ...
        'EdgeColor','none');

end

% ==============================
% Row descriptors
% ==============================

rowLabels = { ...
    'Optimal timing', ...
    'Tumor burden', ...
    'Toxicity'};

for r = 1:3

    yCenter = bottom + ...
        (3-r)*(axHeight + gapY) + ...
        axHeight/2;

    annotation(fig,'textbox', ...
        [0.004,yCenter-0.02,0.075,0.04], ...
        'String',rowLabels{r}, ...
        'FontName','Times New Roman', ...
        'FontSize',headerSize, ...
        'FontWeight','normal', ...
        'HorizontalAlignment','center', ...
        'VerticalAlignment','middle', ...
        'Color','k', ...
        'EdgeColor','none');

end

% ==============================
% Shared colorbar — Optimal timing
% ==============================

cb1 = colorbar(ax1, ...
    'Position',[0.915 0.709 0.018 0.18]);

cb1.FontName = 'Times New Roman';
cb1.FontSize = fontSize;
cb1.Label.String = 'Anti-IL-6 administration time (hours)';
cb1.Label.FontName = 'Times New Roman';
cb1.Label.FontSize = labelSize;
cb1.Color = 'k';
cb1.Label.Color = 'k';

% ==============================
% Shared colorbar — Tumor
% ==============================

cb2 = colorbar(ax2, ...
    'Position',[0.915 0.423 0.018 0.18]);

cb2.FontName = 'Times New Roman';
cb2.FontSize = fontSize;
cb2.Label.String = 'Final tumor burden';
cb2.Label.FontName = 'Times New Roman';
cb2.Label.FontSize = labelSize;
cb2.Color = 'k';
cb2.Label.Color = 'k';

% ==============================
% Shared colorbar — Toxicity
% ==============================

cb3 = colorbar(ax3, ...
    'Position',[0.915 0.145 0.018 0.18]);

cb3.FontName = 'Times New Roman';
cb3.FontSize = fontSize;
cb3.Label.String = 'Normalized peak toxicity';
cb3.Label.FontName = 'Times New Roman';
cb3.Label.FontSize = labelSize;
cb3.Color = 'k';
cb3.Label.Color = 'k';

% ==============================
% Row descriptor
% ==============================

annotation(fig,'textbox', ...
    [0.005,bottom+0.22,0.025,0.45], ...
    'String',{'Outcome'}, ...
    'FontName','Times New Roman', ...
    'FontSize',11, ...
    'HorizontalAlignment','center', ...
    'VerticalAlignment','middle', ...
    'Rotation',90, ...
    'Color','k', ...
    'EdgeColor','none');

% ==============================
% Export
% ==============================

drawnow;

set(fig, ...
    'PaperUnits','inches', ...
    'PaperPosition',[0 0 12 9], ...
    'PaperSize',[12 9]);

print(fig, ...
    'Optimal_Anti_IL6_Timing_Outcomes.pdf', ...
    '-dpdf', ...
    '-painters');

exportgraphics(fig, ...
    'Optimal_Anti_IL6_Timing_Outcomes.png', ...
    'Resolution',600, ...
    'BackgroundColor','white');

%% Combined GSA Figure
% (A) Tumor Reduction Sobol Sensitivity
% (B) L0 vs L1 Sobol Sensitivity vs Initial Tumor Size
%
% Uses functions/results from:
% IL6DosingSobolGSA.m
% ============================================================

% Run GSA

% Main tumor-reduction and toxicity GSA
[S_tumor, ST_tumor, paramNames] = ...
    IL6DosingSobolGSA();

% L0 vs L1 GSA

theta_base = IL6DosingParameters();
theta_base.useIDosing.val = true;

C_initial_vals = linspace(0.01,0.8,15);
nC = numel(C_initial_vals);

S_L0  = zeros(nC,1);
ST_L0 = zeros(nC,1);
S_L1  = zeros(nC,1);
ST_L1 = zeros(nC,1);

N = 2;

for iC = 1:nC

    theta_base.C_initial.val = C_initial_vals(iC);

    [S,ST,paramNames] = ...
        runSobol(theta_base,1,N,'C_initial');

    idx_L0 = find(strcmp(paramNames,'L0'));
    idx_L1 = find(strcmp(paramNames,'L1'));

    S_L0(iC)  = S(idx_L0);
    ST_L0(iC) = ST(idx_L0);

    S_L1(iC)  = S(idx_L1);
    ST_L1(iC) = ST(idx_L1);

end

% Prevent small negative Sobol estimates
S_L0  = max(S_L0,0);
ST_L0 = max(ST_L0,0);
S_L1  = max(S_L1,0);
ST_L1 = max(ST_L1,0);


fontSize   = 11;
labelSize  = 13;
lineWidth  = 1.6;
markerSize = 5;

fig = figure( ...
    'Position',[300 100 850 850], ...
    'Color','w');

tiledlayout(2,1, ...
    'TileSpacing','compact', ...
    'Padding','compact');

% (A) Tumor Reduction
% ============================================================

ax1 = nexttile;

b = bar([S_tumor ST_tumor], ...
    'BarWidth',0.8);

b(1).FaceColor = [0.35 0.35 0.35];
b(2).FaceColor = [0.70 0.70 0.70];

b(1).EdgeColor = 'none';
b(2).EdgeColor = 'none';

set(ax1, ...
    'FontName','Times New Roman', ...
    'FontSize',fontSize, ...
    'XColor','k', ...
    'YColor','k', ...
    'LineWidth',1, ...
    'TickDir','out', ...
    'Box','off', ...
    'Layer','top', ...
    'XTick',1:length(paramNames), ...
    'XTickLabel',paramNames, ...
    'TickLabelInterpreter','none');

xtickangle(ax1,45);

xlabel(ax1, ...
    'Model parameter', ...
    'FontName','Times New Roman', ...
    'FontSize',labelSize, ...
    'Color','k');

ylabel(ax1, ...
    'Sobol Sensitivity Index', ...
    'FontName','Times New Roman', ...
    'FontSize',labelSize, ...
    'Color','k');

legend(ax1, ...
    {'First-order','Total-order'}, ...
    'FontName','Times New Roman', ...
    'FontSize',10, ...
    'Location','northeast', ...
    'Box','off');

xlim(ax1,[0.5 length(paramNames)+0.5]);

ylim(ax1, ...
    [0 max([S_tumor;ST_tumor])*1.10]);

grid(ax1,'on');
ax1.GridAlpha = 0.12;
ax1.GridColor = [0.5 0.5 0.5];
ax1.XGrid = 'off';
ax1.YGrid = 'on';

% Panel label
text(ax1, ...
    0.01,0.98,'(A)', ...
    'Units','normalized', ...
    'HorizontalAlignment','left', ...
    'VerticalAlignment','top', ...
    'FontName','Times New Roman', ...
    'FontSize',13, ...
    'Color','k');


% (B) L0 vs L1 Sensitivity
% ============================================================

ax2 = nexttile;

hold(ax2,'on');

p1 = plot(ax2, ...
    C_initial_vals,S_L0, ...
    '-o', ...
    'Color',[0.20 0.20 0.20], ...
    'LineWidth',lineWidth, ...
    'MarkerSize',markerSize, ...
    'MarkerFaceColor','w');

p2 = plot(ax2, ...
    C_initial_vals,S_L1, ...
    '-o', ...
    'Color','k', ...
    'LineWidth',lineWidth, ...
    'MarkerSize',markerSize, ...
    'MarkerFaceColor','w');

p3 = plot(ax2, ...
    C_initial_vals,ST_L0, ...
    '--', ...
    'Color',[0.20 0.20 0.20], ...
    'LineWidth',lineWidth);

p4 = plot(ax2, ...
    C_initial_vals,ST_L1, ...
    '--', ...
    'Color','k', ...
    'LineWidth',lineWidth);

set(ax2, ...
    'FontName','Times New Roman', ...
    'FontSize',fontSize, ...
    'XColor','k', ...
    'YColor','k', ...
    'LineWidth',1, ...
    'TickDir','out', ...
    'Box','off', ...
    'Layer','top');

xlabel(ax2, ...
    'Initial Tumor Density, $C_{\mathrm{initial}}$', ...
    'Interpreter','latex', ...
    'FontName','Times New Roman', ...
    'FontSize',labelSize, ...
    'Color','k');

ylabel(ax2, ...
    'Sobol Sensitivity Index', ...
    'FontName','Times New Roman', ...
    'FontSize',labelSize, ...
    'Color','k');

legend(ax2, ...
    [p1 p2 p3 p4], ...
    {'$S_{L_0}$','$S_{L_1}$','$S_{T,L_0}$','$S_{T,L_1}$'}, ...
    'Interpreter','latex', ...
    'FontName','Times New Roman', ...
    'FontSize',10, ...
    'Location','best', ...
    'Box','off');

xlim(ax2, ...
    [min(C_initial_vals) max(C_initial_vals)]);

ylim(ax2, ...
    [0 max([S_L0;ST_L0;S_L1;ST_L1])*1.10]);

grid(ax2,'on');
ax2.GridAlpha = 0.12;
ax2.GridColor = [0.5 0.5 0.5];
ax2.XGrid = 'off';
ax2.YGrid = 'on';

% Panel label
text(ax2, ...
    0.01,0.98,'(B)', ...
    'Units','normalized', ...
    'HorizontalAlignment','left', ...
    'VerticalAlignment','top', ...
    'FontName','Times New Roman', ...
    'FontSize',13, ...
    'Color','k');

hold(ax2,'off');


%% ============================================================
% Export
% ============================================================

drawnow;

set(fig,'InvertHardcopy','off');

exportgraphics(fig, ...
    'Combined_GSA_Tumor_L0_L1.pdf', ...
    'ContentType','vector', ...
    'BackgroundColor','white');

exportgraphics(fig, ...
    'Combined_GSA_Tumor_L0_L1.png', ...
    'Resolution',600, ...
    'BackgroundColor','white');