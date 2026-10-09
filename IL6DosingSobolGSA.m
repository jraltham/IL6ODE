%% ============================================================
% IL6 Model – Sobol Sensitivity Analysis
% ============================================================

clear; clc;

%% ============================================================
% SECTION 1 — Standard Sobol (All Parameters Free)
% ============================================================

theta_base = IL6DosingParameters();
theta_base.useIDosing.val = true;

N = 10;   % increase to 800–1000 for publication

% Tumor Reduction
[S_tumor, ST_tumor, paramNames] = runSobol(theta_base,1,N,'');

% Peak Toxicity
[S_tox, ST_tox] = runSobol(theta_base,2,N,'');

%% ---- Plot Tumor ----
figure;
bar([S_tumor ST_tumor])
set(gca,'XTickLabel',paramNames,'XTickLabelRotation',45)
ylabel('Sobol Index')
title('Tumor Reduction Sensitivity')
legend('First Order','Total Order')
grid on

% ---- Plot Toxicity ----
figure;
bar([S_tox ST_tox])
set(gca,'XTickLabel',paramNames,'XTickLabelRotation',45)
ylabel('Sobol Index')
title('Peak Toxicity Sensitivity')
legend('First Order','Total Order')
grid on

%% ---- Plot Tumor ----
figure('Position',[100 100 1400 600])   % wider figure

bar([S_tumor ST_tumor])

ax = gca;
ax.XTick = 1:length(paramNames);
ax.XTickLabel = paramNames;
ax.XTickLabelRotation = 45;
ax.FontSize = 11;

ylabel('Sobol Index')
title('Tumor Reduction Sensitivity')
legend('First Order','Total Order')
grid on

% improve spacing
xtickangle(45)
ax.TickLabelInterpreter = 'none';

% optional: tighten axes
xlim([0.5 length(paramNames)+0.5])

%% ---- Plot Toxicity ----
figure('Position',[100 100 1400 600])

bar([S_tox ST_tox])

ax = gca;
ax.XTick = 1:length(paramNames);
ax.XTickLabel = paramNames;
ax.XTickLabelRotation = 45;
ax.FontSize = 11;

ylabel('Sobol Index')
title('Peak Toxicity Sensitivity')
legend('First Order','Total Order')
grid on

xtickangle(45)
ax.TickLabelInterpreter = 'none';

xlim([0.5 length(paramNames)+0.5])

%% Figure 3 A
% Official SOBOL

% Standard Sobol Sensitivity Analysis — Publication Figure
% ============================================================

theta_base = IL6DosingParameters();
theta_base.useIDosing.val = true;

% Sobol sample size
% Use N = 800–1000+ for final publication analysis
N = 20;

% ============================================================
% Run Sobol analyses
% ============================================================

% Tumor reduction
[S_tumor, ST_tumor, paramNames] = runSobol(theta_base,1,N,'');

% Peak toxicity
[S_tox, ST_tox] = runSobol(theta_base,2,N,'');

% ============================================================
% Prevent small negative Sobol estimates
% ============================================================

S_tumor  = max(S_tumor,0);
ST_tumor = max(ST_tumor,0);

S_tox  = max(S_tox,0);
ST_tox = max(ST_tox,0);

% ============================================================
% Figure parameters
% ============================================================

fontSize   = 12;
labelSize  = 14;
titleSize  = 15;
lineWidth  = 1.2;

nParams = length(paramNames);

% ============================================================
% Create figure
% ============================================================

fig = figure('Position',[100 100 1500 650], ...
    'Color','w');

t = tiledlayout(1,2, ...
    'TileSpacing','compact', ...
    'Padding','compact');

% ============================================================
% Panel A — Tumor Reduction
% ============================================================

ax1 = nexttile;

b1 = bar([S_tumor ST_tumor], ...
    'BarWidth',0.8);

% Professional grayscale / muted appearance
b1(1).FaceColor = [0.35 0.35 0.35];
b1(2).FaceColor = [0.70 0.70 0.70];

b1(1).EdgeColor = 'none';
b1(2).EdgeColor = 'none';

% Axis formatting
set(ax1, ...
    'FontName','Times New Roman', ...
    'FontSize',fontSize, ...
    'LineWidth',1, ...
    'TickDir','out', ...
    'Box','off', ...
    'Layer','top', ...
    'XTick',1:nParams, ...
    'XTickLabel',paramNames, ...
    'TickLabelInterpreter','none');

xtickangle(45);

xlabel('Model parameter', ...
    'FontName','Times New Roman', ...
    'FontSize',labelSize);

ylabel('Sobol Sensitivity Index', ...
    'FontName','Times New Roman', ...
    'FontSize',labelSize);

title('Tumor Reduction', ...
    'FontName','Times New Roman', ...
    'FontSize',titleSize, ...
    'FontWeight','normal');

legend({'First-order','Total-order'}, ...
    'FontName','Times New Roman', ...
    'FontSize',11, ...
    'Location','northeast', ...
    'Box','off');

xlim([0.5 nParams+0.5]);

ylim([0 max([S_tumor;ST_tumor])*1.10]);

% Light horizontal reference lines
grid on;
ax1.GridAlpha = 0.12;
ax1.GridColor = [0.5 0.5 0.5];
ax1.XGrid = 'off';
ax1.YGrid = 'on';

% ============================================================
% Panel B — Peak Toxicity
% ============================================================

ax2 = nexttile;

b2 = bar([S_tox ST_tox], ...
    'BarWidth',0.8);

b2(1).FaceColor = [0.35 0.35 0.35];
b2(2).FaceColor = [0.70 0.70 0.70];

b2(1).EdgeColor = 'none';
b2(2).EdgeColor = 'none';

% Axis formatting
set(ax2, ...
    'FontName','Times New Roman', ...
    'FontSize',fontSize, ...
    'LineWidth',1, ...
    'TickDir','out', ...
    'Box','off', ...
    'Layer','top', ...
    'XTick',1:nParams, ...
    'XTickLabel',paramNames, ...
    'TickLabelInterpreter','none');

xtickangle(45);

xlabel('Model parameter', ...
    'FontName','Times New Roman', ...
    'FontSize',labelSize);

ylabel('Sobol Sensitivity Index', ...
    'FontName','Times New Roman', ...
    'FontSize',labelSize);

title('Peak Toxicity', ...
    'FontName','Times New Roman', ...
    'FontSize',titleSize, ...
    'FontWeight','normal');

legend({'First-order','Total-order'}, ...
    'FontName','Times New Roman', ...
    'FontSize',11, ...
    'Location','northeast', ...
    'Box','off');

xlim([0.5 nParams+0.5]);

ylim([0 max([S_tox;ST_tox])*1.10]);

grid on;
ax2.GridAlpha = 0.12;
ax2.GridColor = [0.5 0.5 0.5];
ax2.XGrid = 'off';
ax2.YGrid = 'on';

% ============================================================
% Export
% ============================================================

exportgraphics(fig, ...
    'Sobol_Sensitivity_Tumor_Toxicity.pdf', ...
    'ContentType','vector');

exportgraphics(fig, ...
    'Sobol_Sensitivity_Tumor_Toxicity.png', ...
    'Resolution',600);
%% ============================================================
% SECTION 2 — Sensitivity of L0 and L1 vs Initial Tumor Size
% ============================================================

theta_base = IL6DosingParameters();
theta_base.useIDosing.val = true;

C0_vals = linspace(0.01,0.8,15);
nC = numel(C0_vals);

S_L0  = zeros(nC,1);
ST_L0 = zeros(nC,1);
S_L1  = zeros(nC,1);
ST_L1 = zeros(nC,1);

N = 1200;

for iC = 1:nC

    fprintf('Running C0 = %.4f\n', C0_vals(iC));

    theta_base.C_initial.val = C0_vals(iC);

    [S,ST,paramNames] = runSobol(theta_base,1,N,'C_initial');

    idx_L0 = find(strcmp(paramNames,'L0'));
    idx_L1 = find(strcmp(paramNames,'L1'));

    S_L0(iC)  = S(idx_L0);
    ST_L0(iC) = ST(idx_L0);
    S_L1(iC)  = S(idx_L1);
    ST_L1(iC) = ST(idx_L1);
end

figure('Position',[300 300 800 450])

%semilogx(C0_vals, S_L0, '-o','LineWidth',2); hold on
%semilogx(C0_vals, S_L1, '-o','LineWidth',2);
%semilogx(C0_vals, ST_L0, '--','LineWidth',2);
%semilogx(C0_vals, ST_L1, '--','LineWidth',2);

plot(C0_vals, S_L0, '-o','LineWidth',2); hold on
plot(C0_vals, S_L1, '-o','LineWidth',2);
plot(C0_vals, ST_L0, '--','LineWidth',2);
plot(C0_vals, ST_L1, '--','LineWidth',2);

xlabel('Initial Tumor Size C_0')
ylabel('Sobol Index')
title('Sensitivity of L0 and L1 vs Initial Tumor Size')

legend('S L0','S L1','ST L0','ST L1','Location','best')
grid on

%% Figure 3B
% Sensitivity of L0 and L1 vs Initial Tumor size
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

    fprintf('Running C_initial = %.4f\n', C_initial_vals(iC));

    theta_base.C_initial.val = C_initial_vals(iC);

    [S,ST,paramNames] = runSobol(theta_base,1,N,'C_initial');

    idx_L0 = find(strcmp(paramNames,'L0'));
    idx_L1 = find(strcmp(paramNames,'L1'));

    S_L0(iC)  = S(idx_L0);
    ST_L0(iC) = ST(idx_L0);
    S_L1(iC)  = S(idx_L1);
    ST_L1(iC) = ST(idx_L1);

end

% Clamp negative Sobol indices to zero
S_L0  = max(S_L0,0);
ST_L0 = max(ST_L0,0);
S_L1  = max(S_L1,0);
ST_L1 = max(ST_L1,0);

% Figure parameters
fontSize  = 12;
labelSize = 14;
lineWidth = 1.8;
markerSize = 6;

fig =figure('Position',[300 300 800 450], ...
       'Color','w');

hold on;

% First-order Sobol indices
p1 = plot(C_initial_vals,S_L0, ...
    '-o', ...
    'LineWidth',lineWidth, ...
    'MarkerSize',markerSize);

p2 = plot(C_initial_vals,S_L1, ...
    '-o', ...
    'LineWidth',lineWidth, ...
    'MarkerSize',markerSize);

% Total-order Sobol indices
p3 = plot(C_initial_vals,ST_L0, ...
    '--', ...
    'LineWidth',lineWidth);

p4 = plot(C_initial_vals,ST_L1, ...
    '--', ...
    'LineWidth',lineWidth);

% Axis formatting
ax = gca;

set(ax, ...
    'FontName','Times New Roman', ...
    'FontSize',fontSize, ...
    'LineWidth',1, ...
    'TickDir','out', ...
    'Box','off', ...
    'Layer','top');

% Labels
xlabel('Initial Tumor Density, $C_{\mathrm{initial}}$', ...
    'Interpreter','latex', ...
    'FontName','Times New Roman', ...
    'FontSize',labelSize);

ylabel('Sobol Sensitivity Index', ...
    'FontName','Times New Roman', ...
    'FontSize',labelSize);

% Optional title — remove for manuscript figure
% title('Sensitivity of IL-6 Pathway Parameters', ...
%     'FontName','Times New Roman', ...
%     'FontSize',labelSize, ...
%     'FontWeight','normal');

% Legend
legend([p1 p2 p3 p4], ...
    {'$S_{L_0}$','$S_{L_1}$','$S_{T,L_0}$','$S_{T,L_1}$'}, ...
    'Interpreter','latex', ...
    'FontName','Times New Roman', ...
    'FontSize',11, ...
    'Location','best', ...
    'Box','off');

% Grid
grid on;
ax.GridAlpha = 0.15;
ax.MinorGridAlpha = 0.08;

% Axis limits
xlim([min(C_initial_vals) max(C_initial_vals)]);

% Ensure sensitivity starts at zero
ylim([0 max([S_L0; ST_L0; S_L1; ST_L1])*1.1]);


% Export
exportgraphics(fig,'L0_L1_Sobol_vs_TumorSize.pdf','ContentType','vector');
exportgraphics(fig,'L0_L1_Sobol_vs_TumorSize.png','Resolution',600);



%% ============================================================
% ------------------ CORE SOBOL FUNCTION ---------------------
%% ============================================================

function [S,ST,paramNames] = runSobol(theta_base,outputType,N,excludeParam)

    % ---- Parameter Selection ----
    fields = fieldnames(theta_base);
    paramNames = {};

    for i = 1:length(fields)

        hasRange = isfield(theta_base.(fields{i}),'min') && ...
                   isfield(theta_base.(fields{i}),'max');

        if hasRange
            if isempty(excludeParam) || ~strcmp(fields{i},excludeParam)
                paramNames{end+1} = fields{i};
            end
        end
    end

    d = length(paramNames);

    % ---- Sobol Sampling ----
    sob = sobolset(2*d);
    sob = scramble(sob,'MatousekAffineOwen');
    X = net(sob,N);

    A = X(:,1:d);
    B = X(:,d+1:2*d);

    for j = 1:d
        low  = theta_base.(paramNames{j}).min;
        high = theta_base.(paramNames{j}).max;

        A(:,j) = low + A(:,j)*(high-low);
        B(:,j) = low + B(:,j)*(high-low);
    end

    % ---- Evaluate Model ----
    YA = zeros(N,1);
    YB = zeros(N,1);

    for i = 1:N
        YA(i) = modelOutput(A(i,:),paramNames,theta_base,outputType);
        YB(i) = modelOutput(B(i,:),paramNames,theta_base,outputType);
    end

    YAB = zeros(N,d);

    for j = 1:d
        AB = A;
        AB(:,j) = B(:,j);

        for i = 1:N
            YAB(i,j) = modelOutput(AB(i,:),paramNames,theta_base,outputType);
        end
    end

    if any(isnan(YA)) || any(isinf(YA))
        error('YA contains NaN or Inf')
    end
    
    if any(isnan(YB)) || any(isinf(YB))
        error('YB contains NaN or Inf')
    end
    
    VY = var([YA;YB],1);

    if VY < 1e-10
        S  = zeros(d,1);
        ST = zeros(d,1);
        return
    end
    % ---- Sobol Indices ----
    [S,ST] = sobolIndices(YA,YB,YAB);
end



%% ============================================================
% ------------------ MODEL OUTPUT -----------------------------
%% ============================================================

function y = modelOutput(paramSet,paramNames,theta_template,outputType)

    theta_i = theta_template;

    for j = 1:length(paramNames)
        theta_i.(paramNames{j}).val = paramSet(j);
    end

    params_i = IL6ParamUnpack(theta_i);
    [~, x_i] = IL6RunSim(params_i);

    C = x_i(:,1);
    X = x_i(:,6);

    switch outputType
        case 1
            y = (C(1) - C(end)) / C(1);
        case 2
            y = max(X);
    end

    if any(isnan(x_i(:))) || any(isinf(x_i(:)))
        y = 0;
        return
    end
end



%% ============================================================
% ------------------ SOBOL INDEX FORMULAS ---------------------
%% ============================================================

function [S,ST] = sobolIndices(YA,YB,YAB)

    VY = var([YA;YB],1);

    d = size(YAB,2);

    S  = zeros(d,1);
    ST = zeros(d,1);

    for j = 1:d
        S(j)  = mean(YB .* (YAB(:,j) - YA)) / VY;
        ST(j) = mean((YA - YAB(:,j)).^2) / (2*VY);
    end
end


%% ============================================================
% SECTION 2 — Sensitivity of L Initial vs Initial Tumor Size
% ============================================================
clear;
clc;
theta_base = IL6DosingParameters();
theta_base.useIDosing.val = true;

C0_vals = linspace(0.01,0.8,10);
nC = numel(C0_vals);

S_Linit  = zeros(nC,1);
ST_Linit = zeros(nC,1);

N = 2000;

for iC = 1:nC

    fprintf('Running C0 = %.4f\n', C0_vals(iC));

    theta_base.C_initial.val = C0_vals(iC);

    [S,ST,paramNames] = runSobol(theta_base,1,N,'C_initial');

    idx_L_init = find(strcmp(paramNames,'L_initial'));

    S_L_init(iC)  = S(idx_L_init);
    ST_L_init(iC) = ST(idx_L_init);
end
% ---- Plot ----
figure('Position',[300 300 800 450])

%semilogx(C0_vals, S_L0, '-o','LineWidth',2); hold on
%semilogx(C0_vals, S_L1, '-o','LineWidth',2);
%semilogx(C0_vals, ST_L0, '--','LineWidth',2);
%semilogx(C0_vals, ST_L1, '--','LineWidth',2);

plot(C0_vals, S_L_init, '-o','LineWidth',2); hold on
plot(C0_vals, ST_L_init, '--','LineWidth',2);

xlabel('Initial Tumor Size C_0')
ylabel('Sobol Index')
title('Sensitivity of Initial IL-6 vs Initial Tumor Size')

legend('S L Initial', 'ST L Initial','Location','best')
grid on

%% 
% Sensitivity of Initial IL6 vs Initial Tumor Size

theta_base = IL6DosingParameters();
theta_base.useIDosing.val = true;

C_initial_vals = linspace(0.01,0.8,10);
nC = numel(C_initial_vals);

S_L_init  = zeros(nC,1);
ST_L_init = zeros(nC,1);

N = 8;

for iC = 1:nC

    fprintf('Running C_initial = %.4f\n', C_initial_vals(iC));

    theta_base.C_initial.val = C_initial_vals(iC);

    [S,ST,paramNames] = runSobol(theta_base,1,N,'C_initial');

    idx_L_init = find(strcmp(paramNames,'L_initial'));

    S_L_init(iC)  = S(idx_L_init);
    ST_L_init(iC) = ST(idx_L_init);

end

% Clamp negative Sobol estimates to zero
S_L_init  = max(S_L_init,0);
ST_L_init = max(ST_L_init,0);

% Figure parameters
fontSize   = 12;
labelSize  = 14;
lineWidth  = 1.8;
markerSize = 6;

fig = figure('Position',[300 300 800 450], ...
             'Color','w');

hold on;

% First-order Sobol index
p1 = plot(C_initial_vals,S_L_init, ...
    '-o', ...
    'LineWidth',lineWidth, ...
    'MarkerSize',markerSize);

% Total-order Sobol index
p2 = plot(C_initial_vals,ST_L_init, ...
    '--', ...
    'LineWidth',lineWidth);

% Axis formatting
ax = gca;

set(ax, ...
    'FontName','Times New Roman', ...
    'FontSize',fontSize, ...
    'LineWidth',1, ...
    'TickDir','out', ...
    'Box','off', ...
    'Layer','top');

% Labels
xlabel('Initial Tumor Density, $C_{\mathrm{initial}}$', ...
    'Interpreter','latex', ...
    'FontName','Times New Roman', ...
    'FontSize',labelSize);

ylabel('Sobol Sensitivity Index', ...
    'FontName','Times New Roman', ...
    'FontSize',labelSize);

% Optional title — remove for manuscript figure
% title('Sensitivity of Initial IL-6 Concentration', ...
%     'FontName','Times New Roman', ...
%     'FontSize',labelSize, ...
%     'FontWeight','normal');

% Legend
legend([p1 p2], ...
    {'$S_{L_{\mathrm{initial}}}$','$S_{T,L_{\mathrm{initial}}}$'}, ...
    'Interpreter','latex', ...
    'FontName','Times New Roman', ...
    'FontSize',11, ...
    'Location','best', ...
    'Box','off');

% Grid
grid on;
ax.GridAlpha = 0.15;
ax.MinorGridAlpha = 0.08;

% Axis limits
xlim([min(C_initial_vals) max(C_initial_vals)]);

ymax = max([S_L_init; ST_L_init]);

if ymax > 0
    ylim([0 ymax*1.1]);
else
    ylim([0 1]);
end

hold off;

% Export
exportgraphics(fig, ...
    'L_initial_Sobol_vs_TumorSize.pdf', ...
    'ContentType','vector');

exportgraphics(fig, ...
    'L_initial_Sobol_vs_TumorSize.png', ...
    'Resolution',600);


%% Figure 3 (complete)
% Combined GSA Figure
% (A) Tumor Reduction Sobol Sensitivity
% (B) L0 vs L1 Sobol Sensitivity vs Initial Tumor Size
% ============================================================

% PART A — Tumor Reduction Sobol Sensitivity
% ============================================================

theta_base = IL6DosingParameters();
theta_base.useIDosing.val = true;

% Sobol sample size
N_tumor = 512 %2048;

% Run Sobol analysis for tumor reduction
[S_tumor, ST_tumor, paramNames] = ...
    runSobol(theta_base,1,N_tumor,'');

% Prevent small negative Sobol estimates
S_tumor  = max(S_tumor,0);
ST_tumor = max(ST_tumor,0);

nParams = length(paramNames);


% PART B — L0 and L1 Sensitivity vs Initial Tumor Size
% ============================================================

C_initial_vals = linspace(0.01,0.8,15);
nC = numel(C_initial_vals);

S_L0  = zeros(nC,1);
ST_L0 = zeros(nC,1);
S_L1  = zeros(nC,1);
ST_L1 = zeros(nC,1);

N_L0L1 = 1024;

for iC = 1:nC

    fprintf('Running C_initial = %.4f\n',C_initial_vals(iC));

    theta_base.C_initial.val = C_initial_vals(iC);

    [S,ST,paramNames_L0L1] = ...
        runSobol(theta_base,1,N_L0L1,'C_initial');

    idx_L0 = find(strcmp(paramNames_L0L1,'L0'));
    idx_L1 = find(strcmp(paramNames_L0L1,'L1'));

    S_L0(iC)  = S(idx_L0);
    ST_L0(iC) = ST(idx_L0);

    S_L1(iC)  = S(idx_L1);
    ST_L1(iC) = ST(idx_L1);

end

% Clamp negative Sobol indices to zero
S_L0  = max(S_L0,0);
ST_L0 = max(ST_L0,0);
S_L1  = max(S_L1,0);
ST_L1 = max(ST_L1,0);


% FIGURE FORMATTING
% ============================================================

fontSize   = 11;
labelSize  = 13;
lineWidth  = 1.6;
markerSize = 5;

fig = figure( ...
    'Position',[300 100 850 850], ...
    'Color','w', ...
    'Units','pixels', ...
    'InvertHardcopy','off');

t = tiledlayout(fig,2,1, ...
    'TileSpacing','compact', ...
    'Padding','compact');


% PANEL A — TUMOR REDUCTION
% ============================================================

ax1 = nexttile(t);

b = bar(ax1,[S_tumor ST_tumor], ...
    'BarWidth',0.80);

% Grayscale publication appearance
b(1).FaceColor = [0.35 0.35 0.35];
b(2).FaceColor = [0.70 0.70 0.70];

b(1).EdgeColor = 'none';
b(2).EdgeColor = 'none';

% White plot background
ax1.Color = 'w';

% Axis formatting
set(ax1, ...
    'FontName','Times New Roman', ...
    'FontSize',fontSize, ...
    'FontWeight','normal', ...
    'XColor','k', ...
    'YColor','k', ...
    'LineWidth',1, ...
    'TickDir','out', ...
    'Box','off', ...
    'Layer','top', ...
    'XTick',1:nParams, ...
    'XTickLabel',paramNames, ...
    'TickLabelInterpreter','none');

xtickangle(ax1,45);

xlabel(ax1, ...
    'Model parameter', ...
    'FontName','Times New Roman', ...
    'FontSize',labelSize, ...
    'FontWeight','normal', ...
    'Color','k');

ylabel(ax1, ...
    'Sobol Sensitivity Index', ...
    'FontName','Times New Roman', ...
    'FontSize',labelSize, ...
    'FontWeight','normal', ...
    'Color','k');

title(ax1, ...
    'Global Sensitivity on Tumor Reduction', ...
    'FontName','Times New Roman', ...
    'FontSize',16, ...
    'FontWeight','normal', ...
    'Color','k');

xlim(ax1,[0.5 nParams+0.5]);

ylim(ax1, ...
    [0 max([S_tumor;ST_tumor])*1.10]);

% Legend
lg1 = legend(ax1, ...
    {'First-order','Total-order'}, ...
    'FontName','Times New Roman', ...
    'FontSize',10, ...
    'FontWeight','normal', ...
    'Location','northeast', ...
    'Box','off');

lg1.TextColor = 'k';

% Grid
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
    'FontWeight','normal', ...
    'Color','k');


% PANEL B — L0 vs L1 SENSITIVITY
% ============================================================

ax2 = nexttile(t);

hold(ax2,'on');

% ------------------------------------------------------------
% Color palette
% ------------------------------------------------------------

% L0 = blue
% L1 = red
color_L0 = [0.12 0.35 0.60];
color_L1 = [0.70 0.20 0.20];


% ------------------------------------------------------------
% First-order Sobol indices
% ------------------------------------------------------------

p1 = plot(ax2, ...
    C_initial_vals,S_L0, ...
    '-o', ...
    'Color',color_L0, ...
    'LineWidth',lineWidth, ...
    'MarkerSize',markerSize, ...
    'MarkerFaceColor','w', ...
    'MarkerEdgeColor',color_L0);

p2 = plot(ax2, ...
    C_initial_vals,S_L1, ...
    '-o', ...
    'Color',color_L1, ...
    'LineWidth',lineWidth, ...
    'MarkerSize',markerSize, ...
    'MarkerFaceColor','w', ...
    'MarkerEdgeColor',color_L1);


% ------------------------------------------------------------
% Total-order Sobol indices
% ------------------------------------------------------------

p3 = plot(ax2, ...
    C_initial_vals,ST_L0, ...
    '--', ...
    'Color',color_L0, ...
    'LineWidth',lineWidth);

p4 = plot(ax2, ...
    C_initial_vals,ST_L1, ...
    '--', ...
    'Color',color_L1, ...
    'LineWidth',lineWidth);


% White plot background
ax2.Color = 'w';

% Axis formatting
set(ax2, ...
    'FontName','Times New Roman', ...
    'FontSize',fontSize, ...
    'FontWeight','normal', ...
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
    'FontWeight','normal', ...
    'Color','k');

ylabel(ax2, ...
    'Sobol Sensitivity Index', ...
    'FontName','Times New Roman', ...
    'FontSize',labelSize, ...
    'FontWeight','normal', ...
    'Color','k');

% Legend
lg2 = legend(ax2, ...
    [p1 p2 p3 p4], ...
    {'$S_{L_0}$','$S_{L_1}$','$S_{T,L_0}$','$S_{T,L_1}$'}, ...
    'Interpreter','latex', ...
    'FontName','Times New Roman', ...
    'FontSize',10, ...
    'FontWeight','normal', ...
    'Location','best', ...
    'Box','off');

lg2.TextColor = 'k';


% Grid
grid(ax2,'on');

ax2.GridAlpha = 0.12;
ax2.GridColor = [0.5 0.5 0.5];
ax2.XGrid = 'off';
ax2.YGrid = 'on';


% Axis limits
xlim(ax2, ...
    [min(C_initial_vals) max(C_initial_vals)]);

ylim(ax2, ...
    [0 max([S_L0;ST_L0;S_L1;ST_L1])*1.10]);


% Panel label
text(ax2, ...
    0.01,0.98,'(B)', ...
    'Units','normalized', ...
    'HorizontalAlignment','left', ...
    'VerticalAlignment','top', ...
    'FontName','Times New Roman', ...
    'FontSize',13, ...
    'FontWeight','normal', ...
    'Color','k');

hold(ax2,'off');


%% ============================================================
% EXPORT
% ============================================================

drawnow;

% Vector PDF
exportgraphics(fig, ...
    'Combined_GSA_Tumor_L0_L1.pdf', ...
    'ContentType','vector', ...
    'BackgroundColor','white');

% High-resolution PNG
exportgraphics(fig, ...
    'Combined_GSA_Tumor_L0_L1.png', ...
    'Resolution',600, ...
    'BackgroundColor','white');