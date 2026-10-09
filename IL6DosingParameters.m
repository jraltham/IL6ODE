
function theta = IL6DosingParameters()
    
    theta = struct();
    
    % Dosing schedule
    theta.tEnd.val = 4000; %(hours)
    
    % --- Anti-CTLA-4 dosing ---
    theta.useIDosing.val   = false;
    theta.doseStart_I.val  = 1;
    theta.doseInterval_I.val = 24*7*3; %every 3 weeks
    theta.nDoses_I.val    = 4;
    theta.doseI.val        = 1;
    
    % --- Anti-IL-6 dosing ---
    theta.useADosing.val   = false;
    theta.doseStart_A.val  = 500;        % delayed start (hours)
    theta.doseInterval_A.val = 24*7*2;
    theta.nDoses_A.val     = 1;
    theta.doseA.val        = 1;
    
    % Initial Conditions
    theta.C_initial.val = 0.1;
    %theta.C_initial.min   = 0.05;
    %theta.C_initial.max   = 0.3;
    theta.E_initial.val = 0.1;
    theta.E_initial.min   = 0.01;
    theta.E_initial.max   = 0.5;
    theta.L_initial.val = 1;
    theta.L_initial.min   = 0.1;
    theta.L_initial.max   = 5;
    theta.A_initial.val = 0;
    theta.I_initial.val = 0;
    theta.X_initial.val = 0;
    
    % Equation Parameters
    
    theta.rT.val   = 0.000387; %5e-4%1.2535e-4;          % tumor growth rate (per day)
    theta.kT.val   = 1;          % tumor carrying capacity (cells)
    theta.kE.val   = 50e-5;   % Range (4.16e-5, 83.2e-5)   % killing efficiency per effector cell
    %theta.kE.min   = 4.16e-5;
    %theta.kE.max   = 83.2e-5;
   
    
    theta.kA.val   = 1.56e-3;     % Range (1.56e-3 , 2e-2)  % baseline T cell activation rate (per hour)
    %theta.kA.min   = 1.56e-3;
    %theta.kA.max   = 2e-2;
    theta.eta.val = 1;         % Range (0,10) IL-6 modulatory strength on activation
    
    theta.L0.val   = 0.01;            % Range (0.01,10) IL-6 half-max activation (pg/mL)
    theta.L0.min   = 0.1;
    theta.L0.max   = 10;
    theta.L1.val   = 5;           % Range (5, 10e3) IL-6 half-max suppression (pg/mL)
    theta.L1.min   = 2;
    theta.L1.max   = 10e3;
    
    theta.n.val    = 1;            % Range (1,3) Hill coefficient for activation
    %theta.n.min   = 1;
    %theta.n.max   = 3;
    theta.m.val    = 1;            % Range (1,3) Hill coefficient for suppression
    %theta.m.min   = 1;
    %theta.m.max   = 3;
    
    theta.delta.val = 0.00404%.0041165;       % Range (0.000833, 0.0074) effector cell death rate (ub 0.0074)
    %theta.delta.min   = .000833;
    %theta.delta.max   = 0.0074;
    
    theta.pL.val   = .2499;         % Range (0, 958.33) IL-6 production rate
    theta.pL.min   = 0;
    theta.pL.max   = 958.33;
    theta.phi.val  = 10;         % Range (0, 325) tumor contribution relative to T cells
    theta.phi.min   = 0;
    theta.phi.max   = 325;
    
    theta.kC6.val  = 0.0463;        % IL-6 clearance (per hour)
    theta.kB.val   = 50;         % Range (27.43, 2743) IL-6–anti-IL6 binding rate
    theta.kB.min   = 27.43;
    theta.kB.max   = 2743;
    theta.KIL6.val = 161.7;         % IL-6 saturation rate
    
    theta.kCE.val  = 0.00222;          % anti-IL-6 clearance
    theta.kCI.val  = 0.00222;          % anti-CTLA4 clearance
    
    theta.kICI.val = 0.0728;          % half-maximal CTLA-4 immune activation
    
    theta.s.val = 1;        %Immune suppression coefficient (0 implies complete immune shutdown)
    
    theta.g.val = 0.0202;     %Tumor size when effector cell recruitment is half maximal
    theta.p.val = 0.00518;    %Maximal activation of effector cells in response to Tumor cell Lysis
    theta.q.val = 1.43e-2;    %CD8 cell inactivation rate by tumor cells.
    
    theta.rx.val = .866 % .662336;          %Rate at which inflammatory conditions translate into toxicity
    theta.rx.min = 0.4 %4e-6;
    theta.rx.max = 2;
    theta.rxres.val =  0.7539 %9.02116e-04;    % Range (0.000541, 0.00138) Toxicity remission rate
    theta.rxres.min = 0.5%1e-6; %5e-4;
    theta.rxres.max = 1;%1.4e-3;
    
    theta.Emax.val = 0;
    theta.beta.val = 0;
    theta.ltox.val = 20.4 %9.83;        %IL-6 concentration at which toxicity induction is half maximal
    theta.ltox.min = 1;
    theta.ltox.max = 30;
    theta.etox.val = 2.63 %1;           %Effector cell concentration at which toxicity induction is half maximal
    theta.etox.min = 0.1;
    theta.etox.max = 5;
end
