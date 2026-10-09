%% ============================
% ODE system
% ============================
function dx = IL6DosingODE(t, x, params)

    % States
    C = x(1);
    E = x(2);
    L = x(3);
    A = x(4);
    I = x(5);
    X = x(6);
   
    % Parameters

    E_initial=params.E_initial;
    C_initial=params.C_initial;
    L_initial=params.L_initial;
    
    rT = params.rT;
    kT = params.kT;
    kE = params.kE;
    kA = params.kA;
    L0 = params.L0;
    L1 = params.L1;
    n = params.n;
    m = params.m;
    delta = params.delta;
    pL = params.pL;
    phi = params.phi;
    kC6 = params.kC6;
    kB = params.kB;
    KIL6 = params.KIL6;
    kCE = params.kCE;
    kCI = params.kCI;
    kICI = params.kICI;
    s = params.s;
    g = params.g;
    p = params.p;
    q = params.q;
    rx = params.rx;
    rxres = params.rxres;
    ltox = params.ltox;
    etox = params.etox;
    
   
    % Drug effects
    Fi = I / (kICI + I);
    A_L = (L/L0).^n ./ (1 + (L/L0).^n);   % activation
    S_L = (L/L1).^m ./ (1 + (L/L1).^m);   % suppression
    
    Fl_raw = A_L .* (1 - S_L);   % ∈ [0,1]

    Edrive = (E - E_initial) + phi*(C - C_initial);
    Edrive_pos = max(Edrive, 0);

    dL_act = max(L - L_initial, 0);
    dE_act = max(E - E_initial, 0);
    
    % Differential equations
    dC = rT * C * (1- (C/kT)) - kE * E * C;     % Tumor cell population rate of change
    dE = s * kA * (Fl_raw) * (0.1 + Fi) + (p*E*C)/(g+C)  - q*E*C - delta * E;   % * (1 + eta * Fl) * (1- E/Emax) * g*Fi * (1 + eta * Fl)  * C (mass action for base recruitment) Effector cell population rate of change
    dL = pL * Edrive_pos *(1- L/KIL6) - kC6 * L - kB * A * L + kC6*L_initial;   % IL-6 concentration rate of change
    dA = - kCE * A;   % anti-IL-6 concentration rate of change
    dI = - kCI * I;   % anti-CTLA-4 concentration rate of change
    dX = rx * (dL_act/(dL_act + ltox)) * (dE_act/(dE_act + etox)) - rxres * X;  % toxicity approximation

    dx = [dC; dE; dL; dA; dI; dX;];
end