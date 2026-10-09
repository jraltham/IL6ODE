function [tAll, xAll] = IL6RunSim(params)

    % Initial conditions
    x0 = [
        params.C_initial
        params.E_initial
        params.L_initial
        params.A_initial
        params.I_initial
        params.X_initial
    ];

    t0 = 0;
    tAll = [];
    xAll = [];

    opts = odeset('RelTol',1e-6,'AbsTol',1e-9, 'MaxStep', 50);

    % ==============================
    % Build dosing schedules
    % ==============================

    allDoseTimes = [];

    if params.useIDosing
        doseTimes_I = params.doseStart_I + ...
                      params.doseInterval_I * (0:params.nDoses_I-1);
        allDoseTimes = [allDoseTimes doseTimes_I];
    else
        doseTimes_I = [];
    end

    if params.useADosing
        doseTimes_A = params.doseStart_A + ...
                      params.doseInterval_A * (0:params.nDoses_A-1);
        allDoseTimes = [allDoseTimes doseTimes_A];
    else
        doseTimes_A = [];
    end

    % Sort all dosing times
    allDoseTimes = sort(allDoseTimes);

    % ==============================
    % Piecewise integration
    % ==============================

    for k = 1:length(allDoseTimes)

        [tSol, xSol] = ode15s(@(t,x) IL6DosingODE(t,x,params), ...
                              [t0 allDoseTimes(k)], x0, opts);

        tAll = [tAll; tSol];
        xAll = [xAll; xSol];

        % Update state at dosing time
        x0 = xSol(end,:)';

        % --- Apply CTLA-4 bolus ---
        if ismember(allDoseTimes(k), doseTimes_I)
            x0(5) = x0(5) + params.doseI;
        end

        % --- Apply Anti-IL-6 bolus ---
        if ismember(allDoseTimes(k), doseTimes_A)
            x0(4) = x0(4) + params.doseA;
        end

        t0 = allDoseTimes(k);
    end

    % ==============================
    % Final segment
    % ==============================

    [tSol, xSol] = ode15s(@(t,x) IL6DosingODE(t,x,params), ...
                          [t0 params.tEnd], x0, opts);

    tAll = [tAll; tSol];
    xAll = [xAll; xSol];
end
