function params = IL6ParamUnpack(theta)
    params.tEnd = theta.tEnd.val;
    
    params.useIDosing = theta.useIDosing.val;
    params.doseInterval_I = theta.doseInterval_I.val;
    params.doseI = theta.doseI.val;
    params.doseStart_I = theta.doseStart_I.val;
    params.nDoses_I = round(theta.nDoses_I.val);

    params.useADosing = theta.useADosing.val;
    params.doseInterval_A = theta.doseInterval_A.val;
    params.doseA = theta.doseA.val;
    params.doseStart_A = theta.doseStart_A.val;
    params.nDoses_A = round(theta.nDoses_A.val);
    
    params.C_initial = theta.C_initial.val;
    params.E_initial = theta.E_initial.val;
    params.L_initial = theta.L_initial.val;
    params.A_initial = theta.A_initial.val;
    params.I_initial = theta.I_initial.val;
    params.X_initial = theta.X_initial.val;
    
    params.rT = theta.rT.val;
    params.kT = theta.kT.val;
    params.kE = theta.kE.val;
    params.Emax = theta.Emax.val;
    params.kA = theta.kA.val;
    params.eta = theta.eta.val;
    params.L0 = theta.L0.val;
    params.L1 = theta.L1.val;
    params.n = theta.n.val;
    params.m = theta.m.val;
    params.beta = theta.beta.val;
    params.delta = theta.delta.val;
    params.pL = theta.pL.val;
    params.phi = theta.phi.val;
    params.kC6 = theta.kC6.val;
    params.kB = theta.kB.val;
    params.KIL6 = theta.KIL6.val;
    params.kCE = theta.kCE.val;
    params.kCI = theta.kCI.val;
    params.kICI = theta.kICI.val;
    params.s = theta.s.val;
    params.g = theta.g.val;
    params.p = theta.p.val;
    params.q = theta.q.val;
    params.rx = theta.rx.val;
    params.rxres = theta.rxres.val;
    params.ltox = theta.ltox.val;
    params.etox = theta.etox.val;
end