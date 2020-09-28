function P = idealGasSolveForP(v, n, T)
    % Solves ideal gas equation for P given V, n , T
    % Load Constants
    
    % Solve ideal gas law for pressure
    P = (n*R*T)/v;
    return
    