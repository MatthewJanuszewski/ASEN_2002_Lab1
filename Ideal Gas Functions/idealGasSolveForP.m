function P = idealGasSolveForP(v, n, T)
    % Solves ideal gas equation for P given v, n , T
    % Load Constants
    load('constants.mat', 'universalGasConstant');
    R = universalGasConstant;
    
    % Solve ideal gas law for pressure
    P = (n*R*T)/v;
    return
    