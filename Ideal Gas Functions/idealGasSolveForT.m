function T = idealGasSolveForT(v, n, P)
    % Solves ideal gas equation for T given v, n , P
    % Load Constants
    load('constants.mat', 'universalGasConstant');
    R = universalGasConstant;
    
    % Solve ideal gas law for pressure
    P = (n*R*T)/v;
    return
    