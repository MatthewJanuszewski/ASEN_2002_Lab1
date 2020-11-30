function dydt = getodes(t,y,P_g,V_water_i,C_D,theta)

% function returns change in volume of air within the bottle vs time before
% all water is expelled from the bottle, mass flow rate of air after water 
% is exhausted from the bottle, velocity in x direction, velocity in z 
% direction, acceleration in x direction, and acceleration in z direction

% Thus it will be the case that...
% dydt(:,1) = dvdt = change in volume of air vs time
% dydt(:,2) = dmdt = change in mass of rocket vs time
% dydt(:,3) = dxdt = horizontal velocity over time
% dydt(:,4) = dzdt = vertical velcity over time
% dydt(:,5) = ddxdtdt = horizontal acceleration over time
% dydt(:,6) = ddzdtdt = vertical accleration over time
% dydt(:,7) = dpdt = change in pressure over time 

global Thrust Time 

%% Constants, Parameters, and Variables

% Constants
g = 9.81; % acceleration of gravity [m/s^2]
C_d = 0.8; % discharge coefficient
rho_air = 0.961; % density of air [kg/m^3]
rho_water = 1000; % density of water [kg/m^3]
gamma = 1.4; % specific heat ratio for air
R = 287; % gas constant of air [J/kgK]

% Parameters
V_b = 0.002; % volume of the bottle [m^3]
M_b = 0.15; % mass of empty bottle rocket [kg]
P_a = 12.1; % atmospheric (ambient) pressure [psi]
P_a = P_a*6894.76; % [Pa] convert to SI
D_t = 2.1; % diameter of the throat [cm]
D_t = D_t/100; % [m] convert to SI
D_b = 10.5; % diamerer of the bottle [cm]
D_b = D_b/100; % [m] convert to SI
T_air_i = 300; % initial temperature of air [K]
v0 = 0.0; % initial velocity of rocket [m/s]
x0 = 0.0; % initial horizontal distance [m]
z0 = 0.25; % initial vertical distance [m]
l_s = 0.5; % length of test stand [m]

% cross sectional areas calculated
A_t = pi*(D_t/2)^2; % throat [m^2]
A_b = pi*(D_b/2)^2; % bottle [m^2]

% additional initial conditions/parameters
V_air_i = V_b - V_water_i; % initial volume of air [m^3]
P_air_i = P_a + P_g; % total initial pressure of air inside bottle
M_tot_i = M_b + (rho_water * V_water_i) + (P_air_i/(R*T_air_i))*V_air_i;

% dzdt(1,1) = change in volume of air
% dzdt(1,2) = change in mass of rocket
% dzdt(1,3) = change in pressure over time
% dzdt(1,4) = horizontal velocity over time
% dzdt(1,5) = vertical velcity over time
% dzdt(1,6) = horizontal acceleration over time
% dzdt(1,7) = vertical accleration over time

%% Condition vector assign
Vol = y(1);
M = y(2);
P = y(3);
x = y(4);
z = y(5);
dxdt = y(6);
dzdt = y(7);

if(Vol < V_b)
    %% before water is exhausted
    
    V = sqrt(dxdt^2 + dzdt^2); % velocity magnitude
    if(sqrt(x^2 + z^2) <= l_s)
       % heading components
       h_x = cos(theta);
       h_z = sin(theta);
    else
       % heading components
       h_x = dxdt/V;
       h_z = dzdt/V;
    end
    
    % change in mass of rocket
    dmdt = -C_d*A_t*sqrt(2*rho_water*(P-P_a));
    
    % change in volume
    dvdt = C_d*A_t*sqrt((2/rho_water)*(P_air_i*((V_air_i/Vol)^gamma)-P_a));
    
    % change in pressure
    dpdt = P_air_i*V_air_i^(gamma)*(-gamma*Vol^(-gamma-1))*dvdt;
    
    % exit velocity
    V_e  = sqrt(2*(P-P_a)./rho_water);
    
    % thrust function F(t)
    F = (2*C_d*A_t)*(P-P_a);
    D = 0.5*rho_air*(V^2)*C_D*A_b;
    F_net_x = (F-D)*h_x;
    F_net_z = (F-D)*h_z-(M*g);
    
elseif(P >= P_a)   
    %% after water is exhausted
    
    V = sqrt(dxdt^2 + dzdt^2); % velocity magnitude
    if(sqrt(x.^2+y.^2) <= l_s)
        % heading components
        h_x = cos(theta);
        h_z = sin(theta);
    else
        % heading components
        h_x = dxdt/V;
        h_z = dzdt/V;
    end
    
    P_end = P_air_i*(V_air_i/V_b)^gamma; % pressure after all water expelled
    T_end = T_air_i*(V_air_i/V_b)^(gamma-1); % temperature after all water expelled
    M_air_i = (P_air_i/(R*T_air_i))*V_air_i; % mass of air in bottle at beginning of phase 2
    
    % constant volume
    dvdt = 0;
    
    % changing density of air
    rho_ch = (M-M_b)/V_b;
    
    % Temperature function T(t)
    T = P/(rho_ch*R);

    % critical pressure
    P_cr = P*(2/(gamma+1))^(gamma/(gamma-1));

    if(P_cr > P_a)
        % choked flow
        M_e = 1; % exit mach number
        P_e = P_cr; % exit pressure
        T_e = T*(2/(gamma+1)); % exit temperature
        rho_e = P_e/(R*T_e); % exit density
        V_e = sqrt((gamma*R)*T_e); % exit velocity
    else
        % not choked flow
        M_e = sqrt(abs((2/(gamma-1))*((P/P_a)^((gamma-1)/gamma)-1)));
        P_e = P_a;
        T_e = T/(1+(((gamma-1)/2)*M_e^2));
        rho_e = P_a/(R*T_e);
        V_e = M_e*sqrt((gamma*R)*T_e);
    end
    
    % change in mass of rocket 
    dmdt = -(C_d*rho_e*A_t)*V_e;
   
    % thrust function F(t)
    F = (abs(dmdt)*V_e) + ((P_e-P_a)*A_t);
    D = 0.5*rho_air*(V^2)*C_D*A_b;
    F_net_x = (F-D)*h_x;
    F_net_z = (F-D)*h_z+(M*g);
    
    % change in pressure
    dpdt = P_end*M_air_i^(-gamma)*(gamma*(M-M_b)^(gamma-1))*dmdt;
   
else
   %% Ballistic Phase
   
   V = sqrt(dxdt^2 + dzdt^2); % velocity magnitude
   
    % heading components
        h_x = dxdt/V;
        h_z = dzdt/V;
   
    % constant mass
    dmdt = 0;
        
    % constant volume
    dvdt = 0;
    
    % constant pressure
    dpdt = 0;
    
    % drag
    D = 0.5*rho_air*(V^2)*C_D*A_b;
    F = 0;
    M = M_b;
    F_net_x = -D*h_x;
    F_net_z = (-D*h_z)-(M_b*g);
   
end

%% Thrust, Time, Acceleration, and all differential equations (dydt)
Thrust = [Thrust; real(F)];
Time = [Time; t];

% acceleration
ddxdtdt = F_net_x/M;
ddzdtdt = F_net_z/M;

dydt = [dvdt;dmdt;dpdt;dxdt;dzdt;ddxdtdt;ddzdtdt];

end