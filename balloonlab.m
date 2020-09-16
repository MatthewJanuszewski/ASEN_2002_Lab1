% Max Martinez
% ASEN 2002 Lab 1
clc
clear
% --------------------------------------------------------------- %

% Constants given to us
Psl = 101300;    % Pressure at sea level [Pascals]
rhosl = 1.2250;  % Air Density at sea level [kg/m^3]
Tsl = 288.16;    % Temperature at sea level [Kelvin]
gsl = 9.80;      % Gravity constant at sea level [m/s^2]
R = 287;         % Ideal gas constant [J/(kg*K)]
Pgage = 10;      % Guage pressure [Pa]
h = 25000;       % Target altitude above sea level [meters]
m = 500;         % Mass of the payload [kg]

% Next three are from Appendix A(pg 853)of the Aerodynamics book
Th = 216.66;    % Temperature at target altitude [Kelvin]
Ph = 2527.3;    % Pressure at target altitude [Pa]
rhoh = 0.040639;% Density of air at target altitude [kg/m^3]

% --------------------------------------------------------------- %

% Calculating the gravititional constant at target altitude
rEarth = 6371000;               % [m] Radius of the earth
g = gsl*(rEarth/(rEarth+h));    % value of g at target altitude

% --------------------------------------------------------------- %

% Material for balloon: Polyester Film
rhoMat = 1390;        % Density of material [kg/m^3]  = 1.39 g/cc
tMat = 0.00000254;    % Thickness of material [meters]
matYS = 27600000;     % Yeild Strength of material [Pa]
FS = 1.5;             % Factor of safety

rshell = (matYS*2*tMat)/(FS*Pgage); % [m]   - Radius of balloon
Vshell = 4*pi*(rshell^2)*tMat;      % [m^3] - Volume of balloon shell
mshell = rhoMat*Vshell;             % [kg]  - Mass of the balloon shell

V = (4*pi*(rshell^3))/3;              % [m^3] - Volume of the balloon

% --------------------------------------------------------------- %

% [T,a,P,rho] = atmosisa(h); % gives us values at target altitude:
% % T = temperature at 25km above sea level
% % a = speed of sound at 25km above sea level
% % P = pressure at 25km above sea level
% % rho = air density at 25km above sea level

% --------------------------------------------------------------- %

% Calculating the density for Helium
MHe = 4.002602;  % Molar Mass of Helium ([u] = [kg/mol])
rhoHe = (Ph*MHe)/(R*Th); % [kg/m^3] - Density of Helium at target altitude

totalWeight = g*(m+mshell+(rhoHe*V));
bouyantForce = rhoh*V*g;
