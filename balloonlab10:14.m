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
Pgage = 10;      % Gage pressure [Pa]
h = 25000;       % Target altitude above sea level [meters]
m = 500;         % Mass of the payload [kg]

% % Next three are from Appendix A(pg 853)of the Aerodynamics book
% Th = 216.66;     % Temperature at target altitude [Kelvin]
% Ph = 2527.3;     % Pressure at target altitude [Pa]
% rhoh = 0.040639; % Density of air at target altitude [kg/m^3]

% --------------------------------------------------------------- %

% Calculating the gravititional constant at target altitude
rEarth = 6371000;               % [m] Radius of the earth
g = gsl*(rEarth/(rEarth+h));    % value of g at target altitude

% --------------------------------------------------------------- %

% Material for balloon: Polyester Film
rhoMat = 1255;        % Density of material [kg/m^3]  = 1.39 g/cc
tMat = 0.00000254;    % Thickness of material [meters]
matYS = 27600000;     % Yeild Strength of material [Pa]
FS = 1.5;             % Factor of safety

rshell = (matYS*2*tMat)/(FS*Pgage); % [m]   - Radius of balloon
Vshell = 4*pi*(rshell^2)*tMat;      % [m^3] - Volume of balloon shell
mshell = rhoMat*Vshell;             % [kg]  - Mass of the balloon shell

V = (4*pi*(rshell^3))/3;            % [m^3] - Volume of the balloon
SA = 4*pi*(rshell^2);               % [m^2] - Surface area of balloon

% --------------------------------------------------------------- %

[T,a,P,rho] = atmoscoesa(h); % gives us values at target altitude:
% T = temperature at 25km above sea level
% a = speed of sound at 25km above sea level
% P = pressure at 25km above sea level
% rho = air density at 25km above sea level

% --------------------------------------------------------------- %

% Calculating the density for Helium
MHe = 4.002602;  % Molar Mass of Helium ([u] = [kg/mol])
Rhe = 2076.9;    % Specific gas constant for helium
rhoHe = (P)/(Rhe*T); % [kg/m^3] - Density of Helium at target altitude
mHe = rhoHe*V;   % Mass of the helium
n = mHe/MHe;     % Ammount of helium (mol)

% Thermal Radiation stuff
% --------------------------------------------------------------- %

% Constants:
sigmaSB = 5.670*(10^-8); % [J/(K^4m^2s)] - Stephan Boltzman constant
alphasb = .6;            % Absortivity of the sun-balloon system
epsilonb = .8;           % Emissivity of the balloon
alphaeb = epsilonb;      % Absortivity of the earth-balloon system
qsun = 1353;             % [W/m^2] - Solar irradiance
qearth = 237;            % [W/m^2] - Earth irradiance

% Equations:
% Qearth = alphaeb*qearth*pi*rshell^2; % Rate of heat transfer from Earth
% Qsolar = alphaeb*qsun*pi*rshell^2;   % Rate of heat transfer from the sun
% Qballoon = epsilonb*sigmaSB*4*pi*rshell^2*Tballoon^4
% Qincident = q*A;
% Qabsorbed = alpha*Qincident
% Qblackbody = epsilon*sigmaSB*As*Ts^4

Tnight = ((alphaeb*qearth)/(epsilonb*sigmaSB*4))^.25;
rhoHeNight = (P+Pgage)/(Rhe*Tnight);
radiusNight = (500/((((4*pi)/3)*(rho-rhoHeNight))-rhoMat*((4*pi*FS*Pgage)/(2*matYS))))^(1/3);

mmat = rhoMat*4*pi*(radiusNight^2)*tMat;
Vnight = (4*pi*(radiusNight^3))/3;
mHeNight = Vnight*rhoHeNight;
rhoAirNight = (m+mmat+mHeNight)/Vnight;
% VdispNight = (m+1255*4*pi*((radiusNight)^2)*(tMat))/(rho-rhoHeNight);
% mHeNight = VdispNight*rhoHeNight;
% mmat = rhoMat*4*pi*(radiusNight^2)*tMat;
% rhoAirNight = (m+mmat+mHeNight)/VdispNight;

Tday = (((alphasb*qsun)+(alphaeb*qearth))/(epsilonb*sigmaSB*4))^.25;
Vday = (Vnight/Tnight)*Tday;
radiusDay = ((3/4*pi)*Vday)^(1/3);
rhoHeDay = (P+Pgage)/(Rhe*Tday);
mHeDay = Vday*rhoHeDay;
rhoAirDay = (m+mmat+mHeDay)/Vday;

nday = ((P+Pgage)*Vnight)/(Rhe*Tday)
nnight = ((P+Pgage)*Vnight)/(Rhe*Tnight)

vent = nnight-nday
massHeDay = Vnight*rhoHeDay
massHeNight = Vnight*rhoHeNight

% rDay = (500/((((4*pi)/3)*(rho-rhoHeDay))-rhoMat*((4*pi*FS*Pgage)/(2*matYS))))^(1/3);

% rhoAirDay = (m+(rhoMat*4*pi*(radiusDay^2)*tMat)+(rhoHeDay*VdispDay))/(VdispDay);
% rhoAirNight = (m+(rhoMat*4*pi*(radiusNight^2)*tMat)+(rhoHeNight*VdispNight))/(VdispNight);

% rhoHeDay = (P+Pgage)/(Rhe*Tday);
% radiusDay = (500/((((4*pi)/3)*(rho-rhoHeDay))-rhoMat*((4*pi*FS*Pgage)/(2*matYS))))^(1/3);
% Vday = (4*pi*(radiusDay^3))/3;
% VdispDay = (m+1255*4*pi*((radiusDay)^2)*(tMat))/(rho-rhoHeDay);
% mHeDay = VdispDay*rhoHeDay;
% mmat2 = rhoMat*4*pi*(radiusDay^2)*tMat;
% rhoAirDay = (m+mmat+mHeDay)/VdispDay;

[Thigh,ahigh,Phigh,rhohigh] = atmoscoesa(h+1000); % 26 km
[Tlow,alow,Plow,rholow] = atmoscoesa(h-1000);     % 24 km
