% ASEN 2002 Thermodynamics Design Laboratory Assignment: Atmospheric
% Satellites
clc
clear
close all
%% Constants
% Flight properties
mass_payload = 500; % [kg]
altitude = 25000; %[m]

% Balloon Material Properties 
% (average of values given on matweb for polyester film)
density_material = 1255; % [kg/m^3]
FS = 1.5; % Factor of Safety
YS = 27.6*10^6; % [Pa]

% Gas properties
gage_pressure = 10;      % [Pa]
R = 8.314;               % [Nm/kmol]
molar_mass_helium = 4.0026e-03; % [kg/mol]
R_helium = 2077.1;       % [J/kgK]

% Radiation Constants
sigmaSB = 5.670*(10^-8); % [J/(K^4m^2s)] - Stephan Boltzman constant
alphasb = .6;            % Absortivity of the sun-balloon system
epsilonb = .8;           % Emissivity of the balloon
alphaeb = epsilonb;      % Absortivity of the earth-balloon system
qsun = 1353;             % [W/m^2] - Solar irradiance
qearth = 237;            % [W/m^2] - Earth irradiance

%% Calculate temperature of the balloon at day and night
temp_day = (((alphasb*qsun)+(alphaeb*qearth))/(epsilonb*sigmaSB*4))^.25;
temp_night = ((alphaeb*qearth)/(epsilonb*sigmaSB*4))^.25;

%% Calculate atmospheric conditions based on the 1976 standard atmosphere at 25km
[temp_25km, speed_of_sound_25km, pressure_25km, density_25km] = ...
    atmoscoesa(altitude, 'None'); % [k, m/s, Pa, kg/m^3]

%% Calculate the density of helium at 25km conditions
density_helium_day = (pressure_25km + 10)/(R_helium * temp_day); % [kg/m^3]
density_helium_night = (pressure_25km + 10)/(R_helium * temp_night); % [kg/m^3]

%% Calculate radius of balloon at night (launch time)
radius_night = nthroot((mass_payload/((((4*pi)/3)*density_25km)-(((4*pi)/3)...
    *density_helium_night)-((4*pi*density_material*gage_pressure*FS)/(2*YS)))), 3); % [m]

%% Calculate mass of balloon material
mass_material = (4*pi*density_material*gage_pressure*FS*(radius_night^3))/(2*YS); % [kg]

%% Calculate # of moles  and mass of helium present at night
volume_night = (4/3)*pi*radius_night^3; % [m^3]
moles_night = ((pressure_25km + gage_pressure)*volume_night)/(R*temp_night); % [mol]
mass_helium_night = moles_night * molar_mass_helium; % [kg]

%% Calcualte the volume of the helium at the daytime temperature
% Note that pressure and # of moles remains constant as temperature
% increases, thus volume must increase. We'll vent gas later




calculated_air_density_night = (density_helium_night*volume_night + ((4*pi*density_material*gage_pressure*FS*radius_night^2)/(2*YS)) + mass_payload)/volume_night
calculated_air_density_day = (density_helium_day*volume_day + ((4*pi*density_material*gage_pressure*FS*radius_day^2)/(2*YS)) + mass_payload)/volume_day







