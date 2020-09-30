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
gage_pressure = 10; % [Pa]
R = 8.314; % [Nm/kmol]
molar_mass_helium = 4.0026e-03; % [kg/mol]

%% Calculate atmospheric conditions based on the 1976 standard atmosphere at 25km
[temp_25km, speed_of_sound_25km, pressure_25km, density_25km] = atmoscoesa(altitude, 'None'); % [k, m/s, Pa, kg/m^3]

%% Calculate volume and mass of gas required at 25km
volume_helium = mass_payload/density_25km; % [m^3]
moles_helium = (pressure_25km * volume_helium) / (R * temp_25km); % [mol]
mass_helium = moles_helium * molar_mass_helium; % [kg]
density_helium = mass_helium/volume_helium; % [kg/m^3]

%% Find radius of balloon
radius = nthroot((mass_payload/((((4*pi)/3)*density_25km)-(((4*pi)/3)*density_helium)-(4*pi*density_material*((gage_pressure*FS)/(2*YS))))), 3); % [m]











