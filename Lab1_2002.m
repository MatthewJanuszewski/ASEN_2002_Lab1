% ASEN 2002 Thermodynamics Design Laboratory Assignment: Atmospheric
% Satellites

clc
clear
close all
%% Constants
mass_payload = 500; % [kg]
altitude = 25000; %[m] (+/- 1000m)
FS = 1.5; % Factor of Safety
Gage_Pressure = 10; % [Pa] (diffence of pressure of balloon and atmosphere
R = 8.314; % [Nm/kmol] Universal Gas Constant
thickness = 2.54*10^-6; % [m] thickness of balloon material
density_material = 0;
YS = 27.6*10^6; % [Pa] Yield Strength of Material
molar_mass_helium = 4.0026e-03; % [kg/mol]

%% Calculate atmospheric conditions based on the 1976 standard atmosphere at 25km
[temp_25km, speed_of_sound_25km, pressure_25km, density_25km] = atmoscoesa(altitude, 'None'); % [k, m/s, Pa, kg/m^3]

%% Calculate volume and mass of gas required at 25km
volume_helium = mass_payload/density_25km; % [m^3]
moles_helium = (pressure_25km * volume_helium) / (R * temp_25km); % [mol]
mass_helium = moles_helium * molar_mass_helium; % [kg]

%% Find radius of balloon
radius = nthroot(((3/(4*pi))*volume_helium), 3);

