% ASEN 2002 Thermodynamics Design Laboratory Assignment: Atmospheric
% Satellites

clc
clear all

% Variables

mass_payload = 500; % [kg]
altitude = 25000; %[m] (+/- 1000m)
FS = 1.5; % Factor of Safety
Gage_Pressure = 10; % [Pa] (diffence of pressure of balloon and atmosphere

R = 8.31432; % [Nm/kmol] Universal Gas Constant
R_air = R/0.02897; % [J/kgK]
R_helium = R/(4.003*10^-3); % [J/kgK]

thickness = 2.54*10^-6; % [m] thickness of balloon material
density_material = 0;
YS = 27.6*10^6; % [PA] Yield Strength of Material

sigma_SB = 5.67*10^-8; % [J/k^4m^2s] SB constant
alpha_SB = 0.6; % absorptivity of sun-balloon system
epsilon_b = 0.8; % emissivity of balloon
alpha_eb = 0.8; % absorptivity of Earth-balloon system
q_sun = 1353; % [W/m^2] Solar constant/irradiance
q_earth = 237; % [W/m^2]

% Equations:

Radius_Balloon = (2*YS*thickness)/(Gage_Pressure*FS)

Volume_Night = (4*pi*Radius_Balloon^3)/3

T_Night = (alpha_eb*q_earth/(4*epsilon_b*sigma_SB))^0.25

T_Day = (((alpha_SB*q_sun)+(alpha_eb*q_earth))/(4*epsilon_b*sigma_SB))^0.25

Volume_Day = T_Day*Volume_Night/T_Night
