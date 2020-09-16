% ASEN 2002 Thermodynamics Design Laboratory Assignment: Atmospheric
% Satellites

clc
clear all

% Variables

mass_payload = 500; % [kg]
altitude = 25000; %[m] (+/- 1000m)
FS = 1.5; % Factor of Safety
Gage_Pressure = 10; % [Pa] (diffence of pressure of balloon and atmosphere
K = 8.31432; % [Nm/kmol] Universal Gas Constant
thickness = 2.54*10^-6; % [m] thickness of balloon material
density_material = 0;
YS = 27.6*10^6; % [PA] Yield Strength of Material

% Equations:

Radius_Balloon = (2*YS*thickness)/(Gage_Pressure*FS)

Volume_Balloon = (4*pi*Radius_Balloon^3)/3