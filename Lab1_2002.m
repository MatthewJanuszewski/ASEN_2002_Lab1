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
temp_day = (((alphasb*qsun)+(alphaeb*qearth))/(epsilonb*sigmaSB*4))^.25; % [K]
temp_night = ((alphaeb*qearth)/(epsilonb*sigmaSB*4))^.25; % [K]
fprintf(['Calculated temperatures using stefan boltzmann equations:',...
    '\ntemp_day: %.3f K, temp_night: %.3f K\n\n'], temp_day, temp_night);


%% Calculate atmospheric conditions based on the 1976 standard atmosphere at 25km
[temp_25km, speed_of_sound_25km, pressure_25km, density_25km] = ...
    atmoscoesa(altitude, 'None'); % [k, m/s, Pa, kg/m^3]
fprintf(['Atmospheric conditions found using 1976 standard atmosphere:'...
    '\ntemp_25km: %.3f K, speed_of_sound_25km: %.3f m/s, pressure_25km: '...
    '%.3f Pa, density_25km: %.3f kg/m^3\n\n'],... 
    temp_25km, speed_of_sound_25km, pressure_25km, density_25km);


%% Calculate the density of helium at 25km conditions
density_helium_day = (pressure_25km + 10)/(R_helium * temp_day); % [kg/m^3]
density_helium_night = (pressure_25km + 10)/(R_helium * temp_night); % [kg/m^3]
fprintf(['Calculated density of helium using the ideal gas law\n'...
    'density_helium_day: %.6f kg/m^3, density_helium_night: '...
    '%.6f kg/m^3\n\n'], density_helium_day, density_helium_night);


%% Calculate radius of balloon at night (launch time)
radius_night = nthroot((mass_payload/((((4*pi)/3)*density_25km)-(((4*pi)/3)...
    *density_helium_night)-((4*pi*density_material*gage_pressure*FS)/(2*YS)))), 3); % [m]
fprintf(['Radius at night calculated using the formula derrived in design specifications:\n'...
    'radius_night: %.3f m\n\n'], radius_night);

%% Calculate mass of balloon material
mass_material = (4*pi*density_material*gage_pressure*FS*(radius_night^3))/(2*YS); % [kg]
fprintf('Mass of the balloon material:\nmass_material: %.3f kg\n\n', mass_material);


%% Calculate # of moles  and mass of helium present at night
volume_helium_night = (4/3)*pi*radius_night^3; % [m^3]
moles_helium_night = ((pressure_25km + gage_pressure)*volume_helium_night)/(R*temp_night); % [mol]
mass_helium_night = moles_helium_night * molar_mass_helium; % [kg]
fprintf(['Extent of the helium gas at night:\n',...
    'volume_helium_night: %.3f m^3, moles_helium_night: %.3f mol, mass_helium_night %.3f kg\n\n'],...
    volume_helium_night, moles_helium_night, mass_helium_night);


%% Calculate how much helium we need to vent
mass_helium_day = (mass_payload+mass_material)/((density_25km/density_helium_day)-1);
mass_delta = mass_helium_night - mass_helium_day;
fprintf(['Mass of daytime balloon after venting:\n',...
    'mass_helium_day: %.3f kg, mass_delta: %.3f kg\n\n'],...
    mass_helium_day, mass_delta);


%% Calculate the volume of the daytime balloon
volume_helium_day = mass_helium_day/density_helium_day;
volume_delta = volume_helium_night - volume_helium_day;
fprintf(['Volume of daytime balloon after venting:\n',...
    'volume_helium_day: %.3f m^3, volume_delta: %.3f m^3\n\n'],...
    volume_helium_day, volume_helium_night);

