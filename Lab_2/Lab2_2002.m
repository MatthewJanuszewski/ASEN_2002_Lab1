% ASEN 2002 Aerodynamics Design Laboratory Assignment: Wind Tunnels
clc
clear
close all
warning('off')

%% Part 1: Airspeed Model
% Note: For the pitot static + water maonometer readings I am using data 
% from group 4 (1st row in water monometer readings file), and for
% venturi + water manometer readings I am using data from group
% 6 (8th row in water manometer readings file)

% For atmospheric temperature and pressure readings I am using data from 
% VenturiTubeS301_2

% I am using the assumption that the difference in atmospheric temperature
% and pressure between when S301_2 and S301_1 were measured and when the 
% actual wind tunnel tests were conducted is negliigable compared to the
% accuracy of the sensors involved.

%% Water Manometer 
% Import data from water manometer readings (Responses).xlsx and get rows 1
% and 4
water_manometer_data = readtable('/Users/matthewj/Documents/GitHub/ASEN_2002_Lab1/Lab_2/selected_data/monometer_voltage_observatons.xlsx','VariableNamingRule','preserve');
water_manometer_data_venturi = table2array(water_manometer_data(8,4:13));
water_manometer_data_pitot_static = table2array(water_manometer_data(1,4:13));


% Venturi tube + water manometer array, column 1 is voltage, column 2 is
% manometer height reading. Height reading (pressure) represents static
% pressure differential between test section and settling chamber (before
% test section)
venturi_water = zeros(5,2); % [volts][in H2O]
venturi_water(:,1) = water_manometer_data_venturi(1:2:9);
venturi_water(:,2) = water_manometer_data_venturi(2:2:10);


% Pitot-static probe + water manometer array, column 1 is voltage, column 2 is
% manometer height reading. Height reading (pressure) represents pressure
% differential between total pressure and static pressure in test section,
% which is equal to the dynamic pressure in the test section
pitot_static_water = zeros(5,2); % [volts][in H2O]
pitot_static_water(:,1) = water_manometer_data_pitot_static(1:2:9);
pitot_static_water(:,2) = water_manometer_data_pitot_static(2:2:10);


% Convert inches H2O into Pa for venturi_water and pitot_static_water
% arrays
venturi_water(:,2) = venturi_water(:,2)*249.0889; % [volts][Pa]
pitot_static_water(:,2) = pitot_static_water(:,2)*249.0889; % [volts][Pa]


% Load atmospheric (ouside the wind tunnel) pressure and temperature data
% and average the readings
atmospheric_data = importfile('/Users/matthewj/Documents/GitHub/ASEN_2002_Lab1/Lab_2/selected_data/VenturiTube_S301_2.csv');
atmospheric_data(1,:) = [];
atmospheric_data = table2array(atmospheric_data);
atm_pressure = mean(atmospheric_data(:,1)); % [Pa]
atm_temp = mean(atmospheric_data(:,2)); % [K]




