% ASEN 2002 Lab 2 Part 2
% Boundary Layer thickness

% Housekeeping
clc
clear all
close all

R = 287; % Gas Constant for Air

% Part 1
% Read in data

Port_1A = readtable('BoundaryLayer_S301_1.csv');
Port_1B = readtable('BoundaryLayer_S301_2.csv');
Port_1A = table2array(Port_1A);
Port_1B = table2array(Port_1B);
Port_1_P_atm = mean([Port_1A(:,1); Port_1B(:,1)]); % P_atm average [Pa]
Port_1_T_atm = mean([Port_1A(:,2); Port_1B(:,2)]); % T_atm average [K]
rho = (Port_1_P_atm/(R*Port_1_T_atm)); % density [kg/m^3]
Port_1_Aux_Diff_Pressure = [Port_1A(:,4); Port_1B(:,4)]; % [Pa]
Port_1_Velocity_Probe = sqrt(2.*Port_1_Aux_Diff_Pressure/rho); % [m/s]
Port_1_Free_Stream_Pressure = mean([Port_1A(5501:6000,3); Port_1B(5501:6000,3)]); % [Pa]
Port_1_Free_Stream_Velocity = sqrt(2*Port_1_Free_Stream_Pressure/rho); % [m/s]
Port_1_Boundary_Layer_Velocity = 0.95*Port_1_Free_Stream_Velocity; % [m/s]
Port_1_y_axis = abs([Port_1B(:,6); Port_1A(:,6)]); % vertical height[mm]

x = (Port_1_Velocity_Probe > Port_1_Boundary_Layer_Velocity);
Port_1_Boundary_Layer_Height = Port_1_y_axis(x);
Port_1_Boundary_Layer_Height = Port_1_Boundary_Layer_Height(1); % [mm]


Port_2A = readtable('BoundaryLayer_S302_1.csv');
Port_2B = readtable('BoundaryLayer_S302_2.csv');
Port_2A = table2array(Port_2A);
Port_2B = table2array(Port_2B);
Port_2_P_atm = mean([Port_2A(:,1); Port_2B(:,1)]); % P_atm average [Pa]
Port_2_T_atm = mean([Port_2A(:,2); Port_2B(:,2)]); % T_atm average [K]
rho = (Port_2_P_atm/(R*Port_2_T_atm)); % density [kg/m^3]
Port_2_Aux_Diff_Pressure = [Port_2A(:,4); Port_2B(:,4)]; % [Pa]
Port_2_Velocity_Probe = sqrt(2.*Port_2_Aux_Diff_Pressure/rho); % [m/s]
Port_2_Free_Stream_Pressure = mean([Port_2A(5501:6000,3); Port_2B(5501:6000,3)]); % [Pa]
Port_2_Free_Stream_Velocity = sqrt(2*Port_2_Free_Stream_Pressure/rho); % [m/s]
Port_2_Boundary_Layer_Velocity = 0.95*Port_2_Free_Stream_Velocity; % [m/s]
Port_2_y_axis = abs([Port_2B(:,6); Port_2A(:,6)]); % vertical height[mm]

x = (Port_2_Velocity_Probe > Port_2_Boundary_Layer_Velocity);
Port_2_Boundary_Layer_Height = Port_2_y_axis(x);
Port_2_Boundary_Layer_Height = Port_2_Boundary_Layer_Height(1); % [mm]


Port_3A = readtable('BoundaryLayer_S303_1.csv');
Port_3B = readtable('BoundaryLayer_S303_2.csv');
Port_3A = table2array(Port_3A);
Port_3B = table2array(Port_3B);
Port_3_P_atm = mean([Port_3A(:,1); Port_3B(:,1)]); % P_atm average [Pa]
Port_3_T_atm = mean([Port_3A(:,2); Port_3B(:,2)]); % T_atm average [K]
rho = (Port_3_P_atm/(R*Port_3_T_atm)); % density [kg/m^3]
Port_3_Aux_Diff_Pressure = [Port_3A(:,4); Port_3B(:,4)]; % [Pa]
Port_3_Velocity_Probe = sqrt(2.*Port_3_Aux_Diff_Pressure/rho); % [m/s]
Port_3_Free_Stream_Pressure = mean([Port_3A(5501:6000,3); Port_3B(5501:6000,3)]); % [Pa]
Port_3_Free_Stream_Velocity = sqrt(2*Port_3_Free_Stream_Pressure/rho); % [m/s]
Port_3_Boundary_Layer_Velocity = 0.95*Port_3_Free_Stream_Velocity; % [m/s]
Port_3_y_axis = abs([Port_3B(:,6); Port_3A(:,6)]); % vertical height[mm]

x = (Port_3_Velocity_Probe > Port_3_Boundary_Layer_Velocity);
Port_3_Boundary_Layer_Height = Port_3_y_axis(x);
Port_3_Boundary_Layer_Height = Port_3_Boundary_Layer_Height(1); % [mm]



Port_4A = readtable('BoundaryLayer_S301_3.csv');
Port_4B = readtable('BoundaryLayer_S301_4.csv');
Port_4A = table2array(Port_4A);
Port_4B = table2array(Port_4B);
Port_4_P_atm = mean([Port_4A(:,1); Port_4B(:,1)]); % P_atm average [Pa]
Port_4_T_atm = mean([Port_4A(:,2); Port_4B(:,2)]); % T_atm average [K]
rho = (Port_4_P_atm/(R*Port_4_T_atm)); % density [kg/m^3]
Port_4_Aux_Diff_Pressure = [Port_4A(:,4); Port_4B(:,4)]; % [Pa]
Port_4_Velocity_Probe = sqrt(2.*Port_4_Aux_Diff_Pressure/rho); % [m/s]
Port_4_Free_Stream_Pressure = mean([Port_4A(5501:6000,3); Port_4B(5501:6000,3)]); % [Pa]
Port_4_Free_Stream_Velocity = sqrt(2*Port_4_Free_Stream_Pressure/rho); % [m/s]
Port_4_Boundary_Layer_Velocity = 0.95*Port_4_Free_Stream_Velocity; % [m/s]
Port_4_y_axis = abs([Port_4B(:,6); Port_4A(:,6)]); % vertical height[mm]

x = (Port_4_Velocity_Probe > Port_4_Boundary_Layer_Velocity);
Port_4_Boundary_Layer_Height = Port_4_y_axis(x);
Port_4_Boundary_Layer_Height = Port_4_Boundary_Layer_Height(1); % [mm]


Port_5A = readtable('BoundaryLayer_S302_3.csv');
Port_5B = readtable('BoundaryLayer_S302_4.csv');
Port_5A = table2array(Port_5A);
Port_5B = table2array(Port_5B);
Port_5_P_atm = mean([Port_5A(:,1); Port_5B(:,1)]); % P_atm average [Pa]
Port_5_T_atm = mean([Port_5A(:,2); Port_5B(:,2)]); % T_atm average [K]
rho = (Port_5_P_atm/(R*Port_5_T_atm)); % density [kg/m^3]
Port_5_Aux_Diff_Pressure = [Port_5A(:,4); Port_5B(:,4)]; % [Pa]
Port_5_Velocity_Probe = sqrt(2.*Port_5_Aux_Diff_Pressure/rho); % [m/s]
Port_5_Free_Stream_Pressure = mean([Port_5A(5501:6000,3); Port_5B(5501:6000,3)]); % [Pa]
Port_5_Free_Stream_Velocity = sqrt(2*Port_5_Free_Stream_Pressure/rho); % [m/s]
Port_5_Boundary_Layer_Velocity = 0.95*Port_5_Free_Stream_Velocity; % [m/s]
Port_5_y_axis = abs([Port_5B(:,6); Port_5A(:,6)]); % vertical height[mm]

x = (Port_5_Velocity_Probe > Port_5_Boundary_Layer_Velocity);
Port_5_Boundary_Layer_Height = Port_5_y_axis(x);
Port_5_Boundary_Layer_Height = Port_5_Boundary_Layer_Height(1); % [mm]


Port_6A = readtable('BoundaryLayer_S301_6.csv');
Port_6B = readtable('BoundaryLayer_S301_9.csv');
Port_6A = table2array(Port_6A);
Port_6B = table2array(Port_6B);
Port_6_P_atm = mean([Port_6A(:,1); Port_6B(:,1)]); % P_atm average [Pa]
Port_6_T_atm = mean([Port_6A(:,2); Port_6B(:,2)]); % T_atm average [K]
rho = (Port_6_P_atm/(R*Port_6_T_atm)); % density [kg/m^3]
Port_6_Aux_Diff_Pressure = [Port_6A(:,4); Port_6B(:,4)]; % [Pa]
Port_6_Velocity_Probe = sqrt(2.*Port_6_Aux_Diff_Pressure/rho); % [m/s]
Port_6_Free_Stream_Pressure = mean([Port_6A(5501:6000,3); Port_6B(5501:6000,3)]); % [Pa]
Port_6_Free_Stream_Velocity = sqrt(2*Port_6_Free_Stream_Pressure/rho); % [m/s]
Port_6_Boundary_Layer_Velocity = 0.95*Port_6_Free_Stream_Velocity; % [m/s]
Port_6_y_axis = abs([Port_6A(:,6); Port_6B(:,6)]); % vertical height[mm]

x = (Port_6_Velocity_Probe > Port_6_Boundary_Layer_Velocity);
Port_6_Boundary_Layer_Height = Port_6_y_axis(x);
Port_6_Boundary_Layer_Height = Port_6_Boundary_Layer_Height(1); % [mm]

Total_Boundary_Layer_Heights = [Port_1_Boundary_Layer_Height Port_2_Boundary_Layer_Height Port_3_Boundary_Layer_Height Port_4_Boundary_Layer_Height Port_5_Boundary_Layer_Height Port_6_Boundary_Layer_Height ];
plot(Total_Boundary_Layer_Heights)