% ASEN 2002 Lab 2 Wind Tunnel
% Section 1: Measurement of Airspeed

% Housekeeping
clc
clear all
close all

% Part 1
% Read in data

A = readtable('VelocityVoltage_S301_1.csv');
A = table2array(A);
A_Pressure_atm = A(1:end,1);
A_Temperature_atm = A(1:end,2);
A_Airspeed_Diff_Pressure = A(1:end,3);
A_Aux_Diff_Pressure = A(1:end,4);
A_Voltage = A(1:end,7);

B = readtable('VelocityVoltage_S301_3.csv');
B = table2array(B);
B_Pressure_atm = B(1:end,1);
B_Temperature_atm = B(1:end,2);
B_Airspeed_Diff_Pressure = B(1:end,3);
B_Aux_Diff_Pressure = B(1:end,4);
B_Voltage = B(1:end,7);

C = readtable('VelocityVoltage_S301_5.csv');
C = table2array(C);
C_Pressure_atm = C(1:end,1);
C_Temperature_atm = C(1:end,2);
C_Airspeed_Diff_Pressure = C(1:end,3);
C_Aux_Diff_Pressure = C(1:end,4);
C_Voltage = C(1:end,7);

D = readtable('VelocityVoltage_S301_7.csv');
D = table2array(D);
D_Pressure_atm = D(1:end,1);
D_Temperature_atm = D(1:end,2);
D_Airspeed_Diff_Pressure = D(1:end,3);
D_Aux_Diff_Pressure = D(1:end,4);
D_Voltage = D(1:end,7);

% Sort data by voltage
% ABCD_Voltage = zeros(10000,1);
% ABCD_Temp_atm = zeros(10000,1);
% ABCD_Pressure_atm = zeros(10000,1);
% ABCD_Aux_Diff_Pressure = zeros(10000,1);
% ABCD_Airspeed_Diff_Pressure = zeros(10000,1);
a = 1;
b = 500;
c = 1;
d = 500;
for i = 1:5
    
ABCD_Voltage(a:b,1) = D_Voltage(c:d,1);
ABCD_Temp_atm(a:b,1) = D_Temperature_atm(c:d,1);
ABCD_Pressure_atm(a:b,1) = D_Pressure_atm(c:d,1);
ABCD_Airspeed_Diff_Pressure(a:b,1) = D_Airspeed_Diff_Pressure(c:d,1);
ABCD_Aux_Diff_Pressure(a:b,1) = D_Aux_Diff_Pressure(c:d,1);
    a = a + 500;
    b = b + 500;
ABCD_Voltage(a:b,1) = A_Voltage(c:d,1);
ABCD_Temp_atm(a:b,1) = A_Temperature_atm(c:d,1);
ABCD_Pressure_atm(a:b,1) = A_Pressure_atm(c:d,1);
ABCD_Airspeed_Diff_Pressure(a:b,1) = A_Airspeed_Diff_Pressure(c:d,1);
ABCD_Aux_Diff_Pressure(a:b,1) = A_Aux_Diff_Pressure(c:d,1);
    a = a + 500;
    b = b + 500;
ABCD_Voltage(a:b,1) = C_Voltage(c:d,1);
ABCD_Temp_atm(a:b,1) = C_Temperature_atm(c:d,1);
ABCD_Pressure_atm(a:b,1) = C_Pressure_atm(c:d,1);
ABCD_Airspeed_Diff_Pressure(a:b,1) = C_Airspeed_Diff_Pressure(c:d,1);
ABCD_Aux_Diff_Pressure(a:b,1) = C_Aux_Diff_Pressure(c:d,1);
    a = a + 500;
    b = b + 500;
ABCD_Voltage(a:b,1) = B_Voltage(c:d,1);
ABCD_Temp_atm(a:b,1) = B_Temperature_atm(c:d,1);
ABCD_Pressure_atm(a:b,1) = B_Pressure_atm(c:d,1);
ABCD_Airspeed_Diff_Pressure(a:b,1) = B_Airspeed_Diff_Pressure(c:d,1);
ABCD_Aux_Diff_Pressure(a:b,1) = B_Aux_Diff_Pressure(c:d,1);
    a = a + 500;
    b = b + 500;
    c = c + 500;
    d = d + 500;
    
end

% calculate velocity for pitot-static + Airspeed Pressure Transducer
R = 287; % gas constant for air
ABCD_Velocity_1 = sqrt(2.*ABCD_Airspeed_Diff_Pressure*R.*ABCD_Temp_atm./ABCD_Pressure_atm);

% plot voltage vs. velocity for pitot-static + Airspeed Pressure Transducer
figure(1)
plot(ABCD_Voltage,ABCD_Velocity_1)
hold on
title('Voltage vs. Velocity')
xlabel('Voltage (V)')
ylabel('Velocity (m/s)')

% calculate velocity for Venturi tube + Water Manometer
Ratio = 1/9.5;
ABCD_Velocity_2 = sqrt(2.*ABCD_Airspeed_Diff_Pressure.*R.*ABCD_Temp_atm./(ABCD_Pressure_atm.*(1-(Ratio)^2)));

% plot voltage vs. velocity for Venturi tube + Water manometer
plot(ABCD_Voltage,ABCD_Velocity_2)
% title('Venturi Tube + Water Manometer')
% xlabel('Voltage (V)')
% ylabel('Velocity (m/s)')


% Part 2
% Read in data

E = readtable('VelocityVoltage_S301_2.csv');
E = table2array(E);
E_Pressure_atm = E(1:end,1);
E_Temperature_atm = E(1:end,2);
E_Airspeed_Diff_Pressure = E(1:end,3);
E_Aux_Diff_Pressure = E(1:end,4);
E_Voltage = E(1:end,7);

F = readtable('VelocityVoltage_S301_4.csv');
F = table2array(F);
F_Pressure_atm = F(1:end,1);
F_Temperature_atm = F(1:end,2);
F_Airspeed_Diff_Pressure = F(1:end,3);
F_Aux_Diff_Pressure = F(1:end,4);
F_Voltage = F(1:end,7);

G = readtable('VelocityVoltage_S301_6.csv');
G = table2array(G);
G_Pressure_atm = G(1:end,1);
G_Temperature_atm = G(1:end,2);
G_Airspeed_Diff_Pressure = G(1:end,3);
G_Aux_Diff_Pressure = G(1:end,4);
G_Voltage = G(1:end,7);

H = readtable('VelocityVoltage_S301_8.csv');
H = table2array(H);
H_Pressure_atm = H(1:end,1);
H_Temperature_atm = H(1:end,2);
H_Airspeed_Diff_Pressure = H(1:end,3);
H_Aux_Diff_Pressure = H(1:end,4);
H_Voltage = H(1:end,7);

% Sort data by voltage
% EFGH_Voltage = zeros(10000,1);
% EFGH_Temp_atm = zeros(10000,1);
% EFGH_Pressure_atm = zeros(10000,1);
% EFGH_Aux_Diff_Pressure = zeros(10000,1);
% EFGH_Airspeed_Diff_Pressure = zeros(10000,1);
a = 1;
b = 500;
c = 1;
d = 500;
for j = 1:5
    
EFGH_Voltage(a:b,1) = H_Voltage(c:d,1);
EFGH_Temp_atm(a:b,1) = H_Temperature_atm(c:d,1);
EFGH_Pressure_atm(a:b,1) = H_Pressure_atm(c:d,1);
EFGH_Airspeed_Diff_Pressure(a:b,1) = H_Airspeed_Diff_Pressure(c:d,1);
EFGH_Aux_Diff_Pressure(a:b,1) = H_Aux_Diff_Pressure(c:d,1);
    a = a + 500;
    b = b + 500;
EFGH_Voltage(a:b,1) = E_Voltage(c:d,1);
EFGH_Temp_atm(a:b,1) = E_Temperature_atm(c:d,1);
EFGH_Pressure_atm(a:b,1) = E_Pressure_atm(c:d,1);
EFGH_Airspeed_Diff_Pressure(a:b,1) = E_Airspeed_Diff_Pressure(c:d,1);
EFGH_Aux_Diff_Pressure(a:b,1) = E_Aux_Diff_Pressure(c:d,1);
    a = a + 500;
    b = b + 500;
EFGH_Voltage(a:b,1) = G_Voltage(c:d,1);
EFGH_Temp_atm(a:b,1) = G_Temperature_atm(c:d,1);
EFGH_Pressure_atm(a:b,1) = G_Pressure_atm(c:d,1);
EFGH_Airspeed_Diff_Pressure(a:b,1) = G_Airspeed_Diff_Pressure(c:d,1);
EFGH_Aux_Diff_Pressure(a:b,1) = G_Aux_Diff_Pressure(c:d,1);
    a = a + 500;
    b = b + 500;
EFGH_Voltage(a:b,1) = F_Voltage(c:d,1);
EFGH_Temp_atm(a:b,1) = F_Temperature_atm(c:d,1);
EFGH_Pressure_atm(a:b,1) = F_Pressure_atm(c:d,1);
EFGH_Airspeed_Diff_Pressure(a:b,1) = F_Airspeed_Diff_Pressure(c:d,1);
EFGH_Aux_Diff_Pressure(a:b,1) = F_Aux_Diff_Pressure(c:d,1);
    a = a + 500;
    b = b + 500;
    c = c + 500;
    d = d + 500;
    
end

% calculate velocity for pitot-static + Water Manometer
EFGH_Velocity_1 = sqrt(2.*EFGH_Airspeed_Diff_Pressure*R.*EFGH_Temp_atm./EFGH_Pressure_atm);

% plot voltage vs. velocity for pitot-static + Water Manometer
% figure(3)
plot(EFGH_Voltage,EFGH_Velocity_1)
% title('Pitot-Static + Water Manometer')
% xlabel('Voltage (V)')
% ylabel('Velocity (m/s)')

% calculate velocity for Venturi tube + Airspeed Pressure Transducer
Ratio = 1/9.5;
EFGH_Velocity_2 = sqrt(2.*EFGH_Airspeed_Diff_Pressure.*R.*EFGH_Temp_atm./(EFGH_Pressure_atm.*(1-(Ratio)^2)));

% plot voltage vs. velocity for Venturi tube + Airspeed Pressure Transducer
% figure(4)
plot(EFGH_Voltage,EFGH_Velocity_2)
% title('Venturi Tube + Airspeed Pressure Transducer')
% xlabel('Voltage (V)')
% ylabel('Velocity (m/s)')
legend('Pitot-Static + Airspeed Pressure Transducer', 'Venturi Tube + Water Manometer', 'Pitot-Static + Water Manometer', 'Venturi Tube + Airspeed Pressure Transducer')
hold off