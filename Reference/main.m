%% ASEN 2012 Project 2
% Main function for Project 2

%   Purpose: Determine the bottle rocket thrust as a function of time, and
%   predict the resulting height and range of the rocket. You are then 
%   asked to use the code to explore the parameter space in order to 
%   determine how each of the parameters affect the height and the range 
%   of the rocket, and what combination of parameters will allow the rocket
%   to land within 1 meter of a 80 meter marker.

%   Imputs: (possibly) initial gage pressure, initial volume of water in 
%   bottle, coefficient of drag, and launch angle

% Author(s): Noah Freeland (108201253)
%Zach Mason (107479299)
% Date Created: 12/3/18
% Date Modified: 12/10/19

%% Housekeeping

clear all
close all
clc


%% Constants, Parameters, and Variables

% Constants
g = 9.81; % acceleration of gravity [m/s^2]
C_d = 0.8; % discharge coefficient
rho_air = 0.961; % density of air [kg/m^3]
rho_water = 1000; % density of water [kg/m^3]
gamma = 1.4; % specific heat ratio for air
R = 287; % gas constant of air [J/kgK]
% convec = [g, C_d, rho_air, rho_water, gamma, R];

% Parameters
V_b = 0.002; % volume of the bottle [m^3]
M_b = 0.15; % mass of empty bottle rocket [kg]
P_a = 12.1; % atmospheric (ambient) pressure [psi]
P_a = P_a*6894.76; % [Pa] convert to SI
D_t = 2.1; % diameter of the throat [cm]
D_t = D_t/100; % [m] convert to SI
D_b = 10.5; % diamerer of the bottle [cm]
D_b = D_b/100; % [m] convert to SI
T_air_i = 300; % initial temperature of air [K]
v0 = 0.0; % initial velocity of rocket [m/s]
x0 = 0.0; % initial horizontal distance [m]
z0 = 0.25; % initial vertical distance [m]
l_s = 0.5; % length of test stand [m]

t_span = [0 10]; % duration of flight (integration time) [sec]

% cross sectional areas calculated
A_t = ((D_t/2)^2)*pi; % throat [m^2]
A_b = ((D_b/2)^2)*pi; % bottle [m^2]
% Avec = [A_t, A_bottle];

% thrust 
global Thrust Time
Thrust = 0;
Time = 0;

    % Variables (verification case values currently)
    V_water_i = 0.001; % initial volume of water inside bottle [m^3]
    C_D = 0.3;
    theta = 42; % initial angle of launch [degrees]
    theta = (2*pi)*theta/360; % convert to SI [radians]
    % varvec = [P_g, V_water_i, C_D, theta];
    
    P_g = 54; % initial gage pressure of air in bottle [psi]
    P_g = P_g*6894.76; % [Pa] convert to SI
    
    % additional initial conditions/parameters
    V_air_i = V_b - V_water_i; % initial volume of air [m^3]
    P_air_i = P_a + P_g; % total initial pressure of air inside bottle
    M_tot_i = M_b + (rho_water * (V_b - V_air_i)) + (P_air_i/(R*T_air_i))*V_air_i;
    
    ic = [V_air_i, M_tot_i, P_air_i, x0, z0, 0, 0]; % initial conditions

    options = odeset('RelTol',1e-8);
    [t,y] = ode45(@(t,y) getodes(t,y,P_g,V_water_i,C_D,theta),t_span,ic,options);

    Vol = y(:,1); % volume function V(t)
    M = y(:,2); % mass function M(t)
    P = y(:,3); % pressure function P(t)
    x = y(:,4); % horizontal position function x(t)
    z = y(:,5); % vertical position function y(t)
    dxdt = y(:,6); % horizontal velocity component over time
    dzdt = y(:,7); % vertical velocity component over time

% P_g_vec = 40:1:80;
% figure(1)
% % plot thrust
% plot(P_g_vec,distMat);
% hold on
% yline(80,'r');
% grid on
% xlabel('P_g [Pa]');
% ylabel('max horizontal distance [m]');
% title('Distance vs P_g for V^i_w = 0.001 m^3, C_D = 0.5, \theta^i = 45^o');

[Thrust,TF] = rmoutliers(Thrust,'movmean',[10 10],'ThresholdFactor',1);
Thrust = [0;Thrust];

Time = sort(Time);

r = find(TF);
Time(r) = [];
Time = [0;Time];

figure(1)
% plot thrust
plot(sort(Time),Thrust,'LineWidth',2);
grid on
xlabel('Time (sec)');
ylabel('Thrust (N)');
xlim([0 0.45]);
ylim([0 300]);
title('Thrust vs Time');

% trajectory
figure(2)
plot(x,z,'LineWidth',2);
grid on;
title('Height vs Distance');
xlabel('Distance (m)');
ylabel('Height (m)');
xlim([0 90]);
ylim([0 40]);

% air volume vs time
figure(3)
plot(t,Vol,'LineWidth',2);
grid on;
title('Air Volume vs Time');
xlabel('Distance (m)');
ylabel('Height (m)');
xlim([0 0.25]);
ylim([0 0.004]);

% horizontal velocity
figure(4)
plot(t,dxdt,'LineWidth',2);
grid on;
title('Horizontal Velocity vs Time');
xlabel('Time (sec)');
ylabel('Horizontal [m/s]');

% vertical velocity
figure(5)
plot(t,dzdt,'LineWidth',2);
grid on;
title('Vertical Velocity vs Time');
xlabel('Time (sec)');
ylabel('Vertical [m/s]');

% Air Pressure vs time
figure(6)
plot(t,P,'LineWidth',2);
grid on;
title('Internal Air Pressure vs Time');
xlabel('Time (sec)');
ylabel('Air Pressure [Pa]');
xlim([0 0.5]);

% Rocket Mass vs time
figure(7)
plot(t,M,'LineWidth',2);
grid on;
title('Rocket Mass vs Time');
xlabel('Time (sec)');
ylabel('Rocket Mass [kg]');
xlim([0 0.5]);


%find max distance and max height
value = find(z<0);
maxDist = x(value(1))
maxHeight = max(z)