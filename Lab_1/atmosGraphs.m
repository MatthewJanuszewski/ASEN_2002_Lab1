clc
clear
close all

altitude = 0:500:30000;
pressure = zeros(1, 60);
temperature = zeros(1, 60);
density = zeros(1, 60);

for i = 1:61
    [temperature(i), ~, pressure(i), density(i)] = atmoscoesa((i-1)*500, 'None'); % [k, m/s, Pa, kg/m^3]
end

altitude = altitude/1000;
pressure = pressure/1000;

figure
plot(altitude, pressure,'LineWidth', 2)
title('1976 Standard Atmosphere Pressure vs Altitude')
xlabel('Altitude (km)')
ylabel('Pressure (kPa)')

figure
plot(altitude, density,'LineWidth', 2)
title('1976 Standard Atmosphere Density vs Altitude')
xlabel('Altitude (km)')
ylabel('Density (kg/m^3)')

figure
plot(altitude, temperature,'LineWidth', 2)
title('1976 Standard Atmosphere Temperature vs Altitude')
xlabel('Altitude (km)')
ylabel('Temperature (K)')


% Balloon Material Properties 
% (average of values given on matweb for polyester film)
density_material = 1255; % [kg/m^3]
FS = 1:0.01:3; % Factor of Safety
YS = 27.6*10^6; % [Pa]

t = FS.*((10*15.56)/(2*YS));

figure
plot(FS, t,'LineWidth', 2)
title('Shell Thickness vs Safety Factor Using Calculated Radius of 15.56 m')
xlabel('Factor of Safety')
ylabel('Shell Thickness (m)')

