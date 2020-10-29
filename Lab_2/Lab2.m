% Max Martinez
% ASEN 2002 Lab 2
% Group 5
% ----------------------------------------------------------------------- %
clc
clear

R = 287;  % [J/kg-K]
n = 1;
A1 = 9.5*42.5; % [in?]
A2 = 9.5*24;
VelVol301{1,4} = [];
P_atm{1,4} = [];
T_atm{1,4} = [];
P_air{1,4} = [];
P_aux{1,4} = [];
eldpx{1,4} = [];
eldpy{1,4} = [];
volt{1,4}  = [];
airSpeed{1,4}  = [];
airSpdAlt{1,4} = [];

for i = linspace(2,8,4)
    VelVol301{n} = readmatrix(sprintf...
        ('VenturiTubeToPressureTransducer/VelocityVoltage_S301_%d.csv',i));
    P_atm{n} = VelVol301{1,n}(:,1);   % [Pa]- Atmospheric Pressure
    T_atm{n} = VelVol301{1,n}(:,2);   % [K] - Atmospheric Temperature
    P_air{n} = VelVol301{1,n}(:,3);   % [Pa]- Air Pressure Differential
    P_aux{n} = VelVol301{1,n}(:,4);   % [Pa]- Aux Pressure Differential
    eldpx{n} = VelVol301{1,n}(:,5);   % [mm]- ELD Probe X axis
    eldpy{n} = VelVol301{1,n}(:,6);   % [mm]- ELD Probe Y axis
    volt{n}  = VelVol301{1,n}(:,7);   % [V] - Voltage
    
    % 5.2 Airspeed Formula Comparisons
    airSpeed{n} = sqrt((2.*R.*T_atm{1,n}.*abs((P_air{1,n}-P_aux{1,n})))...
        ./P_atm{1,n});
    airSpdAlt{n} = sqrt((2.*R.*T_atm{1,n}.*abs((P_air{1,n}-P_aux{1,n}))...
        ./P_atm{1,n}.*((1-A2/A1)^2)));
    figure(1)
        subplot(2,2,n);
        plot(linspace(1,2500,2500),airSpeed{:,n},'r')
        xlabel('Time?')
        ylabel('Airspeed (m/s)')
        title(i)
    figure(2)
        subplot(2,2,n);
        plot(linspace(1,2500,2500),airSpdAlt{:,n},'b')
        xlabel('Time?')
        ylabel('Airspeed (m/s)')
        title(i)

    
    % 5.3 Airspeed Model with respect to Voltage
    figure(3)
        subplot(2,2,n)
        plot(volt{n},airSpeed{:,n},'k')
        xlabel('Voltage (V)')
        ylabel('Airspeed (m/s)')
        title(i)
    
    n = n+1;
end

% ----------------------------------------------------------------------- %