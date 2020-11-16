% Max Martinez
% ASEN 2002 Lab 2
% Group 5
% ----------------------------------------------------------------------- %
clc
clear

R = 287;  % [J/kg-K]
A2 = 144*0.0254; % [m]
A1 = A2*9.5; % [m]
n = 1;
m = 1;
r = 1;
q = 1;

Pitot{1,4} = [];
P = zeros(20,7);
VelVol{1,4} = [];
VV = zeros(20,7);

% ----------------------------------------------------------------------- %
% /// Part 5: Airspeed Calculations /////////////////////////////////////
% ----------------------------------------------------------------------- %
% Pitot-Static Probe

for i = [1 3 5 7]% Reading the files
    Pitot{n} = readmatrix(sprintf(...
        'PitotProbeToPressureTransducer/VelocityVoltage_S301_%d.csv',i));
    n = n+1;
end

for I = linspace(1,5,5)% Concatenating and rearanging based on Voltage
    for J = [4 1 3 2]
        P(q,:) = mean(Pitot{1,J}((I*500)-499:I*500,1:7));
        q = q+1;
    end
end

P_atm = P(:,1); % [Pa] - Atmospheric Pressure
T_atm = P(:,2); % [T]  - Atmospheric Temperature
diffP = P(:,3); % [Pa] - Air Pressure Differential
PauxP = P(:,4); % [Pa] - Aux Pressure Differential
eldxP = P(:,5); % [mm] - ELD Probe X axis
eldyP = P(:,6); % [mm] - ELD Probe Y axis
voltP = P(:,7); % [V]  - Voltage

airSpdP = sqrt((2.*R.*T_atm.*diffP)./P_atm); % [m/s]

% ----------------------------------------------------------------------- %
% Venturi-Tube

for k = [2 4 6 8] % Reading the files
    VelVol{m} = readmatrix(sprintf...
        ('VenturiTubeToPressureTransducer/VelocityVoltage_S301_%d.csv',k));
    m = m+1;
end

for K = linspace(1,5,5)% Concatenating and rearanging based on Voltage
    for L = [4 1 3 2]
        VV(r,:) = mean(VelVol{1,L}((K*500)-499:K*500,1:7));
        r = r+1;
    end
end

PatmV = VV(:,1); % [Pa] - Atmospheric Pressure
TatmV = VV(:,2); % [T]  - Atmospheric Temperature 
diffPV= VV(:,3); % [Pa] - Air Pressure Differential
PauxV = VV(:,4); % [Pa] - Aux Pressure Differential
eldxV = VV(:,5); % [mm] - ELD Probe X axis
eldyV = VV(:,6); % [mm] - ELD Probe Y axis
voltV = VV(:,7); % [V]  - Voltage

airSpdV = sqrt((2*R.*TatmV.*diffPV)./(PatmV.*((1-(A2/A1)^2)))); % [m/s]

% ----------------------------------------------------------------------- %
% Uncertainty

% Pitot-Static Tube
dVdP = mean(-(R.*diffP.*T_atm)./...
    (((2^0.5).*((R.*diffP.*T_atm)./P_atm)).*P_atm.^2));
dVdT = mean((diffP.*R)./(((2^0.5).*P_atm).*(diffP.*R.*T_atm)));
dVdp = mean((R.*T_atm)./(((2^0.5).*P_atm).*(diffP.*R.*T_atm)));

deltaP = 68.9476;  % std(P_atm); - [Pa]
deltaT = 0.25;     % std(T_atm); - [K]
deltap = 3450;     % std(diffP); - [Pa]

% uVP = sqrt((dVdP*deltaP)^2 + (dVdT*deltaT)^2 + (dVdp*deltap)^2);
% errP = uVP.*ones(size(airSpdP));

err = std(airSpdP);
Err = err.*ones(size(airSpdP));

figure 
errorbar(voltP,airSpdP,Err,'LineWidth',2)
xlabel('Input Voltage (V)');
ylabel('Resulting Airspeed (m/s)');
title('Airspeed for Pitot-Static Probe');



% Venturi Tube

% dVdTV = mean(sqrt((2.*diffPV.*R)./(TatmV.*PatmV.*(1-(A2/A1)^2)))./2);
% dVdPV = mean(-(sqrt((2.*R.*TatmV.*diffPV)./(1-(A2/A1)^2)).*(PatmV.^(3/2))));
% dVdpV = mean(sqrt((2.*R.*TatmV)./(diffPV.*PatmV))./2);

% dVdTV = mean((diffPV.*R)./(sqrt(2).*(1-(A2/A1)^2).*PatmV.*sqrt(...
%     (diffPV.*R.*TatmV)./((1-(A2/A1)^2).*PatmV))));
% dVdpV = mean((R.*TatmV)./(sqrt(2).*(1-(A2/A1)^2).*PatmV.*sqrt(...
%     (R.*TatmV.*diffPV)./(PatmV.*(1-(A2/A1)^2)))));
% dVdPV = mean(((-1).*R.*TatmV.*diffPV)./(sqrt(2).*(1-(A2/A1)^2).*sqrt(...
%     (R.*TatmV.*diffPV)./((1-(A2/A1)^2).*PatmV)).*(PatmV).^2));

% dVdTV = mean(gradient(airSpdV,TatmV));
% dVdPV = mean(gradient(airSpdV,PatmV));
% dVdpV = mean(gradient(airSpdV,diffPV));

% dVdTV = mean(sqrt((2.*diffPV.*R)./(TatmV.*PatmV.*(1-(A2/A1)^2)))./2);
% dVdPV = mean(sqrt((2.*R.*TatmV.*diffPV)./(1-(A2/A1)^2)).*(-0.5).*(PatmV.^(3/2)));
% dVdpV = mean(sqrt((2.*R.*TatmV)./(diffPV.*PatmV.*(1-(A2/A1)^2)))./2);

% dPV = 68.9476;  % std(PatmV);
% dTV = 0.25;     % std(TatmV);
% dpV = 3450;     % std(diffPV);
% 
% uVV = sqrt((dVdPV*dPV)^2 + (dVdTV*dTV)^2 + (dVdpV*dpV)^2);
% errV = uVV.*ones(size(airSpdV));

e = std(airSpdV);
E = e.*ones(size(airSpdV));

figure 
errorbar(voltV,airSpdV,E,'LineWidth',2)
xlabel('Input Voltage (V)')
ylabel('Resulting Airspeed (m/s)')
title('Airspeed Reading for Venturi Tube')


% ----------------------------------------------------------------------- %
% /// Part 6: Boundary Layer ////////////////////////////////////////////
% ----------------------------------------------------------------------- %

% Distances [mm]
P1d = 9.05*25.4; P2d = 10.03*25.4; P3d = 11.01*25.4; P4d = 11.99*25.4;
P5d = 12.97*25.4; P6d = 13.95*25.4; P7d = 14.93*25.4; P8d = 15.91*25.4;
P9d = 16.89*25.4; P10d = 17.87*25.4; P11d = 18.85*25.4;

% Preallocating cells and matrices
Port1{1,1} = []; Port2{1,1} = []; Port3{1,1} = []; Port4{1,1} = [];
Port5{1,1} = []; Port6{1,1} = []; Port7{1,1} = []; Port8{1,1} = [];
Port9{1,1} = []; Port10{1,1} = []; Port11{1,1} = []; 

P1 = zeros(12,7); P2 = zeros(12,7); P3 = zeros(12,7); P4 = zeros(12,7);
P5 = zeros(12,7); P6 = zeros(12,7); P7 = zeros(12,7); P8 = zeros(12,7);
P9 = zeros(12,7); P10 = zeros(12,7); P11 = zeros(12,7); 

for l = 1 % Loading in the data
    Port1{l} = readmatrix(sprintf(...
        'BoundaryLayerData/Port1/BoundaryLayer_S301_1.csv'));
    Port2{l} = readmatrix(sprintf(...
        'BoundaryLayerData/Port2/BoundaryLayer_S302_1.csv'));
    Port3{l} = readmatrix(sprintf(...
        'BoundaryLayerData/Port3/BoundaryLayer_S303_1.csv'));
    Port4{l} = readmatrix(sprintf(...
        'BoundaryLayerData/Port4/BoundaryLayer_S301_3.csv'));
    Port5{l} = readmatrix(sprintf(...
        'BoundaryLayerData/Port5/BoundaryLayer_S302_3.csv'));
    Port6{l} = readmatrix(sprintf(...
        'BoundaryLayerData/Port6/BoundaryLayer_S301_6.csv'));
    Port7{l} = readmatrix(sprintf(...
        'BoundaryLayerData/Port7/BoundaryLayer_S301_5.csv'));
    Port8{l} = readmatrix(sprintf(...
        'BoundaryLayerData/Port8/BoundaryLayer_S302_5.csv'));
    Port9{l} = readmatrix(sprintf(...
        'BoundaryLayerData/Port9/BoundaryLayer_S303_7.csv'));
    Port10{l} = readmatrix(sprintf(...
        'BoundaryLayerData/Port10/BoundaryLayer_S301_7.csv'));
    Port11{l} = readmatrix(sprintf(...
        'BoundaryLayerData/Port11/BoundaryLayer_S302_7.csv'));
end

y = 1;
for x = linspace(1,12,12) % Concatenating the matrices
    P1(y,:) = mean(Port1{1,1}((x*500)-499:x*500,1:7));
    P2(y,:) = mean(Port2{1,1}((x*500)-499:x*500,1:7));
    P3(y,:) = mean(Port3{1,1}((x*500)-499:x*500,1:7));
    P4(y,:) = mean(Port4{1,1}((x*500)-499:x*500,1:7));
    P5(y,:) = mean(Port5{1,1}((x*500)-499:x*500,1:7));
    P6(y,:) = mean(Port6{1,1}((x*500)-499:x*500,1:7));
    P7(y,:) = mean(Port7{1,1}((x*500)-499:x*500,1:7));
    P8(y,:) = mean(Port8{1,1}((x*500)-499:x*500,1:7));
    P9(y,:) = mean(Port9{1,1}((x*500)-499:x*500,1:7));
    P10(y,:) = mean(Port10{1,1}((x*500)-499:x*500,1:7));
    P11(y,:) = mean(Port11{1,1}((x*500)-499:x*500,1:7));
    y = y+1;
end

P1T = P1(:,2); P1dP = P1(:,4); P1Patm = P1(:,1); P1y = P1(:,6);
velP1 = sqrt((2.*R.*P1T.*P1dP)./P1Patm); % Pitot-Static

P2T = P2(:,2); P2dP = P2(:,4); P2Patm = P2(:,1); P2y = P2(:,6);
velP2 = sqrt((2.*R.*P2T.*P2dP)./P2Patm); % Pitot-Static

P3T = P3(:,2); P3dP = P3(:,4); P3Patm = P3(:,1); P3y = P3(:,6);
velP3 = sqrt((2.*R.*P3T.*P3dP)./P3Patm); % Pitot-Static

P4T = P4(:,2); P4dP = P4(:,4); P4Patm = P4(:,1); P4y = P4(:,6);
velP4 = sqrt((2.*R.*P4T.*P4dP)./P4Patm); % Pitot-Static

P5T = P5(:,2); P5dP = P5(:,4); P5Patm = P5(:,1); P5y = P5(:,6);
velP5 = sqrt((2.*R.*P5T.*P5dP)./P5Patm); % Pitot-Static

P6T = P6(:,2); P6dP = P6(:,4); P6Patm = P6(:,1); P6y = P6(:,6);
velP6 = sqrt((2.*R.*P6T.*P6dP)./P6Patm); % Pitot-Static

P7T = P7(:,2); P7dP = P7(:,4); P7Patm = P7(:,1); P7y = P7(:,6);
velP7 = sqrt((2.*R.*P7T.*P7dP)./P7Patm); % Pitot-Static

P8T = P8(:,2); P8dP = P8(:,4); P8Patm = P8(:,1); P8y = P8(:,6);
velP8 = sqrt((2.*R.*P8T.*P8dP)./P8Patm); % Pitot-Static

P9T = P9(:,2); P9dP = P9(:,4); P9Patm = P9(:,1); P9y = P9(:,6);
velP9 = sqrt((2.*R.*P9T.*P9dP)./P9Patm); % Pitot-Static

P10T = P10(:,2); P10dP = P10(:,4); P10Patm = P10(:,1); P10y = P10(:,6);
velP10 = sqrt((2.*R.*P10T.*P10dP)./P10Patm); % Pitot-Static

P11T = P11(:,2); P11dP = P11(:,4); P11Patm = P11(:,1); P11y = P11(:,6);
velP11 = sqrt((2.*R.*P11T.*P11dP)./P11Patm); % Pitot-Static

figure
hold on
plot(P1y,velP1)
plot(P2y,velP2)
plot(P3y,velP3)
plot(P4y,velP4)
plot(P5y,velP5)
plot(P6y,velP6)
plot(P7y,velP7)
plot(P8y,velP8)
plot(P9y,velP9)
plot(P10y,velP10)
plot(P11y,velP11)
xlabel('ELD Y Distance [mm]');
ylabel('Velocity (m/s)');
title('Airspeed Relative to ELD Probe Y Distance');
hold off

V = [velP1,velP2,velP3,velP4,velP5,velP6,velP7,velP8,velP9,velP10,velP11];
Vfreestream = mean(V(end));
thickness = (0.95)*Vfreestream;


y = thickness./V(end,:);
ports = linspace(1,11,11);
figure
plot(ports,y)
xlabel('Ports')
ylabel('Boundary Layer Thickness [mm]')
title('Boundary Layer Thickness For Each Port')

% Determine if boundary layer is laminar or turbulent!!!!!


% ----------------------------------------------------------------------- %
% /// Part 7: AirFoil ///////////////////////////////////////////////////
% ----------------------------------------------------------------------- %

