% Max Martinez
% ASEN 2002 Lab 2
% Group 5
% ----------------------------------------------------------------------- %
clc
clear
R = 287;  % [J/kg-K]
A2 = 144*0.0254; % [m]
A1 = A2*9.5; % [m]
n = 1; m = 1; r = 1; q = 1;
Pitot{1,4} = []; P = zeros(20,7); VelVol{1,4} = []; VV = zeros(20,7);

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
% Uncertainty for Pressure Transducer

% Pitot-Static Tube
dVdT = (0.5).*sqrt((2*R*diffP)./(P_atm.*T_atm));
dVdP = sqrt(2.*R.*T_atm.*diffP).*(-0.5).*(P_atm.^(-3/2));
dVdp = (0.5).*sqrt((2.*R.*T_atm)./(P_atm.*diffP));

deltaP = 3450;      % - [Pa]
deltaT = 0.25;      % - [K]
deltap = 68.9476;   % - [Pa]
uVP = sqrt((dVdP.*deltaP).^2 + (dVdT.*deltaT).^2 + (dVdp.*deltap).^2);
errP = uVP.*ones(size(airSpdP));
EP = errP.*ones(size(airSpdP));

% Venturi Tube
dVdTV = (sqrt((2.*diffPV.*R)./(TatmV.*PatmV.*(1-(A2/A1)^2)))./2);
dVdPV = (-(sqrt((2.*R.*TatmV.*diffPV)./(1-(A2/A1)^2)).*(PatmV.^(-3/2))));
dVdpV = (sqrt((2.*R.*TatmV)./(diffPV.*PatmV.*(1-(A2/A1)^2)))./2);
dPV = 3450; 
dTV = 0.25;     
dpV = 68.9476;     
uVV = sqrt((dVdPV.*dPV).^2 + (dVdTV.*dTV).^2 + (dVdpV.*dpV).^2);
errV = uVV.*ones(size(airSpdV));
EV = errV.*ones(size(airSpdV));

figure 
hold on
errorbar(voltV,airSpdV,EV,'LineWidth',2)
errorbar(voltP,airSpdP,EP,'LineWidth',2)
xlabel('Input Voltage (V)')
ylabel('Velocity (m/s)')
title('Airspeed for Electronic Pressure Transducer Readings')
legend('Venturi Tube','Pitot-Static Probe')
hold off

% ----------------------------------------------------------------------- %
% Water Monometer and Uncertainty
WM = readmatrix(sprintf(...
    'Water manometer readings (Responses).xlsx'));

T_atmWM = T_atm;
P_atmWM = P_atm;
diffPwm = 248.81.*[WM(26,5) WM(12,5) WM(19,5) WM(18,5)...
    WM(26,7) WM(12,7) WM(19,7) WM(18,7)...
    WM(26,9) WM(12,9) WM(19,9) WM(18,9)...
    WM(26,11) WM(12,11) WM(19,11) WM(18,11)...
    WM(26,13) WM(12,13) WM(19,13) WM(18,13)];
voltPwm = [WM(26,4) WM(12,4) WM(19,4) WM(18,4)...
    WM(26,6) WM(12,6) WM(19,6) WM(18,6)...
    WM(26,8) WM(12,8) WM(19,8) WM(18,8)...
    WM(26,10) WM(12,10) WM(19,10) WM(18,10)...
    WM(26,12) WM(12,12) WM(19,12) WM(18,12)];
airSpdPWM = sqrt((2.*R.*T_atmWM.*diffPwm)./P_atmWM);
airSpdPwm = airSpdPWM(1,:);
dVdTPwm = (sqrt((2.*diffPwm.*R)./(T_atmWM.*P_atmWM.*(1-(A2/A1)^2)))./2);
dVdPPwm = (-(sqrt((2.*R.*T_atmWM.*diffPwm)./(1-(A2/A1)^2)).*(P_atmWM.^(-3/2))));
dVdpPwm = (sqrt((2.*R.*T_atmWM)./(diffPwm.*P_atmWM.*(1-(A2/A1)^2)))./2);
dPVwm = 3450; 
dTVwm = 0.25;     
dpVwm = 68.9476;     
uVVwm = sqrt((dVdPPwm.*dPVwm).^2+(dVdTPwm.*dTVwm).^2+(dVdpPwm.*dpVwm).^2);
errVwm = uVVwm.*ones(size(airSpdPwm));
EVwm = errVwm.*ones(size(airSpdPwm));


TatmWM = TatmV;
PatmWM = PatmV;
diffVWM = 248.81.*[WM(22,5) WM(6,5) WM(10,5) WM(13,5) ...
    WM(22,7) WM(6,7) WM(10,7) WM(13,7)...
    WM(22,9) WM(6,9) WM(10,9) WM(13,9)...
    WM(22,11) WM(6,11) WM(10,11) WM(13,11) ...
    WM(22,13) WM(6,13) WM(10,13) WM(13,13)];
voltVWM = [WM(22,4) WM(6,4) WM(10,4) WM(13,4) ...
    WM(22,6) WM(6,6) WM(10,6) WM(13,6)...
    WM(22,8) WM(6,8) WM(10,8) WM(13,8)...
    WM(22,10) WM(6,10) WM(10,10) WM(13,10)...
    WM(22,12) WM(6,12) WM(10,12) WM(13,12)];
airSpdVwm = sqrt((2.*R.*TatmWM.*diffVWM)./PatmWM);
airSpdVWM = airSpdVwm(1,:);
dVdTVWM = (sqrt((2.*diffVWM.*R)./(TatmWM.*PatmWM.*(1-(A2/A1)^2)))./2);
dVdPVWM = (-(sqrt((2.*R.*TatmWM.*diffVWM)./(1-(A2/A1)^2)).*(PatmWM.^(-3/2))));
dVdpVWM = (sqrt((2.*R.*TatmWM)./(diffVWM.*PatmWM.*(1-(A2/A1)^2)))./2);
dPVWM = 3450; 
dTVWM = 0.25;     
dpVWM = 68.9476;     
uVVWM = sqrt((dVdPVWM.*dPVWM).^2+(dVdTVWM.*dTVWM).^2+(dVdpVWM.*dpVWM).^2);
errVWM = uVVWM.*ones(size(airSpdVWM));
EVWM = errVWM.*ones(size(airSpdVWM));

figure
hold on
errorbar(voltVWM,airSpdVWM,EVWM(1,:),'LineWidth',2)
errorbar(voltPwm,airSpdPwm(1,:),EVwm(1,:),'LineWidth',2)
xlabel('Input Voltage (V)')
ylabel('Velocity (m/s)')
title('Airspeed from Water Manometer Readings')
legend('Venturi Tube','Pitot-Static Probe')
hold on


figure
hold on
plot(voltP,airSpdP,'LineWidth',2)
plot(voltV,airSpdV,'LineWidth',2)
plot(voltPwm,airSpdPwm(1,:),'LineWidth',2)
plot(voltVWM,airSpdVWM(1,:),'LineWidth',2)
ylabel('Velocity (m/s)')
xlabel('Voltage (V)')
title('Airspeed Calulations Using Electronic Pressure Transducer and Water Manometer')
legend('Pitot-Static Pressure Transducer','Venturi Pressure Transducer',...
    'Pitot-Static Water Manometer','Venturi Water Manometer','Location','NorthWest')
hold off

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
% % % plot(velP1,P1y)
% % % plot(velP2,P2y)
% % % plot(velP3,P3y)
% % % plot(velP4,P4y)
% % % plot(velP5,P5y)
% % % plot(velP6,P6y)
% % % plot(velP7,P7y)
% % % plot(velP8,P8y)
% % % plot(velP9,P9y)
% % % plot(velP10,P10y)
% % % plot(velP11,P11y)

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

rho1 = mean(P1Patm./(R.*P1T)); 
rho2 = mean(P2Patm./(R.*P2T));
rho3 = mean(P3Patm./(R.*P3T)); 
rho4 = mean(P4Patm./(R.*P4T));
rho5 = mean(P5Patm./(R.*P5T)); 
rho6 = mean(P6Patm./(R.*P6T));
rho7 = mean(P7Patm./(R.*P7T)); 
rho8 = mean(P8Patm./(R.*P8T));
rho9 = mean(P9Patm./(R.*P9T)); 
rho10 = mean(P10Patm./(R.*P10T));
rho11 = mean(P11Patm./(R.*P11T));

rho =(rho1 + rho2 + rho3 + rho4 + rho5 +rho6 + rho7 + rho8 + rho9 + ...
    rho10 + rho11)/11;

V = [velP1,velP2,velP3,velP4,velP5,velP6,velP7,velP8,velP9,velP10,velP11];
Vfreestream = mean(V(end)); 
vbound = (0.95)*Vfreestream;
ports = linspace(1,11,11);

PortVelocityProbe = [sqrt(2.*P1dP/rho),sqrt(2.*P2dP/rho),...
    sqrt(2.*P3dP/rho),sqrt(2.*P4dP/rho),sqrt(2.*P5dP/rho),...
    sqrt(2.*P6dP/rho),sqrt(2.*P7dP/rho),sqrt(2.*P8dP/rho),...
    sqrt(2.*P9dP/rho),sqrt(2.*P10dP/rho),sqrt(2.*P11dP/rho)];
Pyaxis = [P1y P2y P3y P4y P5y P6y P7y P8y P9y P10y P11y];

a = (PortVelocityProbe > vbound);
blh = [Pyaxis(6,1);Pyaxis(6,2);Pyaxis(6,3);Pyaxis(6,4);Pyaxis(8,5);...
    Pyaxis(9,6);Pyaxis(9,7);Pyaxis(10,8);Pyaxis(9,9);Pyaxis(9,10);...
    Pyaxis(10,11)];

% Determine if boundary layer is laminar or turbulent!!!!!

mu = 1.7894*(10^-5);

Rex1 =((P1d*velP1(end)*rho)/mu)/1000; 
Rex2 =((P2d*velP2(end)*rho)/mu)/1000;
Rex3 =((P3d*velP3(end)*rho)/mu)/1000; 
Rex4 =((P4d*velP4(end)*rho)/mu)/1000;
Rex5 =((P5d*velP5(end)*rho)/mu)/1000; 
Rex6 =((P6d*velP6(end)*rho)/mu)/1000;
Rex7 =((P7d*velP7(end)*rho)/mu)/1000; 
Rex8 =((P8d*velP8(end)*rho)/mu)/1000;
Rex9 =((P9d*velP9(end)*rho)/mu)/1000;
Rex10 =((P10d*velP10(end)*rho)/mu)/1000;
Rex11 =((P11d*velP11(end)*rho)/mu)/1000;

thickTurb = [(0.37*P1d)/(Rex1^(0.2));(0.37*P2d)/(Rex2^(0.2));...
    (0.37*P3d)/(Rex3^(0.2));(0.37*P4d)/(Rex4^(0.2));...
    (0.37*P5d)/(Rex5^(0.2));(0.37*P6d)/(Rex6^(0.2));...
    (0.37*P7d)/(Rex7^(0.2));(0.37*P8d)/(Rex8^(0.2));...
    (0.37*P9d)/(Rex9^(0.2));(0.37*P10d)/(Rex10^(0.2));...
    (0.37*P11d)/(Rex11^(0.2))];

thickLam = [(5.2*P1d)/sqrt(Rex1);(5.2*P2d)/sqrt(Rex2);...
    (5.2*P3d)/sqrt(Rex3);(5.2*P4d)/sqrt(Rex4);(5.2*P5d)/sqrt(Rex5);...
    (5.2*P6d)/sqrt(Rex6);(5.2*P7d)/sqrt(Rex7);(5.2*P8d)/sqrt(Rex8);...
    (5.2*P9d)/sqrt(Rex9);(5.2*P10d)/sqrt(Rex10);(5.2*P11d)/sqrt(Rex11)];

figure
hold on
plot(ports,thickTurb,'LineWidth',2)
plot(ports,thickLam,'LineWidth',2)
plot(ports,blh,'LineWidth',2)
xlabel('Port Number')
ylabel('Boundary Layer Thickness (mm)')
title('Boundary Layer Thickness Turbulent and Laminar')
legend('Turbulent','Laminar','Actual Boundary Layer')
hold off

% ----------------------------------------------------------------------- %
% /// Part 7: AirFoil ///////////////////////////////////////////////////
% ----------------------------------------------------------------------- %
b = 1;
r = 1;
Chord = 88.978; % [mm] - chord length of the airfoil

for k = [1 2 3 4 5 6 7 8] % Reading the files
    AirfoilPressure{b} = readmatrix(sprintf...
        ('GroupData/AirfoilPressure_S301_%d.csv',k));
    b = b+1;
end

for K = linspace(1,12,12) % Concatenating
    for L = [4 3 2 1 8 7 6 5]
        AP(r,:) = mean(AirfoilPressure{1,L}((K*20)-19:K*20,1:28));
        r = r+1;
    end
end

AP9 = [AP(1:8,:);AP(25:32,:);AP(49:56,:);AP(73:80,:)];
AP17 = [AP(9:16,:);AP(33:40,:);AP(57:64,:);AP(81:88,:)];
AP34 = [AP(17:24,:);AP(41:48,:);AP(65:72,:);AP(89:96,:)];

AOA9 = AP9(:,23);
AOA17 = AP17(:,23);
AOA34 = AP34(:,23);

xtop = [0 0.175 0.35 0.7 1.05 1.4 1.75 2.1 2.8 3.5];
xbot = [3.5 2.8 2.1 1.4 1.05 0.7 0.35 0.175 0];
xtot = 25.4.*[0 0.175 0.35 0.7 1.05 1.4 1.75 2.1 2.8 3.5 ...
    2.8 2.1 1.4 1.05 0.7 0.35 0.175];
ytot = 25.4.*[0.14665 0.33075 0.4018 0.476 0.49 0.4774 0.4403 0.38325...
    0.21875 0 0 0 0 0 0.0014 0.0175 0.03885];

% Extrapolating pressure at 'port 11' ----------------------------------- %
m1 = (AP9(:,15)-AP9(:,14))./(19);
m2 = (AP9(:,16)-AP9(:,17))./(19);
b1 = AP9(:,15) - m1.*71.12;
b2 = AP9(:,17) - m2.*71.12;
y1 = b1 + m1.*88.9;
y2 = b2 + m2.*88.9;
Port11Pressure9 = (y1+y2)./2;

m117 = (AP17(:,15)-AP17(:,14))./(19);
m217 = (AP17(:,16)-AP17(:,17))./(19);
b117 = AP17(:,15) - m117.*71.12;
b217 = AP17(:,17) - m217.*71.12;
y117 = b117 + m117.*88.9;
y217 = b217 + m217.*88.9;
Port11Pressure17 = (y117+y217)./2;

m134 = (AP34(:,15)-AP34(:,14))./(19);
m234 = (AP34(:,16)-AP34(:,17))./(19);
b134 = AP34(:,15) - m134.*71.12;
b234 = AP34(:,17) - m234.*71.12;
y134 = b134 + m134.*88.9;
y234 = b234 + m234.*88.9;
Port11Pressure34 = (y134+y234)./2;

PortPressures9 = [AP9(:,7:15),Port11Pressure9,AP9(:,16:22)];
PortPressures17 = [AP17(:,7:15),Port11Pressure17,AP17(:,16:22)];
PortPressures34 = [AP17(:,7:15),Port11Pressure34,AP17(:,16:22)];

% Calculating the Coefficient of Pressure ------------------------------- %
Cp9 = PortPressures9./((AP9(:,3).*AP9(:,4).^2)./2);
Cp17 = PortPressures17./((AP17(:,3).*AP17(:,4).^2)./2);
Cp34 = PortPressures34./((AP34(:,3).*AP34(:,4).^2)./2);

figure
hold on
% % angle of attack = -15^o
% plot(xtot./Chord,Cp9(1,:),'LineWidth',2)
% plot(xtot./Chord,Cp17(1,:),'LineWidth',2)
% plot(xtot./Chord,Cp34(1,:),'LineWidth',2)
% title('Coeeficient of Pressure At -15^o Angle of Attack')

% % angle of attack = -14^o
% plot(xtot./Chord,Cp9(2,:),'LineWidth',2)
% plot(xtot./Chord,Cp17(2,:),'LineWidth',2)
% plot(xtot./Chord,Cp34(2,:),'LineWidth',2)
% title('Coeeficient of Pressure At -14^o Angle of Attack')

% % angle of attack = -13^o
% plot(xtot./Chord,Cp9(3,:),'LineWidth',2)
% plot(xtot./Chord,Cp17(3,:),'LineWidth',2)
% plot(xtot./Chord,Cp34(3,:),'LineWidth',2)
% title('Coeeficient of Pressure At -13^o Angle of Attack')

% % angle of attack = -12^o
% plot(xtot./Chord,Cp9(4,:),'LineWidth',2)
% plot(xtot./Chord,Cp17(4,:),'LineWidth',2)
% plot(xtot./Chord,Cp34(4,:),'LineWidth',2)
% title('Coeeficient of Pressure At -12^o Angle of Attack')

% % angle of attack = -11^o
% plot(xtot./Chord,Cp9(5,:),'LineWidth',2)
% plot(xtot./Chord,Cp17(5,:),'LineWidth',2)
% plot(xtot./Chord,Cp34(5,:),'LineWidth',2)
% title('Coeeficient of Pressure At -11^o Angle of Attack')

% % angle of attack = -10^o
% plot(xtot./Chord,Cp9(6,:),'LineWidth',2)
% plot(xtot./Chord,Cp17(6,:),'LineWidth',2)
% plot(xtot./Chord,Cp34(6,:),'LineWidth',2)
% title('Coeeficient of Pressure At -10^o Angle of Attack')

% % angle of attack = -9^o
% plot(xtot./Chord,Cp9(7,:),'LineWidth',2)
% plot(xtot./Chord,Cp17(7,:),'LineWidth',2)
% plot(xtot./Chord,Cp34(7,:),'LineWidth',2)
% title('Coeeficient of Pressure At -9^o Angle of Attack')

% % angle of attack = -8^o
% plot(xtot./Chord,Cp9(8,:),'LineWidth',2)
% plot(xtot./Chord,Cp17(8,:),'LineWidth',2)
% plot(xtot./Chord,Cp34(8,:),'LineWidth',2)
% title('Coeeficient of Pressure At -8^o Angle of Attack')

% % angle of attack = -7^o
% plot(xtot./Chord,Cp9(9,:),'LineWidth',2)
% plot(xtot./Chord,Cp17(9,:),'LineWidth',2)
% plot(xtot./Chord,Cp34(9,:),'LineWidth',2)
% title('Coeeficient of Pressure At -7^o Angle of Attack')

% % angle of attack = -6^o
% plot(xtot./Chord,Cp9(10,:),'LineWidth',2)
% plot(xtot./Chord,Cp17(10,:),'LineWidth',2)
% plot(xtot./Chord,Cp34(10,:),'LineWidth',2)
% title('Coeeficient of Pressure At -6^o Angle of Attack')

% % angle of attack = -5^o
% plot(xtot./Chord,Cp9(11,:),'LineWidth',2)
% plot(xtot./Chord,Cp17(11,:),'LineWidth',2)
% plot(xtot./Chord,Cp34(11,:),'LineWidth',2)
% title('Coeeficient of Pressure At -5^o Angle of Attack')

% % angle of attack = -4^o
% plot(xtot./Chord,Cp9(12,:),'LineWidth',2)
% plot(xtot./Chord,Cp17(12,:),'LineWidth',2)
% plot(xtot./Chord,Cp34(12,:),'LineWidth',2)
% title('Coeeficient of Pressure At -4^o Angle of Attack')

% % angle of attack = -3^o
% plot(xtot./Chord,Cp9(13,:),'LineWidth',2)
% plot(xtot./Chord,Cp17(13,:),'LineWidth',2)
% plot(xtot./Chord,Cp34(13,:),'LineWidth',2)
% title('Coeeficient of Pressure At -3^o Angle of Attack')

% % angle of attack = -2^o
% plot(xtot./Chord,Cp9(14,:),'LineWidth',2)
% plot(xtot./Chord,Cp17(14,:),'LineWidth',2)
% plot(xtot./Chord,Cp34(14,:),'LineWidth',2)
% title('Coeeficient of Pressure At -2^o Angle of Attack')

% % angle of attack = -1^o
% plot(xtot./Chord,Cp9(15,:),'LineWidth',2)
% plot(xtot./Chord,Cp17(15,:),'LineWidth',2)
% plot(xtot./Chord,Cp34(15,:),'LineWidth',2)
% title('Coeeficient of Pressure At -1^o Angle of Attack')

% angle of attack = 0^o
plot(xtot./Chord,Cp9(16,:),'LineWidth',2)
plot(xtot./Chord,Cp17(16,:),'LineWidth',2)
plot(xtot./Chord,Cp34(16,:),'LineWidth',2)
title('Coeeficient of Pressure At 0^o Angle of Attack')

% % angle of attack = 1^o
% plot(xtot./Chord,Cp9(17,:),'LineWidth',2)
% plot(xtot./Chord,Cp17(17,:),'LineWidth',2)
% plot(xtot./Chord,Cp34(17,:),'LineWidth',2)
% title('Coeeficient of Pressure At 1^o Angle of Attack')

% % angle of attack = 2^o
% plot(xtot./Chord,Cp9(18,:),'LineWidth',2)
% plot(xtot./Chord,Cp17(18,:),'LineWidth',2)
% plot(xtot./Chord,Cp34(18,:),'LineWidth',2)
% title('Coeeficient of Pressure At 2^o Angle of Attack')

% % angle of attack = 3^o
% plot(xtot./Chord,Cp9(19,:),'LineWidth',2)
% plot(xtot./Chord,Cp17(19,:),'LineWidth',2)
% plot(xtot./Chord,Cp34(19,:),'LineWidth',2)
% title('Coeeficient of Pressure At 3^o Angle of Attack')

% % angle of attack = 4^o
% plot(xtot./Chord,Cp9(20,:),'LineWidth',2)
% plot(xtot./Chord,Cp17(20,:),'LineWidth',2)
% plot(xtot./Chord,Cp34(20,:),'LineWidth',2)
% title('Coeeficient of Pressure At 4^o Angle of Attack')

% % angle of attack = 5^o
% plot(xtot./Chord,Cp9(21,:),'LineWidth',2)
% plot(xtot./Chord,Cp17(21,:),'LineWidth',2)
% plot(xtot./Chord,Cp34(21,:),'LineWidth',2)
% title('Coeeficient of Pressure At 5^o Angle of Attack')

% % angle of attack = 6^o
% plot(xtot./Chord,Cp9(22,:),'LineWidth',2)
% plot(xtot./Chord,Cp17(22,:),'LineWidth',2)
% plot(xtot./Chord,Cp34(22,:),'LineWidth',2)
% title('Coeeficient of Pressure At 6^o Angle of Attack')

% % angle of attack = 7^o
% plot(xtot./Chord,Cp9(23,:),'LineWidth',2)
% plot(xtot./Chord,Cp17(23,:),'LineWidth',2)
% plot(xtot./Chord,Cp34(23,:),'LineWidth',2)
% title('Coeeficient of Pressure At 7^o Angle of Attack')

% % angle of attack = 8^o
% plot(xtot./Chord,Cp9(24,:),'LineWidth',2)
% plot(xtot./Chord,Cp17(24,:),'LineWidth',2)
% plot(xtot./Chord,Cp34(24,:),'LineWidth',2)
% title('Coeeficient of Pressure At 8^o Angle of Attack')

% % angle of attack = 9^o
% plot(xtot./Chord,Cp9(25,:),'LineWidth',2)
% plot(xtot./Chord,Cp17(25,:),'LineWidth',2)
% plot(xtot./Chord,Cp34(25,:),'LineWidth',2)
% title('Coeeficient of Pressure At 9^o Angle of Attack')

% % angle of attack = 10^o
% plot(xtot./Chord,Cp9(26,:),'LineWidth',2)
% plot(xtot./Chord,Cp17(26,:),'LineWidth',2)
% plot(xtot./Chord,Cp34(26,:),'LineWidth',2)
% title('Coeeficient of Pressure At 10^o Angle of Attack')

% % angle of attack = 11^o
% plot(xtot./Chord,Cp9(27,:),'LineWidth',2)
% plot(xtot./Chord,Cp17(27,:),'LineWidth',2)
% plot(xtot./Chord,Cp34(27,:),'LineWidth',2)
% title('Coeeficient of Pressure At 11^o Angle of Attack')

% % angle of attack = 12^o
% plot(xtot./Chord,Cp9(28,:),'LineWidth',2)
% plot(xtot./Chord,Cp17(28,:),'LineWidth',2)
% plot(xtot./Chord,Cp34(28,:),'LineWidth',2)
% title('Coeeficient of Pressure At 12^o Angle of Attack')

% % angle of attack = 13^o
% plot(xtot./Chord,Cp9(29,:),'LineWidth',2)
% plot(xtot./Chord,Cp17(29,:),'LineWidth',2)
% plot(xtot./Chord,Cp34(29,:),'LineWidth',2)
% title('Coeeficient of Pressure At 13^o Angle of Attack')

% % angle of attack = 14^o
% plot(xtot./Chord,Cp9(30,:),'LineWidth',2)
% plot(xtot./Chord,Cp17(30,:),'LineWidth',2)
% plot(xtot./Chord,Cp34(30,:),'LineWidth',2)
% title('Coeeficient of Pressure At 14^o Angle of Attack')

% % angle of attack = 15^o
% plot(xtot./Chord,Cp9(31,:),'LineWidth',2)
% plot(xtot./Chord,Cp17(31,:),'LineWidth',2)
% plot(xtot./Chord,Cp34(31,:),'LineWidth',2)
% title('Coeeficient of Pressure At 15^o Angle of Attack')

% % angle of attack = 16^o
% plot(xtot./Chord,Cp9(32,:),'LineWidth',2)
% plot(xtot./Chord,Cp17(32,:),'LineWidth',2)
% plot(xtot./Chord,Cp34(32,:),'LineWidth',2)
% title('Coeeficient of Pressure At 16^o Angle of Attack')

xlabel('Normalized Chord Position (x/c)')
ylabel('Coeeficient of Pressure (C_p)')
legend('~9m/s','~17m/s','~34m/s','Location','SouthEast')
set(gca,'Ydir','reverse');
hold off

% Coefficient of Lift and Coefficent of Drag ---------------------------- %
dx = [4.445 4.445 8.89 8.89 8.89 8.89 8.89 17.78 17.78 -17.78 -17.78...
    -17.78 -8.89 -8.89 -8.89 -4.445 -4.445];
dy = [4.6742 1.8046 1.8847 0.3556 -0.32 -0.9424 -1.4491 -4.1783 -5.5562...
    0 0 0 0 0.0356 0.4089 0.5423 2.73811];

% % alpha = -15
% Cn91 = (Cp9(1,1)+Cp9(1,2)).*(dx(1)/(Chord*2));
% Cn92 = (Cp9(1,2)+Cp9(1,3)).*(dx(2)/(Chord*2));
% Cn93 = (Cp9(1,3)+Cp9(1,4)).*(dx(3)/(Chord*2));
% Cn94 = (Cp9(1,4)+Cp9(1,5)).*(dx(4)/(Chord*2));
% Cn95 = (Cp9(1,5)+Cp9(1,6)).*(dx(5)/(Chord*2));
% Cn96 = (Cp9(1,6)+Cp9(1,7)).*(dx(6)/(Chord*2));
% Cn97 = (Cp9(1,7)+Cp9(1,8)).*(dx(7)/(Chord*2));
% Cn98 = (Cp9(1,8)+Cp9(1,9)).*(dx(8)/(Chord*2));
% Cn99 = (Cp9(1,9)+Cp9(1,10)).*(dx(9)/(Chord*2));
% Cn910 = (Cp9(1,10)+Cp9(1,11)).*(dx(10)/(Chord*2));
% Cn911 = (Cp9(1,11)+Cp9(1,12)).*(dx(11)/(Chord*2));
% Cn912 = (Cp9(1,12)+Cp9(1,13)).*(dx(12)/(Chord*2));
% Cn913 = (Cp9(1,13)+Cp9(1,14)).*(dx(13)/(Chord*2));
% Cn914 = (Cp9(1,14)+Cp9(1,15)).*(dx(14)/(Chord*2));
% Cn915 = (Cp9(1,15)+Cp9(1,16)).*(dx(15)/(Chord*2));
% Cn916 = (Cp9(1,16)+Cp9(1,17)).*(dx(16)/(Chord*2));
% Cn917 = (Cp9(1,17)).*(dx(17)/(Chord*2));

CP9 = [Cp9 zeros(32,1)];
CP17 = [Cp17 zeros(32,1)];
CP34 = [Cp34 zeros(32,1)];

for j = 1:32
    for i = 1:17
        Cn9(j,i) = (CP9(j,i)+CP9(j,i+1)).*(dx(i)/(Chord*2));
    end
end
for k = 1:32
    cn9 = sum(Cn9(k,1:17));
end

for j = 1:32
    for i = 1:17
        Ca9(j,i) = (CP9(j,i)+CP9(j,i+1)).*(dy(i)/(Chord*2));
    end
end
for k = 1:32
    ca9 = sum(Ca9(k,1:17));
end

Cl9 = cn9.*cosd(AOA9) - ca9.*sind(AOA9);
Cd9 = cn9.*sind(AOA9) + ca9.*cosd(AOA9);

% ----------------------------------------------------------------------- %
% NACA Data
nacaAOA = [-8 -4 0 4 8 12 16 20 24 28];
nacaCl = [-0.12 0.18 0.48 0.78 1.08 1.3 1.52 1.58 1.45 1.28];
nacaCd = [0.015 0.015 0.02 0.04 0.08 0.12 0.17 0.24 0.36 0.46];

figure
hold on
plot(AOA9,Cd9)
plot(nacaAOA,nacaCl,'-ko','LineWidth',2,'MarkerSize',10)
plot(nacaAOA,nacaCd,'-ko','LineWidth',2,'MarkerSize',10)
hold off









% dxi = abs(25.4.*[0.175 0.175 0.35 0.35 0.35 0.35 0.35 0.7 0.7 ...
%     -0.7 -0.7 -0.7 -0.35 -0.35 -0.35 -0.175 -0.175]);
% dyi = abs(25.4.*[0.1841 0.07105 0.742 0.014 -0.0126 -0.0371 -0.05705...
%     -0.1645 -0.21875 0 0 0 0 0.0014 0.0161 0.02135 -0.1078]);
% 
% cn91 = ((Cp9(1,1) + Cp9(1,2))/2)*(dxi(1)/Chord);
% cn92 = ((Cp9(1,2) + Cp9(1,3))/2)*(dxi(2)/Chord);
% cn93 = ((Cp9(1,3) + Cp9(1,4))/2)*(dxi(3)/Chord);
% cn94 = ((Cp9(1,4) + Cp9(1,5))/2)*(dxi(4)/Chord);
% cn95 = ((Cp9(1,5) + Cp9(1,6))/2)*(dxi(5)/Chord);
% cn96 = ((Cp9(1,6) + Cp9(1,7))/2)*(dxi(6)/Chord);
% cn97 = ((Cp9(1,7) + Cp9(1,8))/2)*(dxi(7)/Chord);
% cn98 = ((Cp9(1,8) + Cp9(1,9))/2)*(dxi(8)/Chord);
% cn99 = ((Cp9(1,9) + Cp9(1,10))/2)*(dxi(9)/Chord);
% cn910 = ((Cp9(1,10) + Cp9(1,11))/2)*(dxi(10)/Chord);
% cn911 = ((Cp9(1,11) + Cp9(1,12))/2)*(dxi(11)/Chord);
% cn912 = ((Cp9(1,12) + Cp9(1,13))/2)*(dxi(12)/Chord);
% cn913 = ((Cp9(1,13) + Cp9(1,14))/2)*(dxi(13)/Chord);
% cn914 = ((Cp9(1,14) + Cp9(1,15))/2)*(dxi(14)/Chord);
% cn915 = ((Cp9(1,15) + Cp9(1,16))/2)*(dxi(15)/Chord);
% cn916 = ((Cp9(1,16) + Cp9(1,17))/2)*(dxi(16)/Chord);
% cn917 = ((Cp9(1,17))*(dxi(17)/Chord));
% cn9 = sum([cn91 cn92 cn93 cn94 cn95 cn96 cn97 cn98 cn99 cn910 cn911...
%       cn912 cn913 cn914 cn915 cn916 cn917].*(-1));
% 
% ca91 = ((Cp9(1,1) + Cp9(1,2))/2)*(dyi(1)/Chord);
% ca92 = ((Cp9(1,2) + Cp9(1,3))/2)*(dyi(2)/Chord);
% ca93 = ((Cp9(1,3) + Cp9(1,4))/2)*(dyi(3)/Chord);
% ca94 = ((Cp9(1,4) + Cp9(1,5))/2)*(dyi(4)/Chord);
% ca95 = ((Cp9(1,5) + Cp9(1,6))/2)*(dyi(5)/Chord);
% ca96 = ((Cp9(1,6) + Cp9(1,7))/2)*(dyi(6)/Chord);
% ca97 = ((Cp9(1,7) + Cp9(1,8))/2)*(dyi(7)/Chord);
% ca98 = ((Cp9(1,8) + Cp9(1,9))/2)*(dyi(8)/Chord);
% ca99 = ((Cp9(1,9) + Cp9(1,10))/2)*(dyi(9)/Chord);
% ca910 = ((Cp9(1,10) + Cp9(1,11))/2)*(dyi(10)/Chord);
% ca911 = ((Cp9(1,11) + Cp9(1,12))/2)*(dyi(11)/Chord);
% ca912 = ((Cp9(1,12) + Cp9(1,13))/2)*(dyi(12)/Chord);
% ca913 = ((Cp9(1,13) + Cp9(1,14))/2)*(dyi(13)/Chord);
% ca914 = ((Cp9(1,14) + Cp9(1,15))/2)*(dyi(14)/Chord);
% ca915 = ((Cp9(1,15) + Cp9(1,16))/2)*(dyi(15)/Chord);
% ca916 = ((Cp9(1,16) + Cp9(1,17))/2)*(dyi(16)/Chord);
% ca917 = ((Cp9(1,17))/2)*(dyi(17)/Chord);
% ca9 = sum([ca91 ca92 ca93 ca94 ca95 ca96 ca97 ca98 ca99 ca910 ca911...
%       ca912 ca913 ca914 ca915 ca916 ca917]);
% 
% cl9 = cn9.*cosd(AOA9) - ca9.*sind(AOA9);         
% cd9 = cn9.*sind(AOA9) + ca9.*cosd(AOA9); 
% 
% 
% cn171 = ((Cp17(1,1) + Cp17(1,2))/2)*(dxi(1)/Chord);
% cn172 = ((Cp17(1,2) + Cp17(1,3))/2)*(dxi(2)/Chord);
% cn173 = ((Cp17(1,3) + Cp17(1,4))/2)*(dxi(3)/Chord);
% cn174 = ((Cp17(1,4) + Cp17(1,5))/2)*(dxi(4)/Chord);
% cn175 = ((Cp17(1,5) + Cp17(1,6))/2)*(dxi(5)/Chord);
% cn176 = ((Cp17(1,6) + Cp17(1,7))/2)*(dxi(6)/Chord);
% cn177 = ((Cp17(1,7) + Cp17(1,8))/2)*(dxi(7)/Chord);
% cn178 = ((Cp17(1,8) + Cp17(1,9))/2)*(dxi(8)/Chord);
% cn179 = ((Cp17(1,9) + Cp17(1,10))/2)*(dxi(9)/Chord);
% cn1710 = ((Cp17(1,10) + Cp17(1,11))/2)*(dxi(10)/Chord);
% cn1711 = ((Cp17(1,11) + Cp17(1,12))/2)*(dxi(11)/Chord);
% cn1712 = ((Cp17(1,12) + Cp17(1,13))/2)*(dxi(12)/Chord);
% cn1713 = ((Cp17(1,13) + Cp17(1,14))/2)*(dxi(13)/Chord);
% cn1714 = ((Cp17(1,14) + Cp17(1,15))/2)*(dxi(14)/Chord);
% cn1715 = ((Cp17(1,15) + Cp17(1,16))/2)*(dxi(15)/Chord);
% cn1716 = ((Cp17(1,16) + Cp17(1,17))/2)*(dxi(16)/Chord);
% cn1717 = ((Cp17(1,17))*(dxi(17)/Chord));
% cn17 = sum([cn171 cn172 cn173 cn174 cn175 cn176 cn177 cn178 cn179 cn1710...
%     cn1711 cn1712 cn1713 cn1714 cn1715 cn1716 cn1717].*(-1));
% 
% ca171 = ((Cp17(1,1) + Cp17(1,2))/2)*(dyi(1)/Chord);
% ca172 = ((Cp17(1,2) + Cp17(1,3))/2)*(dyi(2)/Chord);
% ca173 = ((Cp17(1,3) + Cp17(1,4))/2)*(dyi(3)/Chord);
% ca174 = ((Cp17(1,4) + Cp17(1,5))/2)*(dyi(4)/Chord);
% ca175 = ((Cp17(1,5) + Cp17(1,6))/2)*(dyi(5)/Chord);
% ca176 = ((Cp17(1,6) + Cp17(1,7))/2)*(dyi(6)/Chord);
% ca177 = ((Cp17(1,7) + Cp17(1,8))/2)*(dyi(7)/Chord);
% ca178 = ((Cp17(1,8) + Cp17(1,9))/2)*(dyi(8)/Chord);
% ca179 = ((Cp17(1,9) + Cp17(1,10))/2)*(dyi(9)/Chord);
% ca1710 = ((Cp17(1,10) + Cp17(1,11))/2)*(dyi(10)/Chord);
% ca1711 = ((Cp17(1,11) + Cp17(1,12))/2)*(dyi(11)/Chord);
% ca1712 = ((Cp17(1,12) + Cp17(1,13))/2)*(dyi(12)/Chord);
% ca1713 = ((Cp17(1,13) + Cp17(1,14))/2)*(dyi(13)/Chord);
% ca1714 = ((Cp17(1,14) + Cp17(1,15))/2)*(dyi(14)/Chord);
% ca1715 = ((Cp17(1,15) + Cp17(1,16))/2)*(dyi(15)/Chord);
% ca1716 = ((Cp17(1,16) + Cp17(1,17))/2)*(dyi(16)/Chord);
% ca1717 = ((Cp17(1,17))/2)*(dyi(17)/Chord);
% ca17 = sum([ca171 ca172 ca173 ca94 ca175 ca176 ca177 ca178 ca179 ca1710...
%     ca1711 ca1712 ca1713 ca1714 ca1715 ca1716 ca1717]);
% 
% cl17 = cn17.*cosd(AOA17) - ca17.*sind(AOA17);         
% cd17 = cn17.*sind(AOA17) + ca17.*cosd(AOA17); 
% 
% 
% cn341 = ((Cp34(1,1) + Cp34(1,2))/2)*(dxi(1)/Chord);
% cn342 = ((Cp34(1,2) + Cp34(1,3))/2)*(dxi(2)/Chord);
% cn343 = ((Cp34(1,3) + Cp34(1,4))/2)*(dxi(3)/Chord);
% cn344 = ((Cp34(1,4) + Cp34(1,5))/2)*(dxi(4)/Chord);
% cn345 = ((Cp34(1,5) + Cp34(1,6))/2)*(dxi(5)/Chord);
% cn346 = ((Cp34(1,6) + Cp34(1,7))/2)*(dxi(6)/Chord);
% cn347 = ((Cp34(1,7) + Cp34(1,8))/2)*(dxi(7)/Chord);
% cn348 = ((Cp34(1,8) + Cp34(1,9))/2)*(dxi(8)/Chord);
% cn349 = ((Cp34(1,9) + Cp34(1,10))/2)*(dxi(9)/Chord);
% cn3410 = ((Cp34(1,10) + Cp34(1,11))/2)*(dxi(10)/Chord);
% cn3411 = ((Cp34(1,11) + Cp34(1,12))/2)*(dxi(11)/Chord);
% cn3412 = ((Cp34(1,12) + Cp34(1,13))/2)*(dxi(12)/Chord);
% cn3413 = ((Cp34(1,13) + Cp34(1,14))/2)*(dxi(13)/Chord);
% cn3414 = ((Cp34(1,14) + Cp34(1,15))/2)*(dxi(14)/Chord);
% cn3415 = ((Cp34(1,15) + Cp34(1,16))/2)*(dxi(15)/Chord);
% cn3416 = ((Cp34(1,16) + Cp34(1,17))/2)*(dxi(16)/Chord);
% cn3417 = ((Cp34(1,17))*(dxi(17)/Chord));
% cn34 = sum([cn341 cn342 cn343 cn344 cn345 cn346 cn347 cn348 cn349 cn3410...
%     cn3411 cn3412 cn3413 cn3414 cn3415 cn3416 cn3417].*(-1));
% 
% ca341 = ((Cp34(1,1) + Cp34(1,2))/2)*(dyi(1)/Chord);
% ca342 = ((Cp34(1,2) + Cp34(1,3))/2)*(dyi(2)/Chord);
% ca343 = ((Cp34(1,3) + Cp34(1,4))/2)*(dyi(3)/Chord);
% ca344 = ((Cp34(1,4) + Cp34(1,5))/2)*(dyi(4)/Chord);
% ca345 = ((Cp34(1,5) + Cp34(1,6))/2)*(dyi(5)/Chord);
% ca346 = ((Cp34(1,6) + Cp34(1,7))/2)*(dyi(6)/Chord);
% ca347 = ((Cp34(1,7) + Cp34(1,8))/2)*(dyi(7)/Chord);
% ca348 = ((Cp34(1,8) + Cp34(1,9))/2)*(dyi(8)/Chord);
% ca349 = ((Cp34(1,9) + Cp34(1,10))/2)*(dyi(9)/Chord);
% ca3410 = ((Cp34(1,10) + Cp34(1,11))/2)*(dyi(10)/Chord);
% ca3411 = ((Cp34(1,11) + Cp34(1,12))/2)*(dyi(11)/Chord);
% ca3412 = ((Cp34(1,12) + Cp34(1,13))/2)*(dyi(12)/Chord);
% ca3413 = ((Cp34(1,13) + Cp34(1,14))/2)*(dyi(13)/Chord);
% ca3414 = ((Cp34(1,14) + Cp34(1,15))/2)*(dyi(14)/Chord);
% ca3415 = ((Cp34(1,15) + Cp34(1,16))/2)*(dyi(15)/Chord);
% ca3416 = ((Cp34(1,16) + Cp34(1,17))/2)*(dyi(16)/Chord);
% ca3417 = ((Cp34(1,17))/2)*(dyi(17)/Chord);
% ca34 = sum([ca341 ca342 ca343 ca344 ca345 ca346 ca347 ca348 ca349 ca3410...
%     ca3411 ca3412 ca3413 ca3414 ca3415 ca3416 ca3417]);
% 
% cl34 = cn34.*cosd(AOA34) - ca34.*sind(AOA34);         
% cd34 = cn34.*sind(AOA34) + ca34.*cosd(AOA34);
% 
% figure
% hold on
% plot(AOA9,cl9,'LineWidth',2)
% plot(AOA17,cl17,'LineWidth',2)
% plot(AOA34,cl34,'LineWidth',2)
% xlabel('Angle of Attack')
% ylabel('Coefficient of Lift (C_l)')
% legend('~9m/s','~17/m/s','~34m/s','Location','NorthWest')
% hold off
% 
% figure
% hold on
% plot(AOA9,cd9,'LineWidth',2)
% plot(AOA17,cd17,'LineWidth',2)
% plot(AOA34,cd34,'LineWidth',2)
% xlabel('Angle of Attack')
% ylabel('Coefficient of Drag (C_d)')
% legend('~9m/s','~17/m/s','~34m/s','Location','NorthWest')
% hold off












% cp9v = 1 - (AP9(:,4)./Vfreestream).^2;
% cp17v = 1 - (AP17(:,4)./Vfreestream).^2;
% cp34v = 1 - (AP34(:,4)./Vfreestream).^2;
% 
% cn9v = -sum(cp9v.*dxi./Chord);
% cn17v = -sum(cp17v).*dxi./Chord;
% cn34v = -sum(cp34v).*dxi./Chord;
% 
% ca9v = sum(cp9v.*dyi./Chord);
% ca17v = sum(cp17v).*dyi./Chord;
% ca34v = sum(cp34v).*dyi./Chord;
% 
% cl9v = cn9v.*cosd(AOA9) - ca9v.*sind(AOA9);         
% cd9v = cn9v.*sind(AOA9) + ca9v.*cosd(AOA9);
% cl17v = cn17v.*cosd(AOA17) - ca17v.*sind(AOA17);         
% cd17v = cn17v.*sind(AOA17) + ca17v.*cosd(AOA17);
% cl34v = cn34v.*cosd(AOA34) - ca34v.*sind(AOA34);         
% cd34v = cn34v.*sind(AOA34) + ca34v.*cosd(AOA34);
% 
% figure
% plot(AOA9,cl9v)






% m1 = (AP(:,15)-AP(:,14))./(19);
% m2 = (AP(:,16)-AP(:,17))./(19);
% b1 = AP(:,15) - m1.*71.12;
% b2 = AP(:,17) - m2.*71.12;
% y1 = b1 + m1.*88.9;
% y2 = b2 + m2.*88.9;
% Port11Pressure = (y1+y2)./2;
% AOA = AP(:,23);
% PortPressures = [AP(:,7:15),Port11Pressure,AP(:,16:22)];
% Cp = PortPressures./((AP(:,3).*AP(:,4).^2)./2);
% 
% cn1 = ((Cp(1,1) + Cp(1,2))/2)*(dxi(1)/Chord);
% cn2 = ((Cp(1,2) + Cp(1,3))/2)*(dxi(2)/Chord);
% cn3 = ((Cp(1,3) + Cp(1,4))/2)*(dxi(3)/Chord);
% cn4 = ((Cp(1,4) + Cp(1,5))/2)*(dxi(4)/Chord);
% cn5 = ((Cp(1,5) + Cp(1,6))/2)*(dxi(5)/Chord);
% cn6 = ((Cp(1,6) + Cp(1,7))/2)*(dxi(6)/Chord);
% cn7 = ((Cp(1,7) + Cp(1,8))/2)*(dxi(7)/Chord);
% cn8 = ((Cp(1,8) + Cp(1,9))/2)*(dxi(8)/Chord);
% cn9 = ((Cp(1,9) + Cp(1,10))/2)*(dxi(9)/Chord);
% cn10 = ((Cp(1,10) + Cp(1,11))/2)*(dxi(10)/Chord);
% cn11 = ((Cp(1,11) + Cp(1,12))/2)*(dxi(11)/Chord);
% cn12 = ((Cp(1,12) + Cp(1,13))/2)*(dxi(12)/Chord);
% cn13 = ((Cp(1,13) + Cp(1,14))/2)*(dxi(13)/Chord);
% cn14 = ((Cp(1,14) + Cp(1,15))/2)*(dxi(14)/Chord);
% cn15 = ((Cp(1,15) + Cp(1,16))/2)*(dxi(15)/Chord);
% cn16 = ((Cp(1,16) + Cp(1,17))/2)*(dxi(16)/Chord);
% cn17 = ((Cp(1,17))*(dxi(17)/Chord));
% cn = sum([cn1 cn2 cn3 cn4 cn5 cn6 cn7 cn8 cn9 cn10 cn11 cn12 cn13 cn14...
%     cn15 cn16 cn17].*(-1));
% 
% ca1 = ((Cp(1,1) + Cp(1,2))/2)*(dyi(1)/Chord);
% ca2 = ((Cp(1,2) + Cp(1,3))/2)*(dyi(2)/Chord);
% ca3 = ((Cp(1,3) + Cp(1,4))/2)*(dyi(3)/Chord);
% ca4 = ((Cp(1,4) + Cp(1,5))/2)*(dyi(4)/Chord);
% ca5 = ((Cp(1,5) + Cp(1,6))/2)*(dyi(5)/Chord);
% ca6 = ((Cp(1,6) + Cp(1,7))/2)*(dyi(6)/Chord);
% ca7 = ((Cp(1,7) + Cp(1,8))/2)*(dyi(7)/Chord);
% ca8 = ((Cp(1,8) + Cp(1,9))/2)*(dyi(8)/Chord);
% ca9 = ((Cp(1,9) + Cp(1,10))/2)*(dyi(9)/Chord);
% ca10 = ((Cp(1,10) + Cp(1,11))/2)*(dyi(10)/Chord);
% ca11 = ((Cp(1,11) + Cp(1,12))/2)*(dyi(11)/Chord);
% ca12 = ((Cp(1,12) + Cp(1,13))/2)*(dyi(12)/Chord);
% ca13 = ((Cp(1,13) + Cp(1,14))/2)*(dyi(13)/Chord);
% ca14 = ((Cp(1,14) + Cp(1,15))/2)*(dyi(14)/Chord);
% ca15 = ((Cp(1,15) + Cp(1,16))/2)*(dyi(15)/Chord);
% ca16 = ((Cp(1,16) + Cp(1,17))/2)*(dyi(16)/Chord);
% ca17 = ((Cp(1,17))/2)*(dyi(17)/Chord);
% ca = sum([ca1 ca2 ca3 ca4 ca5 ca6 ca7 ca8 ca9 ca10 ca11 ca12 ca13 ca14...
%     ca15 ca16 ca17]);
% 
% cosalpha = cosd(AOA); 
% sinalpha = sind(AOA);
% 
% cl = cn.*cosalpha - ca.*sinalpha;         
% cd = cn.*sinalpha + ca.*cosalpha; 
% 
% Cn = -sum(Cp).*(abs(dxi)./(2*Chord));            
% Ca = sum(Cp).*(abs(dyi)./(2*Chord));
% 
% Cl = Cn.*cosd(AOA) - Ca.*sind(AOA);     
% Cd = Cn.*sind(AOA) + Ca.*cosd(AOA);
% 
% figure
% hold on
% % plot(AOA,cn1.*cosalpha - ca1.*sinalpha)
% plot(AOA,Cd)
% xlabel('Angle of Attack')
% % legend('Coefficient of Lift','Coefficient of Drag')
% hold off

