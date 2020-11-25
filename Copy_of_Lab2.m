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
% Uncertainty

% Pitot-Static Tube
dVdP = (-(R.*diffP.*T_atm)./...
    (((2^0.5).*((R.*diffP.*T_atm)./P_atm)).*P_atm.^2));
dVdT = ((diffP.*R)./(((2^0.5).*P_atm).*(diffP.*R.*T_atm)));
dVdp = ((R.*T_atm)./(((2^0.5).*P_atm).*(diffP.*R.*T_atm)));
deltaP = 3450;      % - [Pa]
deltaT = 0.25;      % - [K]
deltap = 68.9476;   % - [Pa]
uVP = sqrt((dVdP.*deltaP).^2 + (dVdT.*deltaT).^2 + (dVdp.*deltap).^2);
errP = uVP.*ones(size(airSpdP));

figure 
errorbar(voltP,airSpdP,errP,'LineWidth',2)
xlabel('Input Voltage (V)');
ylabel('Resulting Airspeed (m/s)');
title('Airspeed for Pitot-Static Probe');

% Venturi Tube

dVdTV = (sqrt((2.*diffPV.*R)./(TatmV.*PatmV.*(1-(A2/A1)^2)))./2);
dVdPV = (-(sqrt((2.*R.*TatmV.*diffPV)./(1-(A2/A1)^2)).*(PatmV.^(-3/2))));
dVdpV = (sqrt((2.*R.*TatmV)./(diffPV.*PatmV.*(1-(A2/A1)^2)))./2);
dPV = 3450; 
dTV = 0.25;     
dpV = 68.9476;     
uVV = sqrt((dVdPV.*dPV).^2 + (dVdTV.*dTV).^2 + (dVdpV.*dpV).^2);
errV = uVV.*ones(size(airSpdV));
E = errV.*ones(size(airSpdV));

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

rho1 = mean(P1Patm./(R.*P1T)); rho2 = mean(P2Patm./(R.*P2T));
rho3 = mean(P3Patm./(R.*P3T)); rho4 = mean(P4Patm./(R.*P4T));
rho5 = mean(P5Patm./(R.*P5T)); rho6 = mean(P6Patm./(R.*P6T));
rho7 = mean(P7Patm./(R.*P7T)); rho8 = mean(P8Patm./(R.*P8T));
rho9 = mean(P9Patm./(R.*P9T)); rho10 = mean(P10Patm./(R.*P10T));
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
blh = [Pyaxis(6,1) Pyaxis(6,2) Pyaxis(6,3) Pyaxis(6,4) Pyaxis(8,5)...
    Pyaxis(9,6) Pyaxis(9,7) Pyaxis(10,8) Pyaxis(9,9) Pyaxis(9,10)...
    Pyaxis(10,11)];

% Determine if boundary layer is laminar or turbulent!!!!!

mu = 1.7894*(10^-5);

Rex1 =((P1d*velP1(end)*rho)/mu)/1000; Rex2 =((P2d*velP2(end)*rho)/mu)/1000;
Rex3 =((P3d*velP3(end)*rho)/mu)/1000; Rex4 =((P4d*velP4(end)*rho)/mu)/1000;
Rex5 =((P5d*velP5(end)*rho)/mu)/1000; Rex6 =((P6d*velP6(end)*rho)/mu)/1000;
Rex7 =((P7d*velP7(end)*rho)/mu)/1000; Rex8 =((P8d*velP8(end)*rho)/mu)/1000;
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

for K = linspace(1,4,4) % Concatenating based on angle of attack
    for L = [4 3 2 1 8 7 6 5]
        AP(r,:) = mean(AirfoilPressure{1,L}((K*60)-59:K*60,1:28));
        r = r+1;
    end
end

AOA = AP(:,23); % Angles of Attack
ScanivalvePorts = AP(:,7:22); % Scanivlave Matrix
xtop = [0 0.175 0.35 0.7 1.05 1.4 1.75 2.1 2.8 3.5];
xbot = [3.5 2.8 2.1 1.4 1.05 0.7 0.35 0.175 0];
xtot = 25.4.*[0 0.175 0.35 0.7 1.05 1.4 1.75 2.1 2.8 3.5 ...
    2.8 2.1 1.4 1.05 0.7 0.35 0.175];
ytot = [0.14665 0.33075 0.4018 0.476 0.49 0.4774 0.4403 0.38325...
    0.21875 0 0 0 0 0 0.0014 0.0175 0.03885];

% Extrapolating pressure at 'port 11' ----------------------------------- %
m1 = (AP(:,15)-AP(:,14))./(19);
m2 = (AP(:,16)-AP(:,17))./(19);
b1 = AP(:,15) - m1.*71.12;
b2 = AP(:,17) - m2.*71.12;
y1 = b1 + m1.*88.9;
y2 = b2 + m2.*88.9;
Port11Pressure = (y1+y2)./2;

PortPressures = [AP(:,7:15),Port11Pressure,AP(:,16:22)];

% Calculating the Coefficient of Pressure ------------------------------- %
% Cp = (AP(:,7:22)./((AP(:,3).*Vfreestream)./2));
Cp = PortPressures./((AP(:,3).*AP(:,4).^2)./2);
Cp1 = (1-((AP(:,4).^2)./(Vfreestream).^2));

xc = linspace(1,(3.5*25.4)/Chord,17);

figure
hold on
% plot(xtot./Chord,Cp(1,:),'LineWidth',2)
% plot(xtot./Chord,Cp(2,:),'LineWidth',2)
% plot(xtot./Chord,Cp(3,:),'LineWidth',2)
% plot(xtot./Chord,Cp(4,:),'LineWidth',2)
% plot(xtot./Chord,Cp(5,:),'LineWidth',2)
% % plot(xtot./Chord,Cp(6,:),'LineWidth',2)
% % plot(xtot./Chord,Cp(7,:),'LineWidth',2)
% % plot(xtot./Chord,Cp(8,:),'LineWidth',2)
% % plot(xtot./Chord,Cp(9,:),'LineWidth',2)
% % plot(xtot./Chord,Cp(10,:),'LineWidth',2)
% % plot(xtot./Chord,Cp(11,:),'LineWidth',2)
% % plot(xtot./Chord,Cp(12,:),'LineWidth',2)
% % plot(xtot./Chord,Cp(13,:),'LineWidth',2)
% % plot(xtot./Chord,Cp(14,:),'LineWidth',2)
plot(xtot./Chord,Cp(15,:),'LineWidth',2)
plot(xtot./Chord,Cp(16,:),'LineWidth',2)
plot(xtot./Chord,Cp(17,:),'LineWidth',2)
% % plot(xtot./Chord,Cp(18,:),'LineWidth',2)
% % plot(xtot./Chord,Cp(19,:),'LineWidth',2)
% % plot(xtot./Chord,Cp(20,:),'LineWidth',2)
% % plot(xtot./Chord,Cp(21,:),'LineWidth',2)
% % plot(xtot./Chord,Cp(22,:),'LineWidth',2)
% % plot(xtot./Chord,Cp(23,:),'LineWidth',2)
% % plot(xtot./Chord,Cp(24,:),'LineWidth',2)
% % plot(xtot./Chord,Cp(25,:),'LineWidth',2)
% % plot(xtot./Chord,Cp(26,:),'LineWidth',2)
% % plot(xtot./Chord,Cp(27,:),'LineWidth',2)
% % plot(xtot./Chord,Cp(28,:),'LineWidth',2)
% plot(xtot./Chord,Cp(29,:),'LineWidth',2)
% plot(xtot./Chord,Cp(30,:),'LineWidth',2)
% plot(xtot./Chord,Cp(31,:),'LineWidth',2)
% plot(xtot./Chord,Cp(32,:),'LineWidth',2)
xlabel('Normalized Chord Position (x/c)')
ylabel('Coeeficient of Pressure (C_p)')
title('Coeeficient of Pressure At -15^o Angle of Attack')
set(gca,'Ydir','reverse');
hold off



% Coefficient of Lift and Coefficent of Drag ---------------------------- %
dxi = abs(25.4.*[0.175 0.175 0.35 0.35 0.35 0.35 0.35 0.7 0.7 ...
    -0.7 -0.7 -0.7 -0.35 -0.35 -0.35 -0.175 -0.175]);
dyi = abs(25.4.*[0.1841 0.07105 0.742 0.014 -0.0126 -0.0371 -0.05705 -0.1645 -0.21875...
    0 0 0 0 0.0014 0.0161 0.02135 -0.1078]);

Cn = -sum(Cp).*(dxi./(2*Chord));
Ca = sum(Cp).*(dyi./(2*Chord));

cn1 = ((Cp(1,1) + Cp(1,2))/2)*(dxi(1)/Chord);
cn2 = ((Cp(1,2) + Cp(1,3))/2)*(dxi(2)/Chord);
cn3 = ((Cp(1,3) + Cp(1,4))/2)*(dxi(3)/Chord);
cn4 = ((Cp(1,4) + Cp(1,5))/2)*(dxi(4)/Chord);
cn5 = ((Cp(1,5) + Cp(1,6))/2)*(dxi(5)/Chord);
cn6 = ((Cp(1,6) + Cp(1,7))/2)*(dxi(6)/Chord);
cn7 = ((Cp(1,7) + Cp(1,8))/2)*(dxi(7)/Chord);
cn8 = ((Cp(1,8) + Cp(1,9))/2)*(dxi(8)/Chord);
cn9 = ((Cp(1,9) + Cp(1,10))/2)*(dxi(9)/Chord);
cn10 = ((Cp(1,10) + Cp(1,11))/2)*(dxi(10)/Chord);
cn11 = ((Cp(1,11) + Cp(1,12))/2)*(dxi(11)/Chord);
cn12 = ((Cp(1,12) + Cp(1,13))/2)*(dxi(12)/Chord);
cn13 = ((Cp(1,13) + Cp(1,14))/2)*(dxi(13)/Chord);
cn14 = ((Cp(1,14) + Cp(1,15))/2)*(dxi(14)/Chord);
cn15 = ((Cp(1,15) + Cp(1,16))/2)*(dxi(15)/Chord);
cn16 = ((Cp(1,16) + Cp(1,17))/2)*(dxi(16)/Chord);
cn17 = ((Cp(1,17))*(dxi(17)/Chord));
cn = sum([cn1 cn2 cn3 cn4 cn5 cn6 cn7 cn8 cn9 cn10 cn11 cn12 cn13 cn14...
    cn15 cn16 cn17].*(-1));

ca1 = ((Cp(1,1) + Cp(1,2))/2)*(dyi(1)/Chord);
ca2 = ((Cp(1,2) + Cp(1,3))/2)*(dyi(2)/Chord);
ca3 = ((Cp(1,3) + Cp(1,4))/2)*(dyi(3)/Chord);
ca4 = ((Cp(1,4) + Cp(1,5))/2)*(dyi(4)/Chord);
ca5 = ((Cp(1,5) + Cp(1,6))/2)*(dyi(5)/Chord);
ca6 = ((Cp(1,6) + Cp(1,7))/2)*(dyi(6)/Chord);
ca7 = ((Cp(1,7) + Cp(1,8))/2)*(dyi(7)/Chord);
ca8 = ((Cp(1,8) + Cp(1,9))/2)*(dyi(8)/Chord);
ca9 = ((Cp(1,9) + Cp(1,10))/2)*(dyi(9)/Chord);
ca10 = ((Cp(1,10) + Cp(1,11))/2)*(dyi(10)/Chord);
ca11 = ((Cp(1,11) + Cp(1,12))/2)*(dyi(11)/Chord);
ca12 = ((Cp(1,12) + Cp(1,13))/2)*(dyi(12)/Chord);
ca13 = ((Cp(1,13) + Cp(1,14))/2)*(dyi(13)/Chord);
ca14 = ((Cp(1,14) + Cp(1,15))/2)*(dyi(14)/Chord);
ca15 = ((Cp(1,15) + Cp(1,16))/2)*(dyi(15)/Chord);
ca16 = ((Cp(1,16) + Cp(1,17))/2)*(dyi(16)/Chord);
ca17 = ((Cp(1,17))/2)*(dyi(17)/Chord);
ca = sum([ca1 ca2 ca3 ca4 ca5 ca6 ca7 ca8 ca9 ca10 ca11 ca12 ca13 ca14...
    ca15 ca16 ca17]);

cosalpha = cosd(AOA); 
sinalpha = sind(AOA);

cl = cn.*cosalpha - ca.*sinalpha;         
cd = abs(cn.*sinalpha + ca.*cosalpha); 

Cn = -sum(Cp).*(abs(dxi)./(2*Chord));            
Ca = sum(Cp).*(abs(dyi)./(2*Chord));

Cl = Cn.*cosd(AOA) - Ca.*sind(AOA);     
Cd = Cn.*sind(AOA) + Ca.*cosd(AOA);

figure
hold on
% plot(AOA,cn1.*cosalpha - ca1.*sinalpha)
plot(AOA,cd)
xlabel('Angle of Attack')
% legend('Coefficient of Lift','Coefficient of Drag')
hold off

