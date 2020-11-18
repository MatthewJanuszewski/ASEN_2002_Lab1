% Alternate Read-in

clc
clear all
close all

% Port 1

Port_1A = readtable('BoundaryLayer_S301_1.csv');
Port_1B = readtable('BoundaryLayer_S301_2.csv');
Port_1A = table2array(Port_1A);
Port_1B = table2array(Port_1B);
a = 1;
b = 500;
c = 1;
d = 1000;
for i = 1:12
    Port_1(c:d,:) = [Port_1B(a:b,:); Port_1A(a:b,:)];
    a = a + 500;
    b = b + 500;
    c = c + 1000;
    d = d + 1000;
end

Port_1_P_atm = mean(Port_1(:,1)); % P_atm average [Pa]
Port_1_T_atm = mean(Port_1(:,2)); % T_atm average [K]
R = 287;
rho = (Port_1_P_atm/(R*Port_1_T_atm)); % density [kg/m^3]
Port_1_Aux_Diff_Pressure = Port_1(:,4); % [Pa]
Port_1_Velocity_Probe = sqrt(2.*Port_1_Aux_Diff_Pressure/rho); % [m/s]
Port_1_Free_Stream_Pressure = mean(Port_1(11001:12000,3)); % [Pa]
Port_1_Free_Stream_Velocity = sqrt(2*Port_1_Free_Stream_Pressure/rho); % [m/s]
Port_1_Boundary_Layer_Velocity = 0.95*Port_1_Free_Stream_Velocity; % [m/s]
Port_1_y_axis = abs(Port_1(:,6)); % vertical height[mm]

x = (Port_1_Velocity_Probe > Port_1_Boundary_Layer_Velocity);
Port_1_Boundary_Layer_Height = Port_1_y_axis(x);
Boundary_Layer_Height(1) = Port_1_Boundary_Layer_Height(1); % [mm]

e = 1;
f = 500;
for j = 1:22
    y_axis(j) = mean(Port_1_y_axis(e:f));
    e = e + 500;
    f = f + 500;
end

y_axis = y_axis';

e = 1;
f = 500;
for k = 1:22
    Velocity_Probe(k) = mean(Port_1_Velocity_Probe(e:f));
    e = e + 500;
    f = f + 500;
end

Velocity_Probe = Velocity_Probe';
f1 = figure;
plot(Velocity_Probe, y_axis);
title('Port 1')
xlabel('Velocity [m/s]')
ylabel('Height [mm]')

% Port 2

Port_2A = readtable('BoundaryLayer_S302_1.csv');
Port_2B = readtable('BoundaryLayer_S302_2.csv');
Port_2A = table2array(Port_2A);
Port_2B = table2array(Port_2B);
a = 1;
b = 500;
c = 1;
d = 1000;
for i = 1:12
    Port_2(c:d,:) = [Port_2B(a:b,:); Port_2A(a:b,:)];
    a = a + 500;
    b = b + 500;
    c = c + 1000;
    d = d + 1000;
end

Port_2_P_atm = mean(Port_2(:,1)); % P_atm average [Pa]
Port_2_T_atm = mean(Port_2(:,2)); % T_atm average [K]
R = 287;
rho = (Port_2_P_atm/(R*Port_2_T_atm)); % density [kg/m^3]
Port_2_Aux_Diff_Pressure = Port_2(:,4); % [Pa]
Port_2_Velocity_Probe = sqrt(2.*Port_2_Aux_Diff_Pressure/rho); % [m/s]
Port_2_Free_Stream_Pressure = mean(Port_2(11001:12000,3)); % [Pa]
Port_2_Free_Stream_Velocity = sqrt(2*Port_2_Free_Stream_Pressure/rho); % [m/s]
Port_2_Boundary_Layer_Velocity = 0.95*Port_2_Free_Stream_Velocity; % [m/s]
Port_2_y_axis = abs(Port_2(:,6)); % vertical height[mm]

x = (Port_2_Velocity_Probe > Port_2_Boundary_Layer_Velocity);
Port_2_Boundary_Layer_Height = Port_2_y_axis(x);
Boundary_Layer_Height(2) = Port_2_Boundary_Layer_Height(1); % [mm]

e = 1;
f = 500;
for j = 1:22
    y_axis(j) = mean(Port_2_y_axis(e:f));
    e = e + 500;
    f = f + 500;
end

y_axis = y_axis';

e = 1;
f = 500;
for k = 1:22
    Velocity_Probe(k) = mean(Port_2_Velocity_Probe(e:f));
    e = e + 500;
    f = f + 500;
end

Velocity_Probe = Velocity_Probe';
f2 = figure;
plot(Velocity_Probe, y_axis);
title('Port 2')
xlabel('Velocity [m/s]')
ylabel('Height [mm]')

% Port 3

Port_3A = readtable('BoundaryLayer_S303_1.csv');
Port_3B = readtable('BoundaryLayer_S303_2.csv');
Port_3A = table2array(Port_3A);
Port_3B = table2array(Port_3B);
a = 1;
b = 500;
c = 1;
d = 1000;
for i = 1:12
    Port_3(c:d,:) = [Port_3B(a:b,:); Port_3A(a:b,:)];
    a = a + 500;
    b = b + 500;
    c = c + 1000;
    d = d + 1000;
end

Port_3_P_atm = mean(Port_3(:,1)); % P_atm average [Pa]
Port_3_T_atm = mean(Port_3(:,2)); % T_atm average [K]
R = 287;
rho = (Port_3_P_atm/(R*Port_3_T_atm)); % density [kg/m^3]
Port_3_Aux_Diff_Pressure = Port_3(:,4); % [Pa]
Port_3_Velocity_Probe = sqrt(2.*Port_3_Aux_Diff_Pressure/rho); % [m/s]
Port_3_Free_Stream_Pressure = mean(Port_3(11001:12000,3)); % [Pa]
Port_3_Free_Stream_Velocity = sqrt(2*Port_3_Free_Stream_Pressure/rho); % [m/s]
Port_3_Boundary_Layer_Velocity = 0.95*Port_3_Free_Stream_Velocity; % [m/s]
Port_3_y_axis = abs(Port_3(:,6)); % vertical height[mm]

x = (Port_3_Velocity_Probe > Port_3_Boundary_Layer_Velocity);
Port_3_Boundary_Layer_Height = Port_3_y_axis(x);
Boundary_Layer_Height(3) = Port_3_Boundary_Layer_Height(1); % [mm]

e = 1;
f = 500;
for j = 1:22
    y_axis(j) = mean(Port_3_y_axis(e:f));
    e = e + 500;
    f = f + 500;
end

y_axis = y_axis';

e = 1;
f = 500;
for k = 1:22
    Velocity_Probe(k) = mean(Port_3_Velocity_Probe(e:f));
    e = e + 500;
    f = f + 500;
end

Velocity_Probe = Velocity_Probe';
f3 = figure;
plot(Velocity_Probe, y_axis);
title('Port 3')
xlabel('Velocity [m/s]')
ylabel('Height [mm]')

% Port 4

Port_4A = readtable('BoundaryLayer_S301_3.csv');
Port_4B = readtable('BoundaryLayer_S301_4.csv');
Port_4A = table2array(Port_4A);
Port_4B = table2array(Port_4B);
a = 1;
b = 500;
c = 1;
d = 1000;
for i = 1:12
    Port_4(c:d,:) = [Port_4B(a:b,:); Port_4A(a:b,:)];
    a = a + 500;
    b = b + 500;
    c = c + 1000;
    d = d + 1000;
end

Port_4_P_atm = mean(Port_4(:,1)); % P_atm average [Pa]
Port_4_T_atm = mean(Port_4(:,2)); % T_atm average [K]
R = 287;
rho = (Port_4_P_atm/(R*Port_4_T_atm)); % density [kg/m^3]
Port_4_Aux_Diff_Pressure = Port_4(:,4); % [Pa]
Port_4_Velocity_Probe = sqrt(2.*Port_4_Aux_Diff_Pressure/rho); % [m/s]
Port_4_Free_Stream_Pressure = mean(Port_4(11001:12000,3)); % [Pa]
Port_4_Free_Stream_Velocity = sqrt(2*Port_4_Free_Stream_Pressure/rho); % [m/s]
Port_4_Boundary_Layer_Velocity = 0.95*Port_4_Free_Stream_Velocity; % [m/s]
Port_4_y_axis = abs(Port_4(:,6)); % vertical height[mm]

x = (Port_4_Velocity_Probe > Port_4_Boundary_Layer_Velocity);
Port_4_Boundary_Layer_Height = Port_4_y_axis(x);
Boundary_Layer_Height(4) = Port_4_Boundary_Layer_Height(1); % [mm]

e = 1;
f = 500;
for j = 1:22
    y_axis(j) = mean(Port_4_y_axis(e:f));
    e = e + 500;
    f = f + 500;
end

y_axis = y_axis';

e = 1;
f = 500;
for k = 1:22
    Velocity_Probe(k) = mean(Port_4_Velocity_Probe(e:f));
    e = e + 500;
    f = f + 500;
end

Velocity_Probe = Velocity_Probe';
f4 = figure;
plot(Velocity_Probe, y_axis);
title('Port 4')
xlabel('Velocity [m/s]')
ylabel('Height [mm]')

% Port 5

Port_5A = readtable('BoundaryLayer_S302_3.csv');
Port_5B = readtable('BoundaryLayer_S302_4.csv');
Port_5A = table2array(Port_5A);
Port_5B = table2array(Port_5B);
a = 1;
b = 500;
c = 1;
d = 1000;
for i = 1:12
    Port_5(c:d,:) = [Port_5B(a:b,:); Port_5A(a:b,:)];
    a = a + 500;
    b = b + 500;
    c = c + 1000;
    d = d + 1000;
end

Port_5_P_atm = mean(Port_5(:,1)); % P_atm average [Pa]
Port_5_T_atm = mean(Port_5(:,2)); % T_atm average [K]
R = 287;
rho = (Port_5_P_atm/(R*Port_5_T_atm)); % density [kg/m^3]
Port_5_Aux_Diff_Pressure = Port_5(:,4); % [Pa]
Port_5_Velocity_Probe = sqrt(2.*Port_5_Aux_Diff_Pressure/rho); % [m/s]
Port_5_Free_Stream_Pressure = mean(Port_5(11001:12000,3)); % [Pa]
Port_5_Free_Stream_Velocity = sqrt(2*Port_5_Free_Stream_Pressure/rho); % [m/s]
Port_5_Boundary_Layer_Velocity = 0.95*Port_5_Free_Stream_Velocity; % [m/s]
Port_5_y_axis = abs(Port_5(:,6)); % vertical height[mm]

x = (Port_5_Velocity_Probe > Port_5_Boundary_Layer_Velocity);
Port_5_Boundary_Layer_Height = Port_5_y_axis(x);
Boundary_Layer_Height(5) = Port_5_Boundary_Layer_Height(1); % [mm]

e = 1;
f = 500;
for j = 1:22
    y_axis(j) = mean(Port_5_y_axis(e:f));
    e = e + 500;
    f = f + 500;
end

y_axis = y_axis';

e = 1;
f = 500;
for k = 1:22
    Velocity_Probe(k) = mean(Port_5_Velocity_Probe(e:f));
    e = e + 500;
    f = f + 500;
end

Velocity_Probe = Velocity_Probe';
f5 = figure;
plot(Velocity_Probe, y_axis);
title('Port 5')
xlabel('Velocity [m/s]')
ylabel('Height [mm]')

% Port 6

Port_6A = readtable('BoundaryLayer_S301_6.csv');
Port_6B = readtable('BoundaryLayer_S301_9.csv');
Port_6A = table2array(Port_6A);
Port_6B = table2array(Port_6B);
a = 1;
b = 500;
c = 1;
d = 1000;
for i = 1:12
    Port_6(c:d,:) = [Port_6A(a:b,:); Port_6B(a:b,:)];
    a = a + 500;
    b = b + 500;
    c = c + 1000;
    d = d + 1000;
end

Port_6_P_atm = mean(Port_6(:,1)); % P_atm average [Pa]
Port_6_T_atm = mean(Port_6(:,2)); % T_atm average [K]
R = 287;
rho = (Port_6_P_atm/(R*Port_6_T_atm)); % density [kg/m^3]
Port_6_Aux_Diff_Pressure = Port_6(:,4); % [Pa]
Port_6_Velocity_Probe = sqrt(2.*Port_6_Aux_Diff_Pressure/rho); % [m/s]
Port_6_Free_Stream_Pressure = mean(Port_6(11001:12000,3)); % [Pa]
Port_6_Free_Stream_Velocity = sqrt(2*Port_6_Free_Stream_Pressure/rho); % [m/s]
Port_6_Boundary_Layer_Velocity = 0.95*Port_6_Free_Stream_Velocity; % [m/s]
Port_6_y_axis = abs(Port_6(:,6)); % vertical height[mm]

x = (Port_6_Velocity_Probe > Port_6_Boundary_Layer_Velocity);
Port_6_Boundary_Layer_Height = Port_6_y_axis(x);
Boundary_Layer_Height(6) = Port_6_Boundary_Layer_Height(1); % [mm]

e = 1;
f = 500;
for j = 1:22
    y_axis(j) = mean(Port_6_y_axis(e:f));
    e = e + 500;
    f = f + 500;
end

y_axis = y_axis';

e = 1;
f = 500;
for k = 1:22
    Velocity_Probe(k) = mean(Port_6_Velocity_Probe(e:f));
    e = e + 500;
    f = f + 500;
end

Velocity_Probe = Velocity_Probe';
f6 = figure;
plot(Velocity_Probe, y_axis);
title('Port 6')
xlabel('Velocity [m/s]')
ylabel('Height [mm]')

% Port 7

Port_7A = readtable('BoundaryLayer_S303_5.csv');
Port_7B = readtable('BoundaryLayer_S303_6.csv');
Port_7A = table2array(Port_7A);
Port_7B = table2array(Port_7B);
a = 1;
b = 500;
c = 1;
d = 1000;
for i = 1:12
    Port_7(c:d,:) = [Port_7B(a:b,:); Port_7A(a:b,:)];
    a = a + 500;
    b = b + 500;
    c = c + 1000;
    d = d + 1000;
end

Port_7_P_atm = mean(Port_7(:,1)); % P_atm average [Pa]
Port_7_T_atm = mean(Port_7(:,2)); % T_atm average [K]
R = 287;
rho = (Port_7_P_atm/(R*Port_7_T_atm)); % density [kg/m^3]
Port_7_Aux_Diff_Pressure = Port_7(:,4); % [Pa]
Port_7_Velocity_Probe = sqrt(2.*Port_7_Aux_Diff_Pressure/rho); % [m/s]
Port_7_Free_Stream_Pressure = mean(Port_7(11001:12000,3)); % [Pa]
Port_7_Free_Stream_Velocity = sqrt(2*Port_7_Free_Stream_Pressure/rho); % [m/s]
Port_7_Boundary_Layer_Velocity = 0.95*Port_7_Free_Stream_Velocity; % [m/s]
Port_7_y_axis = abs(Port_7(:,6)); % vertical height[mm]

x = (Port_7_Velocity_Probe > Port_7_Boundary_Layer_Velocity);
Port_7_Boundary_Layer_Height = Port_7_y_axis(x);
Boundary_Layer_Height(7) = Port_7_Boundary_Layer_Height(1); % [mm]

e = 1;
f = 500;
for j = 1:22
    y_axis(j) = mean(Port_7_y_axis(e:f));
    e = e + 500;
    f = f + 500;
end

y_axis = y_axis';

e = 1;
f = 500;
for k = 1:22
    Velocity_Probe(k) = mean(Port_7_Velocity_Probe(e:f));
    e = e + 500;
    f = f + 500;
end

Velocity_Probe = Velocity_Probe';
f7 = figure;
plot(Velocity_Probe, y_axis);
title('Port 7')
xlabel('Velocity [m/s]')
ylabel('Height [mm]')

% Port 8

Port_8A = readtable('BoundaryLayer_S302_5.csv');
Port_8B = readtable('BoundaryLayer_S302_6.csv');
Port_8A = table2array(Port_8A);
Port_8B = table2array(Port_8B);
a = 1;
b = 500;
c = 1;
d = 1000;
for i = 1:12
    Port_8(c:d,:) = [Port_8B(a:b,:); Port_8A(a:b,:)];
    a = a + 500;
    b = b + 500;
    c = c + 1000;
    d = d + 1000;
end

Port_8_P_atm = mean(Port_8(:,1)); % P_atm average [Pa]
Port_8_T_atm = mean(Port_8(:,2)); % T_atm average [K]
R = 287;
rho = (Port_8_P_atm/(R*Port_8_T_atm)); % density [kg/m^3]
Port_8_Aux_Diff_Pressure = Port_8(:,4); % [Pa]
Port_8_Velocity_Probe = sqrt(2.*Port_8_Aux_Diff_Pressure/rho); % [m/s]
Port_8_Free_Stream_Pressure = mean(Port_8(11001:12000,3)); % [Pa]
Port_8_Free_Stream_Velocity = sqrt(2*Port_8_Free_Stream_Pressure/rho); % [m/s]
Port_8_Boundary_Layer_Velocity = 0.95*Port_8_Free_Stream_Velocity; % [m/s]
Port_8_y_axis = abs(Port_8(:,6)); % vertical height[mm]

x = (Port_8_Velocity_Probe > Port_8_Boundary_Layer_Velocity);
Port_8_Boundary_Layer_Height = Port_8_y_axis(x);
Boundary_Layer_Height(8) = Port_8_Boundary_Layer_Height(1); % [mm]

e = 1;
f = 500;
for j = 1:22
    y_axis(j) = mean(Port_8_y_axis(e:f));
    e = e + 500;
    f = f + 500;
end

y_axis = y_axis';

e = 1;
f = 500;
for k = 1:22
    Velocity_Probe(k) = mean(Port_8_Velocity_Probe(e:f));
    e = e + 500;
    f = f + 500;
end

Velocity_Probe = Velocity_Probe';
f8 = figure;
plot(Velocity_Probe, y_axis);
title('Port 8')
xlabel('Velocity [m/s]')
ylabel('Height [mm]')

% Port 9

Port_9A = readtable('BoundaryLayer_S303_7.csv');
Port_9B = readtable('BoundaryLayer_S303_8.csv');
Port_9A = table2array(Port_9A);
Port_9B = table2array(Port_9B);
a = 1;
b = 500;
c = 1;
d = 1000;
for i = 1:12
    Port_9(c:d,:) = [Port_9B(a:b,:); Port_9A(a:b,:)];
    a = a + 500;
    b = b + 500;
    c = c + 1000;
    d = d + 1000;
end

Port_9_P_atm = mean(Port_9(:,1)); % P_atm average [Pa]
Port_9_T_atm = mean(Port_9(:,2)); % T_atm average [K]
R = 287;
rho = (Port_9_P_atm/(R*Port_9_T_atm)); % density [kg/m^3]
Port_9_Aux_Diff_Pressure = Port_9(:,4); % [Pa]
Port_9_Velocity_Probe = sqrt(2.*Port_9_Aux_Diff_Pressure/rho); % [m/s]
Port_9_Free_Stream_Pressure = mean(Port_9(11001:12000,3)); % [Pa]
Port_9_Free_Stream_Velocity = sqrt(2*Port_9_Free_Stream_Pressure/rho); % [m/s]
Port_9_Boundary_Layer_Velocity = 0.95*Port_9_Free_Stream_Velocity; % [m/s]
Port_9_y_axis = abs(Port_9(:,6)); % vertical height[mm]

x = (Port_9_Velocity_Probe > Port_9_Boundary_Layer_Velocity);
Port_9_Boundary_Layer_Height = Port_9_y_axis(x);
Boundary_Layer_Height(9) = Port_9_Boundary_Layer_Height(1); % [mm]

e = 1;
f = 500;
for j = 1:22
    y_axis(j) = mean(Port_9_y_axis(e:f));
    e = e + 500;
    f = f + 500;
end

y_axis = y_axis';

e = 1;
f = 500;
for k = 1:22
    Velocity_Probe(k) = mean(Port_9_Velocity_Probe(e:f));
    e = e + 500;
    f = f + 500;
end

Velocity_Probe = Velocity_Probe';
f9 = figure;
plot(Velocity_Probe, y_axis);
title('Port 9')
xlabel('Velocity [m/s]')
ylabel('Height [mm]')

% Port 10

Port_10A = readtable('BoundaryLayer_S301_7.csv');
Port_10B = readtable('BoundaryLayer_S301_8.csv');
Port_10A = table2array(Port_10A);
Port_10B = table2array(Port_10B);
a = 1;
b = 500;
c = 1;
d = 1000;
for i = 1:12
    Port_10(c:d,:) = [Port_10B(a:b,:); Port_10A(a:b,:)];
    a = a + 500;
    b = b + 500;
    c = c + 1000;
    d = d + 1000;
end

Port_10_P_atm = mean(Port_10(:,1)); % P_atm average [Pa]
Port_10_T_atm = mean(Port_10(:,2)); % T_atm average [K]
R = 287;
rho = (Port_10_P_atm/(R*Port_10_T_atm)); % density [kg/m^3]
Port_10_Aux_Diff_Pressure = Port_10(:,4); % [Pa]
Port_10_Velocity_Probe = sqrt(2.*Port_10_Aux_Diff_Pressure/rho); % [m/s]
Port_10_Free_Stream_Pressure = mean(Port_10(11001:12000,3)); % [Pa]
Port_10_Free_Stream_Velocity = sqrt(2*Port_10_Free_Stream_Pressure/rho); % [m/s]
Port_10_Boundary_Layer_Velocity = 0.95*Port_10_Free_Stream_Velocity; % [m/s]
Port_10_y_axis = abs(Port_10(:,6)); % vertical height[mm]

x = (Port_10_Velocity_Probe > Port_10_Boundary_Layer_Velocity);
Port_10_Boundary_Layer_Height = Port_10_y_axis(x);
Boundary_Layer_Height(10) = Port_10_Boundary_Layer_Height(1); % [mm]

e = 1;
f = 500;
for j = 1:22
    y_axis(j) = mean(Port_10_y_axis(e:f));
    e = e + 500;
    f = f + 500;
end

y_axis = y_axis';

e = 1;
f = 500;
for k = 1:22
    Velocity_Probe(k) = mean(Port_10_Velocity_Probe(e:f));
    e = e + 500;
    f = f + 500;
end

Velocity_Probe = Velocity_Probe';
f10 = figure;
plot(Velocity_Probe, y_axis);
title('Port 10')
xlabel('Velocity [m/s]')
ylabel('Height [mm]')

% Port 11

Port_11A = readtable('BoundaryLayer_S302_7.csv');
Port_11B = readtable('BoundaryLayer_S302_8.csv');
Port_11A = table2array(Port_11A);
Port_11B = table2array(Port_11B);
a = 1;
b = 500;
c = 1;
d = 1000;
for i = 1:12
    Port_11(c:d,:) = [Port_11B(a:b,:); Port_11A(a:b,:)];
    a = a + 500;
    b = b + 500;
    c = c + 1000;
    d = d + 1000;
end

Port_11_P_atm = mean(Port_11(:,1)); % P_atm average [Pa]
Port_11_T_atm = mean(Port_11(:,2)); % T_atm average [K]
R = 287;
rho = (Port_11_P_atm/(R*Port_11_T_atm)); % density [kg/m^3]
Port_11_Aux_Diff_Pressure = Port_11(:,4); % [Pa]
Port_11_Velocity_Probe = sqrt(2.*Port_11_Aux_Diff_Pressure/rho); % [m/s]
Port_11_Free_Stream_Pressure = mean(Port_11(11001:12000,3)); % [Pa]
Port_11_Free_Stream_Velocity = sqrt(2*Port_11_Free_Stream_Pressure/rho); % [m/s]
Port_11_Boundary_Layer_Velocity = 0.95*Port_11_Free_Stream_Velocity; % [m/s]
Port_11_y_axis = abs(Port_11(:,6)); % vertical height[mm]

x = (Port_11_Velocity_Probe > Port_11_Boundary_Layer_Velocity);
Port_11_Boundary_Layer_Height = Port_11_y_axis(x);
Boundary_Layer_Height(11) = Port_11_Boundary_Layer_Height(1); % [mm]

e = 1;
f = 500;
for j = 1:22
    y_axis(j) = mean(Port_11_y_axis(e:f));
    e = e + 500;
    f = f + 500;
end

y_axis = y_axis';

e = 1;
f = 500;
for k = 1:22
    Velocity_Probe(k) = mean(Port_11_Velocity_Probe(e:f));
    e = e + 500;
    f = f + 500;
end

Velocity_Probe = Velocity_Probe';
f11 = figure;
plot(Velocity_Probe, y_axis);
title('Port 11')
xlabel('Velocity [m/s]')
ylabel('Height [mm]')

f12 = figure;
plot(Boundary_Layer_Height)