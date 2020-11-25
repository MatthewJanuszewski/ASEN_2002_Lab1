% ASEN 2002 Lab 2 Wind Tunnel
% Section 3: Cambered Airfoil Aerodynamics

clc
clear all
close all

% Read Data

Data1 = readtable('AirfoilPressure_S301_4');
Data1 = table2array(Data1);

% Solve for Pressure at Port 11:
a = 1;
b = 20;

for i = 1:12
    Pressure_Port_1(i) = sum(Data1(a:b,7))/20;
    Pressure_Port_2(i) = sum(Data1(a:b,8))/20;
    Pressure_Port_3(i) = sum(Data1(a:b,9))/20;
    Pressure_Port_4(i) = sum(Data1(a:b,10))/20;
    Pressure_Port_5(i) = sum(Data1(a:b,11))/20;
    Pressure_Port_6(i) = sum(Data1(a:b,12))/20;
    Pressure_Port_7(i) = sum(Data1(a:b,13))/20;
    Pressure_Port_8(i) = sum(Data1(a:b,14))/20;
    Pressure_Port_10(i) = sum(Data1(a:b,15))/20;
    
    Top = Pressure_Port_10(i) - Pressure_Port_8(i);
    Pressure_Port_11_Top = Top + Pressure_Port_10(i);

    Pressure_Port_12(i) = sum(Data1(a:b,16))/20;
    Pressure_Port_14(i) = sum(Data1(a:b,17))/20;
    
    Bottom = Pressure_Port_12(i) - Pressure_Port_14(i);
    Pressure_Port_11_Bottom = Bottom + Pressure_Port_12(i);

    Pressure_Port_11(i) = (Pressure_Port_11_Top + Pressure_Port_11_Bottom)/2;
    Pressure_Port_16(i) = sum(Data1(a:b,18))/20;
    Pressure_Port_17(i) = sum(Data1(a:b,19))/20;
    Pressure_Port_18(i) = sum(Data1(a:b,20))/20;
    Pressure_Port_19(i) = sum(Data1(a:b,21))/20;
    Pressure_Port_20(i) = sum(Data1(a:b,22))/20;
    V_infinity(i) = mean(Data1(a:b,4));
    
    a = a + 20;
    b = b + 20;
end

rho(1) = mean(Data1(:,3));

Data2 = readtable('AirfoilPressure_S301_3.csv');
Data2 = table2array(Data2);

a = 1;
b = 20;

for i = 13:24
    Pressure_Port_1(i) = sum(Data2(a:b,7))/20;
    Pressure_Port_2(i) = sum(Data2(a:b,8))/20;
    Pressure_Port_3(i) = sum(Data2(a:b,9))/20;
    Pressure_Port_4(i) = sum(Data2(a:b,10))/20;
    Pressure_Port_5(i) = sum(Data2(a:b,11))/20;
    Pressure_Port_6(i) = sum(Data2(a:b,12))/20;
    Pressure_Port_7(i) = sum(Data2(a:b,13))/20;
    
    Pressure_Port_8(i) = sum(Data2(a:b,14))/20;
    Pressure_Port_10(i) = sum(Data2(a:b,15))/20;
    
    Top = Pressure_Port_10(i) - Pressure_Port_8(i);
    Pressure_Port_11_Top = Top + Pressure_Port_10(i);

    Pressure_Port_12(i) = sum(Data2(a:b,16))/20;
    Pressure_Port_14(i) = sum(Data2(a:b,17))/20;
    
    Bottom = Pressure_Port_12(i) - Pressure_Port_14(i);
    Pressure_Port_11_Bottom = Bottom + Pressure_Port_12(i);

    Pressure_Port_11(i) = (Pressure_Port_11_Top + Pressure_Port_11_Bottom)/2;
    Pressure_Port_16(i) = sum(Data2(a:b,18))/20;
    Pressure_Port_17(i) = sum(Data2(a:b,19))/20;
    Pressure_Port_18(i) = sum(Data2(a:b,20))/20;
    Pressure_Port_19(i) = sum(Data2(a:b,21))/20;
    Pressure_Port_20(i) = sum(Data2(a:b,22))/20;
    
    V_infinity(i) = mean(Data2(a:b,4));
    
    a = a + 20;
    b = b + 20;
end

rho(2) = mean(Data2(:,3));

Data3 = readtable('AirfoilPressure_S301_2.csv');
Data3 = table2array(Data3);

a = 1;
b = 20;

for i = 25:36
    Pressure_Port_1(i) = sum(Data3(a:b,7))/20;
    Pressure_Port_2(i) = sum(Data3(a:b,8))/20;
    Pressure_Port_3(i) = sum(Data3(a:b,9))/20;
    Pressure_Port_4(i) = sum(Data3(a:b,10))/20;
    Pressure_Port_5(i) = sum(Data3(a:b,11))/20;
    Pressure_Port_6(i) = sum(Data3(a:b,12))/20;
    Pressure_Port_7(i) = sum(Data3(a:b,13))/20;
    
    Pressure_Port_8(i) = sum(Data3(a:b,14))/20;
    Pressure_Port_10(i) = sum(Data3(a:b,15))/20;
    
    Top = Pressure_Port_10(i) - Pressure_Port_8(i);
    Pressure_Port_11_Top = Top + Pressure_Port_10(i);

    Pressure_Port_12(i) = sum(Data3(a:b,16))/20;
    Pressure_Port_14(i) = sum(Data3(a:b,17))/20;
    
    Bottom = Pressure_Port_12(i) - Pressure_Port_14(i);
    Pressure_Port_11_Bottom = Bottom + Pressure_Port_12(i);

    Pressure_Port_11(i) = (Pressure_Port_11_Top + Pressure_Port_11_Bottom)/2;
    Pressure_Port_16(i) = sum(Data3(a:b,18))/20;
    Pressure_Port_17(i) = sum(Data3(a:b,19))/20;
    Pressure_Port_18(i) = sum(Data3(a:b,20))/20;
    Pressure_Port_19(i) = sum(Data3(a:b,21))/20;
    Pressure_Port_20(i) = sum(Data3(a:b,22))/20;
    
    V_infinity(i) = mean(Data3(a:b,4));
    
    a = a + 20;
    b = b + 20;
end

rho(3) = mean(Data3(:,3));

Data4 = readtable('AirfoilPressure_S301_1.csv');
Data4 = table2array(Data4);

a = 1;
b = 20;

for i = 37:48
    Pressure_Port_1(i) = sum(Data4(a:b,7))/20;
    Pressure_Port_2(i) = sum(Data4(a:b,8))/20;
    Pressure_Port_3(i) = sum(Data4(a:b,9))/20;
    Pressure_Port_4(i) = sum(Data4(a:b,10))/20;
    Pressure_Port_5(i) = sum(Data4(a:b,11))/20;
    Pressure_Port_6(i) = sum(Data4(a:b,12))/20;
    Pressure_Port_7(i) = sum(Data4(a:b,13))/20;
    
    Pressure_Port_8(i) = sum(Data4(a:b,14))/20;
    Pressure_Port_10(i) = sum(Data4(a:b,15))/20;
    
    Top = Pressure_Port_10(i) - Pressure_Port_8(i);
    Pressure_Port_11_Top = Top + Pressure_Port_10(i);

    Pressure_Port_12(i) = sum(Data4(a:b,16))/20;
    Pressure_Port_14(i) = sum(Data4(a:b,17))/20;
    
    Bottom = Pressure_Port_12(i) - Pressure_Port_14(i);
    Pressure_Port_11_Bottom = Bottom + Pressure_Port_12(i);

    Pressure_Port_11(i) = (Pressure_Port_11_Top + Pressure_Port_11_Bottom)/2;
    Pressure_Port_16(i) = sum(Data4(a:b,18))/20;
    Pressure_Port_17(i) = sum(Data4(a:b,19))/20;
    Pressure_Port_18(i) = sum(Data4(a:b,20))/20;
    Pressure_Port_19(i) = sum(Data4(a:b,21))/20;
    Pressure_Port_20(i) = sum(Data4(a:b,22))/20;
    
    V_infinity(i) = mean(Data4(a:b,4));
    
    a = a + 20;
    b = b + 20;
end

rho(4) = mean(Data4(:,3));

Data5 = readtable('AirfoilPressure_S301_8.csv');
Data5 = table2array(Data5);

a = 1;
b = 20;

for i = 49:60
    Pressure_Port_1(i) = sum(Data5(a:b,7))/20;
    Pressure_Port_2(i) = sum(Data5(a:b,8))/20;
    Pressure_Port_3(i) = sum(Data5(a:b,9))/20;
    Pressure_Port_4(i) = sum(Data5(a:b,10))/20;
    Pressure_Port_5(i) = sum(Data5(a:b,11))/20;
    Pressure_Port_6(i) = sum(Data5(a:b,12))/20;
    Pressure_Port_7(i) = sum(Data5(a:b,13))/20;
    
    Pressure_Port_8(i) = sum(Data5(a:b,14))/20;
    Pressure_Port_10(i) = sum(Data5(a:b,15))/20;
    
    Top = Pressure_Port_10(i) - Pressure_Port_8(i);
    Pressure_Port_11_Top = Top + Pressure_Port_10(i);

    Pressure_Port_12(i) = sum(Data5(a:b,16))/20;
    Pressure_Port_14(i) = sum(Data5(a:b,17))/20;
    
    Bottom = Pressure_Port_12(i) - Pressure_Port_14(i);
    Pressure_Port_11_Bottom = Bottom + Pressure_Port_12(i);

    Pressure_Port_11(i) = (Pressure_Port_11_Top + Pressure_Port_11_Bottom)/2;
    Pressure_Port_16(i) = sum(Data5(a:b,18))/20;
    Pressure_Port_17(i) = sum(Data5(a:b,19))/20;
    Pressure_Port_18(i) = sum(Data5(a:b,20))/20;
    Pressure_Port_19(i) = sum(Data5(a:b,21))/20;
    Pressure_Port_20(i) = sum(Data5(a:b,22))/20;
    
    V_infinity(i) = mean(Data5(a:b,4));
    
    a = a + 20;
    b = b + 20;
end

rho(5) = mean(Data5(:,3));

Data6 = readtable('AirfoilPressure_S301_7.csv');
Data6 = table2array(Data6);

a = 1;
b = 20;

for i = 61:72
    Pressure_Port_1(i) = sum(Data6(a:b,7))/20;
    Pressure_Port_2(i) = sum(Data6(a:b,8))/20;
    Pressure_Port_3(i) = sum(Data6(a:b,9))/20;
    Pressure_Port_4(i) = sum(Data6(a:b,10))/20;
    Pressure_Port_5(i) = sum(Data6(a:b,11))/20;
    Pressure_Port_6(i) = sum(Data6(a:b,12))/20;
    Pressure_Port_7(i) = sum(Data6(a:b,13))/20;
    
    Pressure_Port_8(i) = sum(Data6(a:b,14))/20;
    Pressure_Port_10(i) = sum(Data6(a:b,15))/20;
    
    Top = Pressure_Port_10(i) - Pressure_Port_8(i);
    Pressure_Port_11_Top = Top + Pressure_Port_10(i);

    Pressure_Port_12(i) = sum(Data6(a:b,16))/20;
    Pressure_Port_14(i) = sum(Data6(a:b,17))/20;
    
    Bottom = Pressure_Port_12(i) - Pressure_Port_14(i);
    Pressure_Port_11_Bottom = Bottom + Pressure_Port_12(i);

    Pressure_Port_11(i) = (Pressure_Port_11_Top + Pressure_Port_11_Bottom)/2;
    Pressure_Port_16(i) = sum(Data6(a:b,18))/20;
    Pressure_Port_17(i) = sum(Data6(a:b,19))/20;
    Pressure_Port_18(i) = sum(Data6(a:b,20))/20;
    Pressure_Port_19(i) = sum(Data6(a:b,21))/20;
    Pressure_Port_20(i) = sum(Data6(a:b,22))/20;
    
    V_infinity(i) = mean(Data6(a:b,4));
    
    a = a + 20;
    b = b + 20;
end

rho(6) = mean(Data6(:,3));

Data7 = readtable('AirfoilPressure_S301_6.csv');
Data7 = table2array(Data7);

a = 1;
b = 20;

for i = 73:84
    Pressure_Port_1(i) = sum(Data7(a:b,7))/20;
    Pressure_Port_2(i) = sum(Data7(a:b,8))/20;
    Pressure_Port_3(i) = sum(Data7(a:b,9))/20;
    Pressure_Port_4(i) = sum(Data7(a:b,10))/20;
    Pressure_Port_5(i) = sum(Data7(a:b,11))/20;
    Pressure_Port_6(i) = sum(Data7(a:b,12))/20;
    Pressure_Port_7(i) = sum(Data7(a:b,13))/20;
    
    Pressure_Port_8(i) = sum(Data7(a:b,14))/20;
    Pressure_Port_10(i) = sum(Data7(a:b,15))/20;
    
    Top = Pressure_Port_10(i) - Pressure_Port_8(i);
    Pressure_Port_11_Top = Top + Pressure_Port_10(i);

    Pressure_Port_12(i) = sum(Data7(a:b,16))/20;
    Pressure_Port_14(i) = sum(Data7(a:b,17))/20;
    
    Bottom = Pressure_Port_12(i) - Pressure_Port_14(i);
    Pressure_Port_11_Bottom = Bottom + Pressure_Port_12(i);

    Pressure_Port_11(i) = (Pressure_Port_11_Top + Pressure_Port_11_Bottom)/2;
    Pressure_Port_16(i) = sum(Data7(a:b,18))/20;
    Pressure_Port_17(i) = sum(Data7(a:b,19))/20;
    Pressure_Port_18(i) = sum(Data7(a:b,20))/20;
    Pressure_Port_19(i) = sum(Data7(a:b,21))/20;
    Pressure_Port_20(i) = sum(Data7(a:b,22))/20;
    
    V_infinity(i) = mean(Data7(a:b,4));
    
    a = a + 20;
    b = b + 20;
end

rho(7) = mean(Data7(:,3));

Data8 = readtable('AirfoilPressure_S301_5.csv');
Data8 = table2array(Data8);

a = 1;
b = 20;

for i = 85:96
    Pressure_Port_1(i) = sum(Data8(a:b,7))/20;
    Pressure_Port_2(i) = sum(Data8(a:b,8))/20;
    Pressure_Port_3(i) = sum(Data8(a:b,9))/20;
    Pressure_Port_4(i) = sum(Data8(a:b,10))/20;
    Pressure_Port_5(i) = sum(Data8(a:b,11))/20;
    Pressure_Port_6(i) = sum(Data8(a:b,12))/20;
    Pressure_Port_7(i) = sum(Data8(a:b,13))/20;
    
    Pressure_Port_8(i) = sum(Data8(a:b,14))/20;
    Pressure_Port_10(i) = sum(Data8(a:b,15))/20;
    Top = Pressure_Port_10(i) - Pressure_Port_8(i);
    Pressure_Port_11_Top = Top + Pressure_Port_10(i);

    Pressure_Port_12(i) = sum(Data8(a:b,16))/20;
    Pressure_Port_14(i) = sum(Data8(a:b,17))/20;
    Bottom = Pressure_Port_12(i) - Pressure_Port_14(i);
    Pressure_Port_11_Bottom = Bottom + Pressure_Port_12(i);

    Pressure_Port_11(i) = (Pressure_Port_11_Top + Pressure_Port_11_Bottom)/2;
    Pressure_Port_16(i) = sum(Data8(a:b,18))/20;
    Pressure_Port_17(i) = sum(Data8(a:b,19))/20;
    Pressure_Port_18(i) = sum(Data8(a:b,20))/20;
    Pressure_Port_19(i) = sum(Data8(a:b,21))/20;
    Pressure_Port_20(i) = sum(Data8(a:b,22))/20;
    
    V_infinity(i) = mean(Data8(a:b,4));
    
    a = a + 20;
    b = b + 20;
end

rho(8) = mean(Data8(:,3));


% Find Cp for angle of attack = -15 and Velocity = 9.1

q_infinity(32) = 0.5*rho(3)*V_infinity(32);

Cp = zeros(17,2);

Cp(1,1) = 0;
Cp(2,1) = 0.175*0.0254;
Cp(3,1) = 0.35*0.0254;
Cp(4,1) = 0.7*0.0254;
Cp(5,1) = 1.05*0.0254;
Cp(6,1) = 1.4*0.0254;
Cp(7,1) = 1.75*0.0254;
Cp(8,1) = 2.1*0.0254;
Cp(9,1) = 2.8*0.0254;
Cp(10,1) = 3.5*0.0254;
Cp(11,1) = 2.8*0.0254;
Cp(12,1) = 2.1*0.0254;
Cp(13,1) = 1.4*0.0254;
Cp(14,1) = 1.05*0.0254;
Cp(15,1) = 0.7*0.0254;
Cp(16,1) = 0.35*0.0254;
Cp(17,1) = 0.175*0.0254;

Cp(1,3) = 0.14665*0.0254;
Cp(2,3) = 0.33075*0.0254;
Cp(3,3) = 0.4018*0.0254;
Cp(4,3) = 0.476*0.0254;
Cp(5,3) = 0.49*0.0254;
Cp(6,3) = 0.4774*0.0254;
Cp(7,3) = 0.4403*0.0254;
Cp(8,3) = 0.38325*0.0254;
Cp(9,3) = 0.21875*0.0254;
Cp(10,3) = 0*0.0254;
Cp(11,3) = 0*0.0254;
Cp(12,3) = 0*0.0254;
Cp(13,3) = 0*0.0254;
Cp(14,3) = 0*0.0254;
Cp(15,3) = 0.0014*0.0254;
Cp(16,3) = 0.0175*0.0254;
Cp(17,3) = 0.03885*0.0254;
    
Cp(1,2) = Pressure_Port_1(32)/(q_infinity(32));
Cp(2,2) = Pressure_Port_2(32)/(q_infinity(32));
Cp(3,2) = Pressure_Port_3(32)/(q_infinity(32));
Cp(4,2) = Pressure_Port_4(32)/(q_infinity(32));
Cp(5,2) = Pressure_Port_5(32)/(q_infinity(32));
Cp(6,2) = Pressure_Port_6(32)/(q_infinity(32));
Cp(7,2) = Pressure_Port_7(32)/(q_infinity(32));
Cp(8,2) = Pressure_Port_8(32)/(q_infinity(32));
Cp(9,2) = Pressure_Port_10(32)/(q_infinity(32));
Cp(10,2) = Pressure_Port_11(32)/(q_infinity(32));
Cp(11,2) = Pressure_Port_12(32)/(q_infinity(32));
Cp(12,2) = Pressure_Port_14(32)/(q_infinity(32));
Cp(13,2) = Pressure_Port_16(32)/(q_infinity(32));
Cp(14,2) = Pressure_Port_17(32)/(q_infinity(32));
Cp(15,2) = Pressure_Port_18(32)/(q_infinity(32));
Cp(16,2) = Pressure_Port_19(32)/(q_infinity(32));
Cp(17,2) = Pressure_Port_20(32)/(q_infinity(32));

Cp_norm = normalize(Cp);

f1 = figure(1);
plot(Cp_norm(:,1),Cp_norm(:,2));
f1.CurrentAxes.YDir = 'Reverse';

c = 3.50307097;

for i = 2:17
    Cn(i) = 0.5*(Cp(i-1,2) + Cp(i,2))*(Cp(i,1) - Cp(i-1,1))/c;
    Ca(i) = 0.5*(Cp(i-1,2) + Cp(i,2))*(Cp(i,3) - Cp(i-1,3))/c;
end

Norm_Coef = sum(Cn)
Axial_Coef = sum(Ca)

Lift_Coef = (Norm_Coef*cos(3))-(Axial_Coef*sin(3))
Drag_Coef = (Norm_Coef*sin(3))+(Axial_Coef*cos(3))