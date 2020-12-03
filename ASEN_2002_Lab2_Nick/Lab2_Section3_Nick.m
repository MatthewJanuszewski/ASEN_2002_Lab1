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
    Pressure_Port_1(i) = mean(Data1(a:b,7));
    Pressure_Port_2(i) = mean(Data1(a:b,8));
    Pressure_Port_3(i) = mean(Data1(a:b,9));
    Pressure_Port_4(i) = mean(Data1(a:b,10));
    Pressure_Port_5(i) = mean(Data1(a:b,11));
    Pressure_Port_6(i) = mean(Data1(a:b,12));
    Pressure_Port_7(i) = mean(Data1(a:b,13));
    Pressure_Port_8(i) = mean(Data1(a:b,14));
    Pressure_Port_10(i) = mean(Data1(a:b,15));
    
    Top = Pressure_Port_10(i) - Pressure_Port_8(i);
    Pressure_Port_11_Top = Top + Pressure_Port_10(i);

    Pressure_Port_12(i) = mean(Data1(a:b,16));
    Pressure_Port_14(i) = mean(Data1(a:b,17));
    
    Bottom = Pressure_Port_12(i) - Pressure_Port_14(i);
    Pressure_Port_11_Bottom = Bottom + Pressure_Port_12(i);

    Pressure_Port_11(i) = (Pressure_Port_11_Top + Pressure_Port_11_Bottom)/2;
    Pressure_Port_16(i) = mean(Data1(a:b,18));
    Pressure_Port_17(i) = mean(Data1(a:b,19));
    Pressure_Port_18(i) = mean(Data1(a:b,20));
    Pressure_Port_19(i) = mean(Data1(a:b,21));
    Pressure_Port_20(i) = mean(Data1(a:b,22));
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
    Pressure_Port_1(i) = mean(Data2(a:b,7));
    Pressure_Port_2(i) = mean(Data2(a:b,8));
    Pressure_Port_3(i) = mean(Data2(a:b,9));
    Pressure_Port_4(i) = mean(Data2(a:b,10));
    Pressure_Port_5(i) = mean(Data2(a:b,11));
    Pressure_Port_6(i) = mean(Data2(a:b,12));
    Pressure_Port_7(i) = mean(Data2(a:b,13));
    
    Pressure_Port_8(i) = mean(Data2(a:b,14));
    Pressure_Port_10(i) = mean(Data2(a:b,15));
    
    Top(i) = Pressure_Port_10(i) - Pressure_Port_8(i);
    Pressure_Port_11_Top = Top(i) + Pressure_Port_10(i);

    Pressure_Port_12(i) = mean(Data2(a:b,16));
    Pressure_Port_14(i) = mean(Data2(a:b,17));
    
    Bottom(i) = Pressure_Port_12(i) - Pressure_Port_14(i);
    Pressure_Port_11_Bottom = Bottom(i) + Pressure_Port_12(i);

    Pressure_Port_11(i) = (Pressure_Port_11_Top + Pressure_Port_11_Bottom)/2;
    Pressure_Port_16(i) = mean(Data2(a:b,18));
    Pressure_Port_17(i) = mean(Data2(a:b,19));
    Pressure_Port_18(i) = mean(Data2(a:b,20));
    Pressure_Port_19(i) = mean(Data2(a:b,21));
    Pressure_Port_20(i) = mean(Data2(a:b,22));
    
    V_infinity(i) = mean(Data2(a:b,4));
    AOA(i) = Data2(i,23);
    
    a = a + 20;
    b = b + 20;
end

rho(2) = mean(Data2(:,3));

Data3 = readtable('AirfoilPressure_S301_2.csv');
Data3 = table2array(Data3);

a = 1;
b = 20;

for i = 25:36
    Pressure_Port_1(i) = mean(Data3(a:b,7));
    Pressure_Port_2(i) = mean(Data3(a:b,8));
    Pressure_Port_3(i) = mean(Data3(a:b,9));
    Pressure_Port_4(i) = mean(Data3(a:b,10));
    Pressure_Port_5(i) = mean(Data3(a:b,11));
    Pressure_Port_6(i) = mean(Data3(a:b,12));
    Pressure_Port_7(i) = mean(Data3(a:b,13));
    
    Pressure_Port_8(i) = mean(Data3(a:b,14));
    Pressure_Port_10(i) = mean(Data3(a:b,15));
    
    Top = Pressure_Port_10(i) - Pressure_Port_8(i);
    Pressure_Port_11_Top = Top + Pressure_Port_10(i);

    Pressure_Port_12(i) = mean(Data3(a:b,16));
    Pressure_Port_14(i) = mean(Data3(a:b,17));
    
    Bottom = Pressure_Port_12(i) - Pressure_Port_14(i);
    Pressure_Port_11_Bottom = Bottom + Pressure_Port_12(i);

    Pressure_Port_11(i) = (Pressure_Port_11_Top + Pressure_Port_11_Bottom)/2;
    Pressure_Port_16(i) = mean(Data3(a:b,18));
    Pressure_Port_17(i) = mean(Data3(a:b,19));
    Pressure_Port_18(i) = mean(Data3(a:b,20));
    Pressure_Port_19(i) = mean(Data3(a:b,21));
    Pressure_Port_20(i) = mean(Data3(a:b,22));
    
    V_infinity(i) = mean(Data3(a:b,4));
    AOA(i) = Data3(i,23);
    
    a = a + 20;
    b = b + 20;
end

rho(3) = mean(Data3(:,3));

Data4 = readtable('AirfoilPressure_S301_1.csv');
Data4 = table2array(Data4);

a = 1;
b = 20;

for i = 37:48
    Pressure_Port_1(i) = mean(Data4(a:b,7));
    Pressure_Port_2(i) = mean(Data4(a:b,8));
    Pressure_Port_3(i) = mean(Data4(a:b,9));
    Pressure_Port_4(i) = mean(Data4(a:b,10));
    Pressure_Port_5(i) = mean(Data4(a:b,11));
    Pressure_Port_6(i) = mean(Data4(a:b,12));
    Pressure_Port_7(i) = mean(Data4(a:b,13));
    
    Pressure_Port_8(i) = mean(Data4(a:b,14));
    Pressure_Port_10(i) = mean(Data4(a:b,15));
    
    Top = Pressure_Port_10(i) - Pressure_Port_8(i);
    Pressure_Port_11_Top = Top + Pressure_Port_10(i);

    Pressure_Port_12(i) = mean(Data4(a:b,16));
    Pressure_Port_14(i) = mean(Data4(a:b,17));
    
    Bottom = Pressure_Port_12(i) - Pressure_Port_14(i);
    Pressure_Port_11_Bottom = Bottom + Pressure_Port_12(i);

    Pressure_Port_11(i) = (Pressure_Port_11_Top + Pressure_Port_11_Bottom)/2;
    Pressure_Port_16(i) = mean(Data4(a:b,18));
    Pressure_Port_17(i) = mean(Data4(a:b,19));
    Pressure_Port_18(i) = mean(Data4(a:b,20));
    Pressure_Port_19(i) = mean(Data4(a:b,21));
    Pressure_Port_20(i) = mean(Data4(a:b,22));
    
    V_infinity(i) = mean(Data4(a:b,4));
    AOA(i) = Data4(i,23);
    
    a = a + 20;
    b = b + 20;
end

rho(4) = mean(Data4(:,3));

Data5 = readtable('AirfoilPressure_S301_8.csv');
Data5 = table2array(Data5);

a = 1;
b = 20;

for i = 49:60
    Pressure_Port_1(i) = mean(Data5(a:b,7));
    Pressure_Port_2(i) = mean(Data5(a:b,8));
    Pressure_Port_3(i) = mean(Data5(a:b,9));
    Pressure_Port_4(i) = mean(Data5(a:b,10));
    Pressure_Port_5(i) = mean(Data5(a:b,11));
    Pressure_Port_6(i) = mean(Data5(a:b,12));
    Pressure_Port_7(i) = mean(Data5(a:b,13));
    
    Pressure_Port_8(i) = mean(Data5(a:b,14));
    Pressure_Port_10(i) = mean(Data5(a:b,15));
    
    Top = Pressure_Port_10(i) - Pressure_Port_8(i);
    Pressure_Port_11_Top = Top + Pressure_Port_10(i);

    Pressure_Port_12(i) = mean(Data5(a:b,16));
    Pressure_Port_14(i) = mean(Data5(a:b,17));
    
    Bottom = Pressure_Port_12(i) - Pressure_Port_14(i);
    Pressure_Port_11_Bottom = Bottom + Pressure_Port_12(i);

    Pressure_Port_11(i) = (Pressure_Port_11_Top + Pressure_Port_11_Bottom)/2;
    Pressure_Port_16(i) = mean(Data5(a:b,18));
    Pressure_Port_17(i) = mean(Data5(a:b,19));
    Pressure_Port_18(i) = mean(Data5(a:b,20));
    Pressure_Port_19(i) = mean(Data5(a:b,21));
    Pressure_Port_20(i) = mean(Data5(a:b,22));
    
    V_infinity(i) = mean(Data5(a:b,4));
    AOA(i) = Data5(i,23);
    
    a = a + 20;
    b = b + 20;
end

rho(5) = mean(Data5(:,3));

Data6 = readtable('AirfoilPressure_S301_7.csv');
Data6 = table2array(Data6);

a = 1;
b = 20;

for i = 61:72
    Pressure_Port_1(i) = mean(Data6(a:b,7));
    Pressure_Port_2(i) = mean(Data6(a:b,8));
    Pressure_Port_3(i) = mean(Data6(a:b,9));
    Pressure_Port_4(i) = mean(Data6(a:b,10));
    Pressure_Port_5(i) = mean(Data6(a:b,11));
    Pressure_Port_6(i) = mean(Data6(a:b,12));
    Pressure_Port_7(i) = mean(Data6(a:b,13));
    
    Pressure_Port_8(i) = mean(Data6(a:b,14));
    Pressure_Port_10(i) = mean(Data6(a:b,15));
    
    Top = Pressure_Port_10(i) - Pressure_Port_8(i);
    Pressure_Port_11_Top = Top + Pressure_Port_10(i);

    Pressure_Port_12(i) = mean(Data6(a:b,16));
    Pressure_Port_14(i) = mean(Data6(a:b,17));
    
    Bottom = Pressure_Port_12(i) - Pressure_Port_14(i);
    Pressure_Port_11_Bottom = Bottom + Pressure_Port_12(i);

    Pressure_Port_11(i) = (Pressure_Port_11_Top + Pressure_Port_11_Bottom)/2;
    Pressure_Port_16(i) = mean(Data6(a:b,18));
    Pressure_Port_17(i) = mean(Data6(a:b,19));
    Pressure_Port_18(i) = mean(Data6(a:b,20));
    Pressure_Port_19(i) = mean(Data6(a:b,21));
    Pressure_Port_20(i) = mean(Data6(a:b,22));
    
    V_infinity(i) = mean(Data6(a:b,4));
    AOA(i) = Data6(i,23);
    
    a = a + 20;
    b = b + 20;
end

rho(6) = mean(Data6(:,3));

Data7 = readtable('AirfoilPressure_S301_6.csv');
Data7 = table2array(Data7);

a = 1;
b = 20;

for i = 73:84
    Pressure_Port_1(i) = mean(Data7(a:b,7));
    Pressure_Port_2(i) = mean(Data7(a:b,8));
    Pressure_Port_3(i) = mean(Data7(a:b,9));
    Pressure_Port_4(i) = mean(Data7(a:b,10));
    Pressure_Port_5(i) = mean(Data7(a:b,11));
    Pressure_Port_6(i) = mean(Data7(a:b,12));
    Pressure_Port_7(i) = mean(Data7(a:b,13));
    
    Pressure_Port_8(i) = mean(Data7(a:b,14));
    Pressure_Port_10(i) = mean(Data7(a:b,15));
    
    Top = Pressure_Port_10(i) - Pressure_Port_8(i);
    Pressure_Port_11_Top = Top + Pressure_Port_10(i);

    Pressure_Port_12(i) = mean(Data7(a:b,16));
    Pressure_Port_14(i) = mean(Data7(a:b,17));
    
    Bottom = Pressure_Port_12(i) - Pressure_Port_14(i);
    Pressure_Port_11_Bottom = Bottom + Pressure_Port_12(i);

    Pressure_Port_11(i) = (Pressure_Port_11_Top + Pressure_Port_11_Bottom)/2;
    Pressure_Port_16(i) = mean(Data7(a:b,18));
    Pressure_Port_17(i) = mean(Data7(a:b,19));
    Pressure_Port_18(i) = mean(Data7(a:b,20));
    Pressure_Port_19(i) = mean(Data7(a:b,21));
    Pressure_Port_20(i) = mean(Data7(a:b,22));
    
    V_infinity(i) = mean(Data7(a:b,4));
    AOA(i) = Data7(i,23);
    
    a = a + 20;
    b = b + 20;
end

rho(7) = mean(Data7(:,3));

Data8 = readtable('AirfoilPressure_S301_5.csv');
Data8 = table2array(Data8);

a = 1;
b = 20;

for i = 85:96
    Pressure_Port_1(i) = mean(Data8(a:b,7));
    Pressure_Port_2(i) = mean(Data8(a:b,8));
    Pressure_Port_3(i) = mean(Data8(a:b,9));
    Pressure_Port_4(i) = mean(Data8(a:b,10));
    Pressure_Port_5(i) = mean(Data8(a:b,11));
    Pressure_Port_6(i) = mean(Data8(a:b,12));
    Pressure_Port_7(i) = mean(Data8(a:b,13));
    
    Pressure_Port_8(i) = mean(Data8(a:b,14));
    Pressure_Port_10(i) = mean(Data8(a:b,15));
    Top = Pressure_Port_10(i) - Pressure_Port_8(i);
    Pressure_Port_11_Top = Top + Pressure_Port_10(i);

    Pressure_Port_12(i) = mean(Data8(a:b,16));
    Pressure_Port_14(i) = mean(Data8(a:b,17));
    Bottom = Pressure_Port_12(i) - Pressure_Port_14(i);
    Pressure_Port_11_Bottom = Bottom + Pressure_Port_12(i);

    Pressure_Port_11(i) = (Pressure_Port_11_Top + Pressure_Port_11_Bottom)/2;
    Pressure_Port_16(i) = mean(Data8(a:b,18));
    Pressure_Port_17(i) = mean(Data8(a:b,19));
    Pressure_Port_18(i) = mean(Data8(a:b,20));
    Pressure_Port_19(i) = mean(Data8(a:b,21));
    Pressure_Port_20(i) = mean(Data8(a:b,22));
    
    V_infinity(i) = mean(Data8(a:b,4));
    AOA(i) = Data8(i,23);
    
    a = a + 20;
    b = b + 20;
end

rho(8) = mean(Data8(:,3));


% Solving for q_infinity for all 96 configurations
for i = 1:96
    if i <= 12
        Rho = rho(1);
        q_infinity(i) = 0.5*Rho*(V_infinity(i)^2);
    elseif (i > 12) && (i < 25)
        Rho = rho(2);
        q_infinity(i) = 0.5*Rho*(V_infinity(i)^2);
    elseif (i > 24) && (i < 37)
        Rho = rho(3);
        q_infinity(i) = 0.5*Rho*(V_infinity(i)^2);
    elseif (i > 36) && (i < 49)
        Rho = rho(4);
        q_infinity(i) = 0.5*Rho*(V_infinity(i)^2);
    elseif (i > 48) && (i < 61)
        Rho = rho(5);
        q_infinity(i) = 0.5*Rho*(V_infinity(i)^2);
    elseif (i > 60) && (i < 73)
        Rho = rho(6);
        q_infinity(i) = 0.5*Rho*(V_infinity(i)^2);
    elseif (i > 72) && (i < 85)
        Rho = rho(7);
        q_infinity(i) = 0.5*Rho*(V_infinity(i)^2);
    elseif (i > 84) && (i < 97)
        Rho = rho(8);
        q_infinity(i) = 0.5*Rho*(V_infinity(i)^2);
    end
end


% Port location: X values (converted to meters)
Cp_location(1,1) = 0;
Cp_location(2,1) = 0.175*0.0254;
Cp_location(3,1) = 0.35*0.0254;
Cp_location(4,1) = 0.7*0.0254;
Cp_location(5,1) = 1.05*0.0254;
Cp_location(6,1) = 1.4*0.0254;
Cp_location(7,1) = 1.75*0.0254;
Cp_location(8,1) = 2.1*0.0254;
Cp_location(9,1) = 2.8*0.0254;
Cp_location(10,1) = 3.5*0.0254;
Cp_location(11,1) = 2.8*0.0254;
Cp_location(12,1) = 2.1*0.0254;
Cp_location(13,1) = 1.4*0.0254;
Cp_location(14,1) = 1.05*0.0254;
Cp_location(15,1) = 0.7*0.0254;
Cp_location(16,1) = 0.35*0.0254;
Cp_location(17,1) = 0.175*0.0254;


% Port location: Y values (converted to meters)
Cp_location(1,2) = 0.14665*0.0254;
Cp_location(2,2) = 0.33075*0.0254;
Cp_location(3,2) = 0.4018*0.0254;
Cp_location(4,2) = 0.476*0.0254;
Cp_location(5,2) = 0.49*0.0254;
Cp_location(6,2) = 0.4774*0.0254;
Cp_location(7,2) = 0.4403*0.0254;
Cp_location(8,2) = 0.38325*0.0254;
Cp_location(9,2) = 0.21875*0.0254;
Cp_location(10,2) = 0*0.0254;
Cp_location(11,2) = 0*0.0254;
Cp_location(12,2) = 0*0.0254;
Cp_location(13,2) = 0*0.0254;
Cp_location(14,2) = 0*0.0254;
Cp_location(15,2) = 0.0014*0.0254;
Cp_location(16,2) = 0.0175*0.0254;
Cp_location(17,2) = 0.03885*0.0254;

c = 0.08897800263; % Chord length (in meters)

% Normalize Cp_location:
Cp_location_norm = (Cp_location/c);



% Finding Cp for each port for each configuration
for i = 1:96
    Cp(1,i) = Pressure_Port_1(i)/(q_infinity(i));
    Cp(2,i) = Pressure_Port_2(i)/(q_infinity(i));
    Cp(3,i) = Pressure_Port_3(i)/(q_infinity(i));
    Cp(4,i) = Pressure_Port_4(i)/(q_infinity(i));
    Cp(5,i) = Pressure_Port_5(i)/(q_infinity(i));
    Cp(6,i) = Pressure_Port_6(i)/(q_infinity(i));
    Cp(7,i) = Pressure_Port_7(i)/(q_infinity(i));
    Cp(8,i) = Pressure_Port_8(i)/(q_infinity(i));
    Cp(9,i) = Pressure_Port_10(i)/(q_infinity(i));
    Cp(10,i) = Pressure_Port_11(i)/(q_infinity(i));
    Cp(11,i) = Pressure_Port_12(i)/(q_infinity(i));
    Cp(12,i) = Pressure_Port_14(i)/(q_infinity(i));
    Cp(13,i) = Pressure_Port_16(i)/(q_infinity(i));
    Cp(14,i) = Pressure_Port_17(i)/(q_infinity(i));
    Cp(15,i) = Pressure_Port_18(i)/(q_infinity(i));
    Cp(16,i) = Pressure_Port_19(i)/(q_infinity(i));
    Cp(17,i) = Pressure_Port_20(i)/(q_infinity(i));
    
    Cp_norm = (Cp/3.158);
    
   
end
f = figure(1);
hold on
 % Plotting Cp for each of 96 configurations
    for i = [30 90 57 24 48]
     plot(Cp_location_norm(:,1),Cp_norm(:,i));
     f.CurrentAxes.YDir = 'Reverse';
     legend('-5','0','5','10','12')
    end
hold off



for i = 1:96
    for j = 2:17
        Cn(j,i) = 0.5*(Cp(j-1,i) + Cp(j,i))*(Cp_location(j,1) - Cp_location(j-1,1))/c;
        Ca(j,i) = 0.5*(Cp(j-1,i) + Cp(j,i))*(Cp_location(j,2) - Cp_location(j-1,2))/c;
    end
    Cn_total(i) = -(sum(Cn(2:17,i)));
    Ca_total(i) = sum(Ca(2:17,i));
end


a = 1;
b = 20;
c = 1;
d = 20;
e = 1;
f = 20;
g = 1;
h = 20;
j = 1;
k = 20;
l = 1;
m = 20;
n = 1;
p = 20;
q = 1;
r = 20;

for i = 1:96
    if (i < 13)
        AOA(i) = mean(Data1(a:b,23));
        a = a + 20;
        b = b + 20;
    elseif (i < 25) && (i > 12)
        AOA(i) = mean(Data2(c:d,23));
        c = c + 20;
        d = d + 20;
    elseif (i > 24) && (i < 37)
        AOA(i) = mean(Data3(e:f,23));
        e = e + 20;
        f = f + 20;
    elseif (i > 36) && (i < 49)
        AOA(i) = mean(Data4(g:h,23));
        g = g + 20;
        h = h + 20;
    elseif (i > 48) && (i < 61)
        AOA(i) = mean(Data5(j:k,23));
        j = j + 20;
        k = k + 20;
    elseif (i > 60) && (i < 73)
        AOA(i) = mean(Data6(l:m,23));
        l = l + 20;
        m = m + 20;
    elseif (i > 72) && (i < 85)
        AOA(i) = mean(Data7(n:p,23));
        n = n + 20;
        p = p + 20;
    elseif (i > 84) && (i < 97)
        AOA(i) = mean(Data8(q:r,23));
        q = q + 20;
        r = r + 20;
    end
    
end

for i = 1:96
    Lift_Coef(i) = (Cn_total(i)*cosd(AOA(i)))-(Ca_total(i)*sind(AOA(i)));
    Drag_Coef(i) = (Cn_total(i)*sind(AOA(i)))+(Ca_total(i)*cosd(AOA(i)));

end


% CL vs AOA for V = 9.1 m/s

j = 1;
for i = 1:3:96
    CL_9(j) = Lift_Coef(i);
    CD_9(j) = Drag_Coef(i);
    AOA_9(j) = AOA(i);
    % q_infinity_9(j) = q_infinity(i);
    j = j + 1;
end

j = 1;
for i = 1:4:32
    AOA_9a(j) = AOA_9(i);
    CL_9a(j) = CL_9(i);
    CD_9a(j) = CD_9(i);
    j = j + 1;
end

j = 9;
for i = 2:4:32
    AOA_9a(j) = AOA_9(i);
    CL_9a(j) = CL_9(i);
    CD_9a(j) = CD_9(i);
    j = j + 1;
end

j = 17;
for i = 3:4:32
    AOA_9a(j) = AOA_9(i);
    CL_9a(j) = CL_9(i);
    CD_9a(j) = CD_9(i);
    j = j + 1;
end

j = 25;
for i = 4:4:32
    AOA_9a(j) = AOA_9(i);
    CL_9a(j) = CL_9(i);
    CD_9a(j) = CD_9(i);
    j = j + 1;
end




% CL vs AOA for V = 17 m/s

j = 1;
for i = 2:3:96
    CL_17(j) = Lift_Coef(i);
    CD_17(j) = Drag_Coef(i);
    AOA_17(j) = AOA(i);
    j = j + 1;
end

j = 1;
for i = 1:4:32
    AOA_17a(j) = AOA_17(i);
    CL_17a(j) = CL_17(i);
    CD_17a(j) = CD_17(i);
    j = j + 1;
end

j = 9;
for i = 2:4:32
    AOA_17a(j) = AOA_17(i);
    CL_17a(j) = CL_17(i);
    CD_17a(j) = CD_17(i);
    j = j + 1;
end

j = 17;
for i = 3:4:32
    AOA_17a(j) = AOA_17(i);
    CL_17a(j) = CL_17(i);
    CD_17a(j) = CD_17(i);
    j = j + 1;
end

j = 25;
for i = 4:4:32
    AOA_17a(j) = AOA_17(i);
    CL_17a(j) = CL_17(i);
    CD_17a(j) = CD_17(i);
    j = j + 1;
end



% CL vs AOA for V = 34 m/s

j = 1;
for i = 3:3:96
    CL_34(j) = Lift_Coef(i);
    CD_34(j) = Drag_Coef(i);
    AOA_34(j) = AOA(i);
    j = j + 1;
end

j = 1;
for i = 1:4:32
    AOA_34a(j) = AOA_34(i);
    CL_34a(j) = CL_34(i);
    CD_34a(j) = CD_34(i);
    j = j + 1;
end

j = 9;
for i = 2:4:32
    AOA_34a(j) = AOA_34(i);
    CL_34a(j) = CL_34(i);
    CD_34a(j) = CD_34(i);
    j = j + 1;
end

j = 17;
for i = 3:4:32
    AOA_34a(j) = AOA_34(i);
    CL_34a(j) = CL_34(i);
    CD_34a(j) = CD_34(i);
    j = j + 1;
end

j = 25;
for i = 4:4:32
    AOA_34a(j) = AOA_34(i);
    CL_34a(j) = CL_34(i);
    CD_34a(j) = CD_34(i);
    j = j + 1;
end

NACA_Cl = [-8 -0.12; -4 0.18; 0 0.48; 4 0.78; 8 1.08; 12 1.3; 16 1.52; 20 1.58; 24 1.45; 28 1.28];
NACA_Cd = [-8 0.015; -4 0.015; 0 0.02; 4 0.04; 8 0.08; 12 0.12; 16 0.17; 20 0.24; 24 0.36; 28 0.46];

figure(2);
plot(AOA_9a,CL_9a)
hold on
plot(AOA_17a,CL_17a)
plot(AOA_34a,CL_34a)
plot(NACA_Cl(:,1),NACA_Cl(:,2))
title('CL vs. AOA')
xlabel('Angle of Attack (in degrees)')
ylabel('Coefficient of Lift')
legend('9 m/s','17 m/s', '34 m/s','NACA Data')
hold off

figure(3);
plot(AOA_9a,CD_9a)
hold on
plot(AOA_17a,CD_17a)
plot(AOA_34a,CD_34a)
plot(NACA_Cd(:,1),NACA_Cd(:,2))
title('CD vs. AOA')
xlabel('Angle of Attack (in degrees)')
ylabel('Coefficient of Drag')
legend('9 m/s','17 m/s', '34 m/s','NACA Data')
hold off