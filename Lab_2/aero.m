
clear, clc, close all
n = 1;
airspdArr = [];
for i = 1:2:7
    vv{n} = readmatrix(sprintf('VelocityVoltage_S301_%d.csv',i));
    n = n + 1;
end

for i = 1:4
    R = 8.314;
    p_atm = vv{1,i}(:,1);
    T_atm = vv{1,i}(:,2);
    p_airspeeddiff = vv{1,i}(:,3);
    p_auxdiff = vv{1,i}(:,4);
    eldx = vv{1,i}(:,5);
    eldy = vv{1,i}(:,6);
    voltg = vv{1,i}(:,7);
    
    calcAirSpd = sqrt(p_airspeeddiff.*T_atm.*R./p_atm);
    for j = 1:5
        airspdArr = [airspdArr, mean(calcAirSpd([((j-1)*500)+1:j*500]))];
    end
end
airspdArr = sort(airspdArr);
voltArr = 0.5:0.5:10;
%voltArr = [1,3,5,7,9,2,4,6,8,10,1.5,3.5,5.5,7.5,9.5,0.5,2.5,4.5,6.5,8.5];
figure
plot(voltArr,airspdArr,'*-')
title('Airspeed vs Voltage for Pitot-Static and Airspeed Pressure Transducer')
xlabel('Voltage')
ylabel('Airspeed')
