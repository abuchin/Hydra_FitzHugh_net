
clear
close all

%% Numerics
dt=0.01;
T=500;


%% PARAMETERs

I=0.2;

%% ICs
V=-1.1994*ones(1,round(T/dt));
W=-0.6243*ones(1,round(T/dt));


%% TIME INTEGRATION LOOP
for t=1:1:round(T/dt)
    
% Voltage
dVdt=V(t)-(V(t)^3)/3-W(t)+I;
dWdt=0.08*(V(t)+0.7-0.8*W(t));

V(t+1)=dVdt*dt + V(t);
W(t+1)=dWdt*dt + W(t);

end


%% PLOT part
time=(1:1:round(T/dt))*dt;

subplot(2,1,1);
plot(time,V(1:end-1));
set(gca,'Fontsize',30);
title('V(t)');

subplot(2,1,2);
plot(time,W(1:end-1));
set(gca,'Fontsize',30);
title('W(t)');
xlabel('Time(ms)');
