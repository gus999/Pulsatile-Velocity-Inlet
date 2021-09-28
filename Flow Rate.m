% FFT formulation for aortic inlet 
ow rate
clear all;
clc;
%Period (72 strokes/min)
T=60/72;
%Diameter (m)
D=0.025;
%Transversal Area (m2)
Area=4.90875e-4;
%Density
ro=1050;
%dynamic viscosity (Pa)
mu=0.0035;
%conversion factor from Velocity to ml/s (multiplied by transversal area)
CF=(1e6)*Area;
%scale factor to obtain 70ml/stroke due to previous data is 63.7377 ml/s
scale= 1.09825;
time = T *[0:0.02:1]; % this interval must have the same correspondence with vel.
vel =scale*[0.00146551 0.0630498 0.129036 0.203825 0.27421 0.329928 0.390043 0.441359
0.478005 0.50585 0.523424 0.533664 0.535101 0.524801 0.507168 0.482201 0.439628 0.399991 0.353017 0.294308 0.228261 0.163685 0.0932393 0.0242605 0.001*ones(1,27)];
plot(time,vel,'ro','Linewidth',2)
xlim([0 T]);
xlabel('time (seconds)');
ylabel('velocity (m/s)');
grid on;
d = t(vel);
m = length(vel); M = floor((m + 1)=2); a0 = d(1)/m;
an = 2*real(d(2:M))/m; alast = d(M + 1)=m; bn = -2*imag(d(2:M))/m;
hold on
t = 0:0.01:2;
n = 1:length(an);
y = a0 + an*cos(2*pi*n'*t/T) ...
+ bn*sin(2*pi*n'*t/T) ...
+ alast*cos(2*pi*(length(an)+1)*t/T);
plot(t,y,'Linewidth',2)
legend('Data','DFT Interpolant')
title('Velocity through ascending aorta inlet');
% Print out formulation
fprintf(['Volume flow rate through aorta inlet/n Interpolant: v=' num2str(a0)]);
% Integral over period
mlPerBeat = CF*(a0*(T)-bn*(ones(size(n))./(2*pi*n'/T)')'+ bn*(ones(size(n))./(2*pi*n'/T)')');
litresPerMinute = mlPerBeat/T*60/1000;
Vmax=max (vel);
Vmin=min(vel);
Reynolds=