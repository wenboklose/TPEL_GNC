clear all;
% close all;
%%
% when using matlab linearization tool, better first simulate the circuit
% to find operation point, then linerize it, some circuit's operation point
% is related to operation point (AFE), some not (VSI).
R=15;
L=1.1e-3; RL=150e-3;Rl=150e-3; 
C=41.3e-6;
Cdc=0.15e-3;
Rcdc=1e-6;
Ldc=0.99e-3;
Rsdc=1e-6;

f=400; 
w=2*pi*f; 
fsw=20e+3; 
Tsw=1/fsw; 
Tdelay = 0.5*Tsw;
Tdelay1= 0.6*Tsw;
wsw=2*pi*fsw;
DT=1e-6;

J=[0 -1; 1 0];

Vdc=270; Vdcref=270;
Vse=57.5;
Vsm=Vse*sqrt(2);
Vsdq=[sqrt(3/2)*Vsm; 0];
Vdref=Vsdq(1); Vqref=Vsdq(2);
omega = 400*2*pi;
Rac =  R;
Vod = Vdref;
Voq = Vqref;
% Rl=RL;
Dd = -(omega*L*Voq-Rac*Vod-Rl*Vod+Rl*omega*C*Voq*Rac+omega^2*L*C*Vod*Rac)/(Rac*Vdc)
Dq = -(-Rl*Voq-Rac*Voq-omega*L*Vod+omega^2*L*C*Voq*Rac-Rl*omega*C*Vod*Rac)/(Rac*Vdc)
Idref = -(-Vod+omega*C*Voq*Rac)/Rac
Iqref = (Voq+omega*C*Vod*Rac)/Rac

fv=200; fi=800;

kpi=2*pi*fi*L/Vdcref*0.5
kii=2*pi*fi*(Rl)/Vdcref*13
kpv=2*pi*fv*C
kiv=2*pi*fv/R

% new controller
