%% AFE parameters and Control to output transfer functions:
clc
clear all
close all

P=bodeoptions;
P.XLabel.FontSize=16;
P.YLabel.FontSize=16;
P.TickLabel.FontSize=16;
P.Title.FontSize=16;
P.title.String=' ';
P.Grid='on';
P.XLim={[20 1e4]};
P.XLimMode={'manual'};
P.FreqUnits='Hz';
P.PhaseWrapping='off';

%% for phase to netural derivation
s=tf([1,0],[0,1]);
omega=400*2*pi;
L=1200e-6;
Vdc=270;
Vod=99.6;
Voq=0;
Rac=15;
C=22e-6;
Cdc=150e-6;
Rl=0.1;
Iq=0;
fsw=20e3;
%% old controller
fv=100; 
fi=1000;
Kpi=2*pi*fi*L/Vdc; 
Kii=2*pi*fi*(Rl)/Vdc;
%% new controller
Kpi=94.221*0.00022
Kii=94.221
% Kpi=0.0115; 
% Kii=217.3372;
% Kpv=2*pi*fv*C;
% Kiv=2*pi*fv/Rac;
Tdelay = 1/fsw;
tfdelay = (1-Tdelay*s)/(1+Tdelay*s);
% Rac = 15e10;
tfididref = (C^2*Kpi*Rac^2*Vdc*s^3+C^2*Kpi*Rac^2*Vdc*s*omega^2+C^2*L*Rac^2*s^4+C^2*L*Rac^2*s^2*omega^2+C^2*Kii*Rac^2*Vdc*s^2+C^2*Kii*Rac^2*Vdc*omega^2+C^2*Rac^2*Rl*s^3+C^2*Rac^2*Rl*s*omega^2+2*Kpi*Vdc*s^2*C*Rac+2*C*L*Rac*s^3+2*C*Kii*Rac*Vdc*s+C*Rac^2*s^2+2*s^2*C*Rac*Rl+s*Vdc*Kpi+s^2*L+Kii*Vdc+s*Rac+s*Rl)*Vdc*(s*Kpi+Kii)/(C^2*Kpi^2*Rac^2*Vdc^2*s^4+C^2*Kpi^2*Rac^2*Vdc^2*s^2*omega^2+2*C^2*Kpi*L*Rac^2*Vdc*s^5+2*C^2*Kpi*L*Rac^2*Vdc*s^3*omega^2+C^2*L^2*Rac^2*s^6+2*C^2*L^2*Rac^2*s^4*omega^2+C^2*L^2*Rac^2*s^2*omega^4+2*C^2*Kii*Kpi*Rac^2*Vdc^2*s^3+2*C^2*Kii*Kpi*Rac^2*Vdc^2*s*omega^2+2*C^2*Kii*L*Rac^2*Vdc*s^4+2*C^2*Kii*L*Rac^2*Vdc*s^2*omega^2+2*C^2*Kpi*Rac^2*Rl*Vdc*s^4+2*C^2*Kpi*Rac^2*Rl*Vdc*s^2*omega^2+2*C^2*L*Rac^2*Rl*s^5+2*C^2*L*Rac^2*Rl*s^3*omega^2+C^2*Kii^2*Rac^2*Vdc^2*s^2+C^2*Kii^2*Rac^2*Vdc^2*omega^2+2*C^2*Kii*Rac^2*Rl*Vdc*s^3+2*C^2*Kii*Rac^2*Rl*Vdc*s*omega^2+C^2*Rac^2*Rl^2*s^4+C^2*Rac^2*Rl^2*s^2*omega^2+2*C*Kpi^2*Rac*Vdc^2*s^3+4*C*Kpi*L*Rac*Vdc*s^4+2*C*L^2*Rac*s^5+2*C*L^2*Rac*s^3*omega^2+4*C*Kii*Kpi*Rac*Vdc^2*s^2+4*C*Kii*L*Rac*Vdc*s^3+2*C*Kpi*Rac^2*Vdc*s^3+4*C*Kpi*Rac*Rl*Vdc*s^3+2*C*L*Rac^2*s^4-2*C*L*Rac^2*s^2*omega^2+4*C*L*Rac*Rl*s^4+2*C*Kii^2*Rac*Vdc^2*s+2*C*Kii*Rac^2*Vdc*s^2+4*C*Kii*Rac*Rl*Vdc*s^2+2*C*Rac^2*Rl*s^3+2*C*Rac*Rl^2*s^3+Kpi^2*Vdc^2*s^2+2*Kpi*L*Vdc*s^3+L^2*s^4+L^2*s^2*omega^2+2*Kii*Kpi*Vdc^2*s+2*Kii*L*Vdc*s^2+2*Kpi*Rac*Vdc*s^2+2*Kpi*Rl*Vdc*s^2+2*Rac*Rl*s^2+2*L*Rac*s^3+2*L*Rl*s^3+Kii^2*Vdc^2+2*Kii*Rac*Vdc*s+2*Kii*Rl*Vdc*s+Rl^2*s^2+Rac^2*s^2);

tfidiqref = (C^2*L*Rac^2*s^2+C^2*L*Rac^2*omega^2-C*Rac^2+2*L*s*C*Rac+L)*s*omega*Vdc*(s*Kpi+Kii)/(C^2*Kpi^2*Rac^2*Vdc^2*s^4+C^2*Kpi^2*Rac^2*Vdc^2*s^2*omega^2+2*C^2*Kpi*L*Rac^2*Vdc*s^5+2*C^2*Kpi*L*Rac^2*Vdc*s^3*omega^2+C^2*L^2*Rac^2*s^6+2*C^2*L^2*Rac^2*s^4*omega^2+C^2*L^2*Rac^2*s^2*omega^4+2*C^2*Kii*Kpi*Rac^2*Vdc^2*s^3+2*C^2*Kii*Kpi*Rac^2*Vdc^2*s*omega^2+2*C^2*Kii*L*Rac^2*Vdc*s^4+2*C^2*Kii*L*Rac^2*Vdc*s^2*omega^2+2*C^2*Kpi*Rac^2*Rl*Vdc*s^4+2*C^2*Kpi*Rac^2*Rl*Vdc*s^2*omega^2+2*C^2*L*Rac^2*Rl*s^5+2*C^2*L*Rac^2*Rl*s^3*omega^2+C^2*Kii^2*Rac^2*Vdc^2*s^2+C^2*Kii^2*Rac^2*Vdc^2*omega^2+2*C^2*Kii*Rac^2*Rl*Vdc*s^3+2*C^2*Kii*Rac^2*Rl*Vdc*s*omega^2+C^2*Rac^2*Rl^2*s^4+C^2*Rac^2*Rl^2*s^2*omega^2+2*C*Kpi^2*Rac*Vdc^2*s^3+4*C*Kpi*L*Rac*Vdc*s^4+2*C*L^2*Rac*s^5+2*C*L^2*Rac*s^3*omega^2+4*C*Kii*Kpi*Rac*Vdc^2*s^2+4*C*Kii*L*Rac*Vdc*s^3+2*C*Kpi*Rac^2*Vdc*s^3+4*C*Kpi*Rac*Rl*Vdc*s^3+2*C*L*Rac^2*s^4-2*C*L*Rac^2*s^2*omega^2+4*C*L*Rac*Rl*s^4+2*C*Kii^2*Rac*Vdc^2*s+2*C*Kii*Rac^2*Vdc*s^2+4*C*Kii*Rac*Rl*Vdc*s^2+2*C*Rac^2*Rl*s^3+2*C*Rac*Rl^2*s^3+Kpi^2*Vdc^2*s^2+2*Kpi*L*Vdc*s^3+L^2*s^4+L^2*s^2*omega^2+2*Kii*Kpi*Vdc^2*s+2*Kii*L*Vdc*s^2+2*Kpi*Rac*Vdc*s^2+2*Kpi*Rl*Vdc*s^2+2*Rac*Rl*s^2+2*L*Rac*s^3+2*L*Rl*s^3+Kii^2*Vdc^2+2*Kii*Rac*Vdc*s+2*Kii*Rl*Vdc*s+Rl^2*s^2+Rac^2*s^2);

tfiqidref = -(C^2*L*Rac^2*s^2+C^2*L*Rac^2*omega^2-C*Rac^2+2*L*s*C*Rac+L)*s*omega*Vdc*(s*Kpi+Kii)/(C^2*Kpi^2*Rac^2*Vdc^2*s^4+C^2*Kpi^2*Rac^2*Vdc^2*s^2*omega^2+2*C^2*Kpi*L*Rac^2*Vdc*s^5+2*C^2*Kpi*L*Rac^2*Vdc*s^3*omega^2+C^2*L^2*Rac^2*s^6+2*C^2*L^2*Rac^2*s^4*omega^2+C^2*L^2*Rac^2*s^2*omega^4+2*C^2*Kii*Kpi*Rac^2*Vdc^2*s^3+2*C^2*Kii*Kpi*Rac^2*Vdc^2*s*omega^2+2*C^2*Kii*L*Rac^2*Vdc*s^4+2*C^2*Kii*L*Rac^2*Vdc*s^2*omega^2+2*C^2*Kpi*Rac^2*Rl*Vdc*s^4+2*C^2*Kpi*Rac^2*Rl*Vdc*s^2*omega^2+2*C^2*L*Rac^2*Rl*s^5+2*C^2*L*Rac^2*Rl*s^3*omega^2+C^2*Kii^2*Rac^2*Vdc^2*s^2+C^2*Kii^2*Rac^2*Vdc^2*omega^2+2*C^2*Kii*Rac^2*Rl*Vdc*s^3+2*C^2*Kii*Rac^2*Rl*Vdc*s*omega^2+C^2*Rac^2*Rl^2*s^4+C^2*Rac^2*Rl^2*s^2*omega^2+2*C*Kpi^2*Rac*Vdc^2*s^3+4*C*Kpi*L*Rac*Vdc*s^4+2*C*L^2*Rac*s^5+2*C*L^2*Rac*s^3*omega^2+4*C*Kii*Kpi*Rac*Vdc^2*s^2+4*C*Kii*L*Rac*Vdc*s^3+2*C*Kpi*Rac^2*Vdc*s^3+4*C*Kpi*Rac*Rl*Vdc*s^3+2*C*L*Rac^2*s^4-2*C*L*Rac^2*s^2*omega^2+4*C*L*Rac*Rl*s^4+2*C*Kii^2*Rac*Vdc^2*s+2*C*Kii*Rac^2*Vdc*s^2+4*C*Kii*Rac*Rl*Vdc*s^2+2*C*Rac^2*Rl*s^3+2*C*Rac*Rl^2*s^3+Kpi^2*Vdc^2*s^2+2*Kpi*L*Vdc*s^3+L^2*s^4+L^2*s^2*omega^2+2*Kii*Kpi*Vdc^2*s+2*Kii*L*Vdc*s^2+2*Kpi*Rac*Vdc*s^2+2*Kpi*Rl*Vdc*s^2+2*Rac*Rl*s^2+2*L*Rac*s^3+2*L*Rl*s^3+Kii^2*Vdc^2+2*Kii*Rac*Vdc*s+2*Kii*Rl*Vdc*s+Rl^2*s^2+Rac^2*s^2);

tfiqiqref = (C^2*Kpi*Rac^2*Vdc*s^3+C^2*Kpi*Rac^2*Vdc*s*omega^2+C^2*L*Rac^2*s^4+C^2*L*Rac^2*s^2*omega^2+C^2*Kii*Rac^2*Vdc*s^2+C^2*Kii*Rac^2*Vdc*omega^2+C^2*Rac^2*Rl*s^3+C^2*Rac^2*Rl*s*omega^2+2*Kpi*Vdc*s^2*C*Rac+2*C*L*Rac*s^3+2*C*Kii*Rac*Vdc*s+C*Rac^2*s^2+2*s^2*C*Rac*Rl+s*Vdc*Kpi+s^2*L+Kii*Vdc+s*Rac+s*Rl)*Vdc*(s*Kpi+Kii)/(C^2*Kpi^2*Rac^2*Vdc^2*s^4+C^2*Kpi^2*Rac^2*Vdc^2*s^2*omega^2+2*C^2*Kpi*L*Rac^2*Vdc*s^5+2*C^2*Kpi*L*Rac^2*Vdc*s^3*omega^2+C^2*L^2*Rac^2*s^6+2*C^2*L^2*Rac^2*s^4*omega^2+C^2*L^2*Rac^2*s^2*omega^4+2*C^2*Kii*Kpi*Rac^2*Vdc^2*s^3+2*C^2*Kii*Kpi*Rac^2*Vdc^2*s*omega^2+2*C^2*Kii*L*Rac^2*Vdc*s^4+2*C^2*Kii*L*Rac^2*Vdc*s^2*omega^2+2*C^2*Kpi*Rac^2*Rl*Vdc*s^4+2*C^2*Kpi*Rac^2*Rl*Vdc*s^2*omega^2+2*C^2*L*Rac^2*Rl*s^5+2*C^2*L*Rac^2*Rl*s^3*omega^2+C^2*Kii^2*Rac^2*Vdc^2*s^2+C^2*Kii^2*Rac^2*Vdc^2*omega^2+2*C^2*Kii*Rac^2*Rl*Vdc*s^3+2*C^2*Kii*Rac^2*Rl*Vdc*s*omega^2+C^2*Rac^2*Rl^2*s^4+C^2*Rac^2*Rl^2*s^2*omega^2+2*C*Kpi^2*Rac*Vdc^2*s^3+4*C*Kpi*L*Rac*Vdc*s^4+2*C*L^2*Rac*s^5+2*C*L^2*Rac*s^3*omega^2+4*C*Kii*Kpi*Rac*Vdc^2*s^2+4*C*Kii*L*Rac*Vdc*s^3+2*C*Kpi*Rac^2*Vdc*s^3+4*C*Kpi*Rac*Rl*Vdc*s^3+2*C*L*Rac^2*s^4-2*C*L*Rac^2*s^2*omega^2+4*C*L*Rac*Rl*s^4+2*C*Kii^2*Rac*Vdc^2*s+2*C*Kii*Rac^2*Vdc*s^2+4*C*Kii*Rac*Rl*Vdc*s^2+2*C*Rac^2*Rl*s^3+2*C*Rac*Rl^2*s^3+Kpi^2*Vdc^2*s^2+2*Kpi*L*Vdc*s^3+L^2*s^4+L^2*s^2*omega^2+2*Kii*Kpi*Vdc^2*s+2*Kii*L*Vdc*s^2+2*Kpi*Rac*Vdc*s^2+2*Kpi*Rl*Vdc*s^2+2*Rac*Rl*s^2+2*L*Rac*s^3+2*L*Rl*s^3+Kii^2*Vdc^2+2*Kii*Rac*Vdc*s+2*Kii*Rl*Vdc*s+Rl^2*s^2+Rac^2*s^2);

tfvdidref = (Kpi*Vdc*s^2*C*Rac+C*L*Rac*s^3-C*L*Rac*s*omega^2+C*Kii*Rac*Vdc*s+s^2*C*Rac*Rl+s*Vdc*Kpi+s^2*L+Kii*Vdc+s*Rac+s*Rl)*Rac*Vdc*(s*Kpi+Kii)/(C^2*Kpi^2*Rac^2*Vdc^2*s^4+C^2*Kpi^2*Rac^2*Vdc^2*s^2*omega^2+2*C^2*Kpi*L*Rac^2*Vdc*s^5+2*C^2*Kpi*L*Rac^2*Vdc*s^3*omega^2+C^2*L^2*Rac^2*s^6+2*C^2*L^2*Rac^2*s^4*omega^2+C^2*L^2*Rac^2*s^2*omega^4+2*C^2*Kii*Kpi*Rac^2*Vdc^2*s^3+2*C^2*Kii*Kpi*Rac^2*Vdc^2*s*omega^2+2*C^2*Kii*L*Rac^2*Vdc*s^4+2*C^2*Kii*L*Rac^2*Vdc*s^2*omega^2+2*C^2*Kpi*Rac^2*Rl*Vdc*s^4+2*C^2*Kpi*Rac^2*Rl*Vdc*s^2*omega^2+2*C^2*L*Rac^2*Rl*s^5+2*C^2*L*Rac^2*Rl*s^3*omega^2+C^2*Kii^2*Rac^2*Vdc^2*s^2+C^2*Kii^2*Rac^2*Vdc^2*omega^2+2*C^2*Kii*Rac^2*Rl*Vdc*s^3+2*C^2*Kii*Rac^2*Rl*Vdc*s*omega^2+C^2*Rac^2*Rl^2*s^4+C^2*Rac^2*Rl^2*s^2*omega^2+2*C*Kpi^2*Rac*Vdc^2*s^3+4*C*Kpi*L*Rac*Vdc*s^4+2*C*L^2*Rac*s^5+2*C*L^2*Rac*s^3*omega^2+4*C*Kii*Kpi*Rac*Vdc^2*s^2+4*C*Kii*L*Rac*Vdc*s^3+2*C*Kpi*Rac^2*Vdc*s^3+4*C*Kpi*Rac*Rl*Vdc*s^3+2*C*L*Rac^2*s^4-2*C*L*Rac^2*s^2*omega^2+4*C*L*Rac*Rl*s^4+2*C*Kii^2*Rac*Vdc^2*s+2*C*Kii*Rac^2*Vdc*s^2+4*C*Kii*Rac*Rl*Vdc*s^2+2*C*Rac^2*Rl*s^3+2*C*Rac*Rl^2*s^3+Kpi^2*Vdc^2*s^2+2*Kpi*L*Vdc*s^3+L^2*s^4+L^2*s^2*omega^2+2*Kii*Kpi*Vdc^2*s+2*Kii*L*Vdc*s^2+2*Kpi*Rac*Vdc*s^2+2*Kpi*Rl*Vdc*s^2+2*Rac*Rl*s^2+2*L*Rac*s^3+2*L*Rl*s^3+Kii^2*Vdc^2+2*Kii*Rac*Vdc*s+2*Kii*Rl*Vdc*s+Rl^2*s^2+Rac^2*s^2);

tfvdiqref = (s*C*Rac*Vdc*Kpi+2*s^2*C*Rac*L+C*Kii*Rac*Vdc+s*C*Rac*Rl+s*L)*Vdc*omega*Rac*(s*Kpi+Kii)/(C^2*Kpi^2*Rac^2*Vdc^2*s^4+C^2*Kpi^2*Rac^2*Vdc^2*s^2*omega^2+2*C^2*Kpi*L*Rac^2*Vdc*s^5+2*C^2*Kpi*L*Rac^2*Vdc*s^3*omega^2+C^2*L^2*Rac^2*s^6+2*C^2*L^2*Rac^2*s^4*omega^2+C^2*L^2*Rac^2*s^2*omega^4+2*C^2*Kii*Kpi*Rac^2*Vdc^2*s^3+2*C^2*Kii*Kpi*Rac^2*Vdc^2*s*omega^2+2*C^2*Kii*L*Rac^2*Vdc*s^4+2*C^2*Kii*L*Rac^2*Vdc*s^2*omega^2+2*C^2*Kpi*Rac^2*Rl*Vdc*s^4+2*C^2*Kpi*Rac^2*Rl*Vdc*s^2*omega^2+2*C^2*L*Rac^2*Rl*s^5+2*C^2*L*Rac^2*Rl*s^3*omega^2+C^2*Kii^2*Rac^2*Vdc^2*s^2+C^2*Kii^2*Rac^2*Vdc^2*omega^2+2*C^2*Kii*Rac^2*Rl*Vdc*s^3+2*C^2*Kii*Rac^2*Rl*Vdc*s*omega^2+C^2*Rac^2*Rl^2*s^4+C^2*Rac^2*Rl^2*s^2*omega^2+2*C*Kpi^2*Rac*Vdc^2*s^3+4*C*Kpi*L*Rac*Vdc*s^4+2*C*L^2*Rac*s^5+2*C*L^2*Rac*s^3*omega^2+4*C*Kii*Kpi*Rac*Vdc^2*s^2+4*C*Kii*L*Rac*Vdc*s^3+2*C*Kpi*Rac^2*Vdc*s^3+4*C*Kpi*Rac*Rl*Vdc*s^3+2*C*L*Rac^2*s^4-2*C*L*Rac^2*s^2*omega^2+4*C*L*Rac*Rl*s^4+2*C*Kii^2*Rac*Vdc^2*s+2*C*Kii*Rac^2*Vdc*s^2+4*C*Kii*Rac*Rl*Vdc*s^2+2*C*Rac^2*Rl*s^3+2*C*Rac*Rl^2*s^3+Kpi^2*Vdc^2*s^2+2*Kpi*L*Vdc*s^3+L^2*s^4+L^2*s^2*omega^2+2*Kii*Kpi*Vdc^2*s+2*Kii*L*Vdc*s^2+2*Kpi*Rac*Vdc*s^2+2*Kpi*Rl*Vdc*s^2+2*Rac*Rl*s^2+2*L*Rac*s^3+2*L*Rl*s^3+Kii^2*Vdc^2+2*Kii*Rac*Vdc*s+2*Kii*Rl*Vdc*s+Rl^2*s^2+Rac^2*s^2);

tfvqidref = -(s*C*Rac*Vdc*Kpi+2*s^2*C*Rac*L+C*Kii*Rac*Vdc+s*C*Rac*Rl+s*L)*Vdc*omega*Rac*(s*Kpi+Kii)/(C^2*Kpi^2*Rac^2*Vdc^2*s^4+C^2*Kpi^2*Rac^2*Vdc^2*s^2*omega^2+2*C^2*Kpi*L*Rac^2*Vdc*s^5+2*C^2*Kpi*L*Rac^2*Vdc*s^3*omega^2+C^2*L^2*Rac^2*s^6+2*C^2*L^2*Rac^2*s^4*omega^2+C^2*L^2*Rac^2*s^2*omega^4+2*C^2*Kii*Kpi*Rac^2*Vdc^2*s^3+2*C^2*Kii*Kpi*Rac^2*Vdc^2*s*omega^2+2*C^2*Kii*L*Rac^2*Vdc*s^4+2*C^2*Kii*L*Rac^2*Vdc*s^2*omega^2+2*C^2*Kpi*Rac^2*Rl*Vdc*s^4+2*C^2*Kpi*Rac^2*Rl*Vdc*s^2*omega^2+2*C^2*L*Rac^2*Rl*s^5+2*C^2*L*Rac^2*Rl*s^3*omega^2+C^2*Kii^2*Rac^2*Vdc^2*s^2+C^2*Kii^2*Rac^2*Vdc^2*omega^2+2*C^2*Kii*Rac^2*Rl*Vdc*s^3+2*C^2*Kii*Rac^2*Rl*Vdc*s*omega^2+C^2*Rac^2*Rl^2*s^4+C^2*Rac^2*Rl^2*s^2*omega^2+2*C*Kpi^2*Rac*Vdc^2*s^3+4*C*Kpi*L*Rac*Vdc*s^4+2*C*L^2*Rac*s^5+2*C*L^2*Rac*s^3*omega^2+4*C*Kii*Kpi*Rac*Vdc^2*s^2+4*C*Kii*L*Rac*Vdc*s^3+2*C*Kpi*Rac^2*Vdc*s^3+4*C*Kpi*Rac*Rl*Vdc*s^3+2*C*L*Rac^2*s^4-2*C*L*Rac^2*s^2*omega^2+4*C*L*Rac*Rl*s^4+2*C*Kii^2*Rac*Vdc^2*s+2*C*Kii*Rac^2*Vdc*s^2+4*C*Kii*Rac*Rl*Vdc*s^2+2*C*Rac^2*Rl*s^3+2*C*Rac*Rl^2*s^3+Kpi^2*Vdc^2*s^2+2*Kpi*L*Vdc*s^3+L^2*s^4+L^2*s^2*omega^2+2*Kii*Kpi*Vdc^2*s+2*Kii*L*Vdc*s^2+2*Kpi*Rac*Vdc*s^2+2*Kpi*Rl*Vdc*s^2+2*Rac*Rl*s^2+2*L*Rac*s^3+2*L*Rl*s^3+Kii^2*Vdc^2+2*Kii*Rac*Vdc*s+2*Kii*Rl*Vdc*s+Rl^2*s^2+Rac^2*s^2);

tfvqiqref = (Kpi*Vdc*s^2*C*Rac+C*L*Rac*s^3-C*L*Rac*s*omega^2+C*Kii*Rac*Vdc*s+s^2*C*Rac*Rl+s*Vdc*Kpi+s^2*L+Kii*Vdc+s*Rac+s*Rl)*Rac*Vdc*(s*Kpi+Kii)/(C^2*Kpi^2*Rac^2*Vdc^2*s^4+C^2*Kpi^2*Rac^2*Vdc^2*s^2*omega^2+2*C^2*Kpi*L*Rac^2*Vdc*s^5+2*C^2*Kpi*L*Rac^2*Vdc*s^3*omega^2+C^2*L^2*Rac^2*s^6+2*C^2*L^2*Rac^2*s^4*omega^2+C^2*L^2*Rac^2*s^2*omega^4+2*C^2*Kii*Kpi*Rac^2*Vdc^2*s^3+2*C^2*Kii*Kpi*Rac^2*Vdc^2*s*omega^2+2*C^2*Kii*L*Rac^2*Vdc*s^4+2*C^2*Kii*L*Rac^2*Vdc*s^2*omega^2+2*C^2*Kpi*Rac^2*Rl*Vdc*s^4+2*C^2*Kpi*Rac^2*Rl*Vdc*s^2*omega^2+2*C^2*L*Rac^2*Rl*s^5+2*C^2*L*Rac^2*Rl*s^3*omega^2+C^2*Kii^2*Rac^2*Vdc^2*s^2+C^2*Kii^2*Rac^2*Vdc^2*omega^2+2*C^2*Kii*Rac^2*Rl*Vdc*s^3+2*C^2*Kii*Rac^2*Rl*Vdc*s*omega^2+C^2*Rac^2*Rl^2*s^4+C^2*Rac^2*Rl^2*s^2*omega^2+2*C*Kpi^2*Rac*Vdc^2*s^3+4*C*Kpi*L*Rac*Vdc*s^4+2*C*L^2*Rac*s^5+2*C*L^2*Rac*s^3*omega^2+4*C*Kii*Kpi*Rac*Vdc^2*s^2+4*C*Kii*L*Rac*Vdc*s^3+2*C*Kpi*Rac^2*Vdc*s^3+4*C*Kpi*Rac*Rl*Vdc*s^3+2*C*L*Rac^2*s^4-2*C*L*Rac^2*s^2*omega^2+4*C*L*Rac*Rl*s^4+2*C*Kii^2*Rac*Vdc^2*s+2*C*Kii*Rac^2*Vdc*s^2+4*C*Kii*Rac*Rl*Vdc*s^2+2*C*Rac^2*Rl*s^3+2*C*Rac*Rl^2*s^3+Kpi^2*Vdc^2*s^2+2*Kpi*L*Vdc*s^3+L^2*s^4+L^2*s^2*omega^2+2*Kii*Kpi*Vdc^2*s+2*Kii*L*Vdc*s^2+2*Kpi*Rac*Vdc*s^2+2*Kpi*Rl*Vdc*s^2+2*Rac*Rl*s^2+2*L*Rac*s^3+2*L*Rl*s^3+Kii^2*Vdc^2+2*Kii*Rac*Vdc*s+2*Kii*Rl*Vdc*s+Rl^2*s^2+Rac^2*s^2);
Rac = 15e10;
Zoutddcl = (C*Kpi^2*Rac*Vdc^2*s^3+2*C*Kpi*L*Rac*Vdc*s^4+C*L^2*Rac*s^5+C*L^2*Rac*s^3*omega^2+2*C*Kii*Kpi*Rac*Vdc^2*s^2+2*C*Kii*L*Rac*Vdc*s^3+2*C*Kpi*Rac*Rl*Vdc*s^3+2*C*L*Rac*Rl*s^4+C*Kii^2*Rac*Vdc^2*s+2*C*Kii*Rac*Rl*Vdc*s^2+C*Rac*Rl^2*s^3+Kpi^2*Vdc^2*s^2+2*Kpi*L*Vdc*s^3+L^2*s^4+L^2*s^2*omega^2+2*Kii*Kpi*Vdc^2*s+2*Kii*L*Vdc*s^2+Kpi*Rac*Vdc*s^2+2*Kpi*Rl*Vdc*s^2+Rac*Rl*s^2+L*Rac*s^3+2*L*Rl*s^3+Kii^2*Vdc^2+Kii*Rac*Vdc*s+2*Kii*Rl*Vdc*s+Rl^2*s^2)*Rac/(C^2*Kpi^2*Rac^2*Vdc^2*s^4+C^2*Kpi^2*Rac^2*Vdc^2*s^2*omega^2+2*C^2*Kpi*L*Rac^2*Vdc*s^5+2*C^2*Kpi*L*Rac^2*Vdc*s^3*omega^2+C^2*L^2*Rac^2*s^6+2*C^2*L^2*Rac^2*s^4*omega^2+C^2*L^2*Rac^2*s^2*omega^4+2*C^2*Kii*Kpi*Rac^2*Vdc^2*s^3+2*C^2*Kii*Kpi*Rac^2*Vdc^2*s*omega^2+2*C^2*Kii*L*Rac^2*Vdc*s^4+2*C^2*Kii*L*Rac^2*Vdc*s^2*omega^2+2*C^2*Kpi*Rac^2*Rl*Vdc*s^4+2*C^2*Kpi*Rac^2*Rl*Vdc*s^2*omega^2+2*C^2*L*Rac^2*Rl*s^5+2*C^2*L*Rac^2*Rl*s^3*omega^2+C^2*Kii^2*Rac^2*Vdc^2*s^2+C^2*Kii^2*Rac^2*Vdc^2*omega^2+2*C^2*Kii*Rac^2*Rl*Vdc*s^3+2*C^2*Kii*Rac^2*Rl*Vdc*s*omega^2+C^2*Rac^2*Rl^2*s^4+C^2*Rac^2*Rl^2*s^2*omega^2+2*C*Kpi^2*Rac*Vdc^2*s^3+4*C*Kpi*L*Rac*Vdc*s^4+2*C*L^2*Rac*s^5+2*C*L^2*Rac*s^3*omega^2+4*C*Kii*Kpi*Rac*Vdc^2*s^2+4*C*Kii*L*Rac*Vdc*s^3+2*C*Kpi*Rac^2*Vdc*s^3+4*C*Kpi*Rac*Rl*Vdc*s^3+2*C*L*Rac^2*s^4-2*C*L*Rac^2*s^2*omega^2+4*C*L*Rac*Rl*s^4+2*C*Kii^2*Rac*Vdc^2*s+2*C*Kii*Rac^2*Vdc*s^2+4*C*Kii*Rac*Rl*Vdc*s^2+2*C*Rac^2*Rl*s^3+2*C*Rac*Rl^2*s^3+Kpi^2*Vdc^2*s^2+2*Kpi*L*Vdc*s^3+L^2*s^4+L^2*s^2*omega^2+2*Kii*Kpi*Vdc^2*s+2*Kii*L*Vdc*s^2+2*Kpi*Rac*Vdc*s^2+2*Kpi*Rl*Vdc*s^2+2*Rac*Rl*s^2+2*L*Rac*s^3+2*L*Rl*s^3+Kii^2*Vdc^2+2*Kii*Rac*Vdc*s+2*Kii*Rl*Vdc*s+Rl^2*s^2+Rac^2*s^2);

Zoutdqcl = omega*Rac^2*(C*Kpi^2*Vdc^2*s^2+2*C*Kpi*L*Vdc*s^3+C*L^2*s^4+C*L^2*s^2*omega^2+2*C*Kii*Kpi*Vdc^2*s+2*C*Kii*L*Vdc*s^2+2*C*Kpi*Rl*Vdc*s^2+2*C*L*Rl*s^3+C*Kii^2*Vdc^2+2*C*Kii*Rl*Vdc*s+C*Rl^2*s^2-s^2*L)/(C^2*Kpi^2*Rac^2*Vdc^2*s^4+C^2*Kpi^2*Rac^2*Vdc^2*s^2*omega^2+2*C^2*Kpi*L*Rac^2*Vdc*s^5+2*C^2*Kpi*L*Rac^2*Vdc*s^3*omega^2+C^2*L^2*Rac^2*s^6+2*C^2*L^2*Rac^2*s^4*omega^2+C^2*L^2*Rac^2*s^2*omega^4+2*C^2*Kii*Kpi*Rac^2*Vdc^2*s^3+2*C^2*Kii*Kpi*Rac^2*Vdc^2*s*omega^2+2*C^2*Kii*L*Rac^2*Vdc*s^4+2*C^2*Kii*L*Rac^2*Vdc*s^2*omega^2+2*C^2*Kpi*Rac^2*Rl*Vdc*s^4+2*C^2*Kpi*Rac^2*Rl*Vdc*s^2*omega^2+2*C^2*L*Rac^2*Rl*s^5+2*C^2*L*Rac^2*Rl*s^3*omega^2+C^2*Kii^2*Rac^2*Vdc^2*s^2+C^2*Kii^2*Rac^2*Vdc^2*omega^2+2*C^2*Kii*Rac^2*Rl*Vdc*s^3+2*C^2*Kii*Rac^2*Rl*Vdc*s*omega^2+C^2*Rac^2*Rl^2*s^4+C^2*Rac^2*Rl^2*s^2*omega^2+2*C*Kpi^2*Rac*Vdc^2*s^3+4*C*Kpi*L*Rac*Vdc*s^4+2*C*L^2*Rac*s^5+2*C*L^2*Rac*s^3*omega^2+4*C*Kii*Kpi*Rac*Vdc^2*s^2+4*C*Kii*L*Rac*Vdc*s^3+2*C*Kpi*Rac^2*Vdc*s^3+4*C*Kpi*Rac*Rl*Vdc*s^3+2*C*L*Rac^2*s^4-2*C*L*Rac^2*s^2*omega^2+4*C*L*Rac*Rl*s^4+2*C*Kii^2*Rac*Vdc^2*s+2*C*Kii*Rac^2*Vdc*s^2+4*C*Kii*Rac*Rl*Vdc*s^2+2*C*Rac^2*Rl*s^3+2*C*Rac*Rl^2*s^3+Kpi^2*Vdc^2*s^2+2*Kpi*L*Vdc*s^3+L^2*s^4+L^2*s^2*omega^2+2*Kii*Kpi*Vdc^2*s+2*Kii*L*Vdc*s^2+2*Kpi*Rac*Vdc*s^2+2*Kpi*Rl*Vdc*s^2+2*Rac*Rl*s^2+2*L*Rac*s^3+2*L*Rl*s^3+Kii^2*Vdc^2+2*Kii*Rac*Vdc*s+2*Kii*Rl*Vdc*s+Rl^2*s^2+Rac^2*s^2);

Zoutqdcl = -omega*Rac^2*(C*Kpi^2*Vdc^2*s^2+2*C*Kpi*L*Vdc*s^3+C*L^2*s^4+C*L^2*s^2*omega^2+2*C*Kii*Kpi*Vdc^2*s+2*C*Kii*L*Vdc*s^2+2*C*Kpi*Rl*Vdc*s^2+2*C*L*Rl*s^3+C*Kii^2*Vdc^2+2*C*Kii*Rl*Vdc*s+C*Rl^2*s^2-s^2*L)/(C^2*Kpi^2*Rac^2*Vdc^2*s^4+C^2*Kpi^2*Rac^2*Vdc^2*s^2*omega^2+2*C^2*Kpi*L*Rac^2*Vdc*s^5+2*C^2*Kpi*L*Rac^2*Vdc*s^3*omega^2+C^2*L^2*Rac^2*s^6+2*C^2*L^2*Rac^2*s^4*omega^2+C^2*L^2*Rac^2*s^2*omega^4+2*C^2*Kii*Kpi*Rac^2*Vdc^2*s^3+2*C^2*Kii*Kpi*Rac^2*Vdc^2*s*omega^2+2*C^2*Kii*L*Rac^2*Vdc*s^4+2*C^2*Kii*L*Rac^2*Vdc*s^2*omega^2+2*C^2*Kpi*Rac^2*Rl*Vdc*s^4+2*C^2*Kpi*Rac^2*Rl*Vdc*s^2*omega^2+2*C^2*L*Rac^2*Rl*s^5+2*C^2*L*Rac^2*Rl*s^3*omega^2+C^2*Kii^2*Rac^2*Vdc^2*s^2+C^2*Kii^2*Rac^2*Vdc^2*omega^2+2*C^2*Kii*Rac^2*Rl*Vdc*s^3+2*C^2*Kii*Rac^2*Rl*Vdc*s*omega^2+C^2*Rac^2*Rl^2*s^4+C^2*Rac^2*Rl^2*s^2*omega^2+2*C*Kpi^2*Rac*Vdc^2*s^3+4*C*Kpi*L*Rac*Vdc*s^4+2*C*L^2*Rac*s^5+2*C*L^2*Rac*s^3*omega^2+4*C*Kii*Kpi*Rac*Vdc^2*s^2+4*C*Kii*L*Rac*Vdc*s^3+2*C*Kpi*Rac^2*Vdc*s^3+4*C*Kpi*Rac*Rl*Vdc*s^3+2*C*L*Rac^2*s^4-2*C*L*Rac^2*s^2*omega^2+4*C*L*Rac*Rl*s^4+2*C*Kii^2*Rac*Vdc^2*s+2*C*Kii*Rac^2*Vdc*s^2+4*C*Kii*Rac*Rl*Vdc*s^2+2*C*Rac^2*Rl*s^3+2*C*Rac*Rl^2*s^3+Kpi^2*Vdc^2*s^2+2*Kpi*L*Vdc*s^3+L^2*s^4+L^2*s^2*omega^2+2*Kii*Kpi*Vdc^2*s+2*Kii*L*Vdc*s^2+2*Kpi*Rac*Vdc*s^2+2*Kpi*Rl*Vdc*s^2+2*Rac*Rl*s^2+2*L*Rac*s^3+2*L*Rl*s^3+Kii^2*Vdc^2+2*Kii*Rac*Vdc*s+2*Kii*Rl*Vdc*s+Rl^2*s^2+Rac^2*s^2);

Zoutqqcl = (C*Kpi^2*Rac*Vdc^2*s^3+2*C*Kpi*L*Rac*Vdc*s^4+C*L^2*Rac*s^5+C*L^2*Rac*s^3*omega^2+2*C*Kii*Kpi*Rac*Vdc^2*s^2+2*C*Kii*L*Rac*Vdc*s^3+2*C*Kpi*Rac*Rl*Vdc*s^3+2*C*L*Rac*Rl*s^4+C*Kii^2*Rac*Vdc^2*s+2*C*Kii*Rac*Rl*Vdc*s^2+C*Rac*Rl^2*s^3+Kpi^2*Vdc^2*s^2+2*Kpi*L*Vdc*s^3+L^2*s^4+L^2*s^2*omega^2+2*Kii*Kpi*Vdc^2*s+2*Kii*L*Vdc*s^2+Kpi*Rac*Vdc*s^2+2*Kpi*Rl*Vdc*s^2+Rac*Rl*s^2+L*Rac*s^3+2*L*Rl*s^3+Kii^2*Vdc^2+Kii*Rac*Vdc*s+2*Kii*Rl*Vdc*s+Rl^2*s^2)*Rac/(C^2*Kpi^2*Rac^2*Vdc^2*s^4+C^2*Kpi^2*Rac^2*Vdc^2*s^2*omega^2+2*C^2*Kpi*L*Rac^2*Vdc*s^5+2*C^2*Kpi*L*Rac^2*Vdc*s^3*omega^2+C^2*L^2*Rac^2*s^6+2*C^2*L^2*Rac^2*s^4*omega^2+C^2*L^2*Rac^2*s^2*omega^4+2*C^2*Kii*Kpi*Rac^2*Vdc^2*s^3+2*C^2*Kii*Kpi*Rac^2*Vdc^2*s*omega^2+2*C^2*Kii*L*Rac^2*Vdc*s^4+2*C^2*Kii*L*Rac^2*Vdc*s^2*omega^2+2*C^2*Kpi*Rac^2*Rl*Vdc*s^4+2*C^2*Kpi*Rac^2*Rl*Vdc*s^2*omega^2+2*C^2*L*Rac^2*Rl*s^5+2*C^2*L*Rac^2*Rl*s^3*omega^2+C^2*Kii^2*Rac^2*Vdc^2*s^2+C^2*Kii^2*Rac^2*Vdc^2*omega^2+2*C^2*Kii*Rac^2*Rl*Vdc*s^3+2*C^2*Kii*Rac^2*Rl*Vdc*s*omega^2+C^2*Rac^2*Rl^2*s^4+C^2*Rac^2*Rl^2*s^2*omega^2+2*C*Kpi^2*Rac*Vdc^2*s^3+4*C*Kpi*L*Rac*Vdc*s^4+2*C*L^2*Rac*s^5+2*C*L^2*Rac*s^3*omega^2+4*C*Kii*Kpi*Rac*Vdc^2*s^2+4*C*Kii*L*Rac*Vdc*s^3+2*C*Kpi*Rac^2*Vdc*s^3+4*C*Kpi*Rac*Rl*Vdc*s^3+2*C*L*Rac^2*s^4-2*C*L*Rac^2*s^2*omega^2+4*C*L*Rac*Rl*s^4+2*C*Kii^2*Rac*Vdc^2*s+2*C*Kii*Rac^2*Vdc*s^2+4*C*Kii*Rac*Rl*Vdc*s^2+2*C*Rac^2*Rl*s^3+2*C*Rac*Rl^2*s^3+Kpi^2*Vdc^2*s^2+2*Kpi*L*Vdc*s^3+L^2*s^4+L^2*s^2*omega^2+2*Kii*Kpi*Vdc^2*s+2*Kii*L*Vdc*s^2+2*Kpi*Rac*Vdc*s^2+2*Kpi*Rl*Vdc*s^2+2*Rac*Rl*s^2+2*L*Rac*s^3+2*L*Rl*s^3+Kii^2*Vdc^2+2*Kii*Rac*Vdc*s+2*Kii*Rl*Vdc*s+Rl^2*s^2+Rac^2*s^2);
Rac = 15;
Zincl = -Rac^2*Vdc^2*(2*Rac*Rl*s^2+Rac^2*s^2+C^2*Kpi^2*Rac^2*Vdc^2*s^4+C^2*Kpi^2*Rac^2*Vdc^2*s^2*omega^2+2*C^2*Kpi*L*Rac^2*Vdc*s^5+2*C^2*Kpi*L*Rac^2*Vdc*s^3*omega^2+C^2*L^2*Rac^2*s^6+2*C^2*L^2*Rac^2*s^4*omega^2+C^2*L^2*Rac^2*s^2*omega^4+2*C^2*Kii*Kpi*Rac^2*Vdc^2*s^3+2*C^2*Kii*Kpi*Rac^2*Vdc^2*s*omega^2+2*C^2*Kii*L*Rac^2*Vdc*s^4+2*C^2*Kii*L*Rac^2*Vdc*s^2*omega^2+2*C^2*Kpi*Rac^2*Rl*Vdc*s^4+2*C^2*Kpi*Rac^2*Rl*Vdc*s^2*omega^2+2*C^2*L*Rac^2*Rl*s^5+2*C^2*L*Rac^2*Rl*s^3*omega^2+C^2*Kii^2*Rac^2*Vdc^2*s^2+C^2*Kii^2*Rac^2*Vdc^2*omega^2+2*C^2*Kii*Rac^2*Rl*Vdc*s^3+2*C^2*Kii*Rac^2*Rl*Vdc*s*omega^2+C^2*Rac^2*Rl^2*s^4+C^2*Rac^2*Rl^2*s^2*omega^2+2*C*Kpi^2*Rac*Vdc^2*s^3+4*C*Kpi*L*Rac*Vdc*s^4+2*C*L^2*Rac*s^5+2*C*L^2*Rac*s^3*omega^2+4*C*Kii*Kpi*Rac*Vdc^2*s^2+4*C*Kii*L*Rac*Vdc*s^3+2*C*Kpi*Rac^2*Vdc*s^3+4*C*Kpi*Rac*Rl*Vdc*s^3+2*C*L*Rac^2*s^4-2*C*L*Rac^2*s^2*omega^2+4*C*L*Rac*Rl*s^4+2*C*Kii^2*Rac*Vdc^2*s+2*C*Kii*Rac^2*Vdc*s^2+4*C*Kii*Rac*Rl*Vdc*s^2+2*C*Rac^2*Rl*s^3+2*C*Rac*Rl^2*s^3+Kpi^2*Vdc^2*s^2+2*Kpi*L*Vdc*s^3+L^2*s^4+L^2*s^2*omega^2+2*Kii*Kpi*Vdc^2*s+2*Kii*L*Vdc*s^2+2*Kpi*Rac*Vdc*s^2+2*Kpi*Rl*Vdc*s^2+2*L*Rac*s^3+2*L*Rl*s^3+Kii^2*Vdc^2+2*Kii*Rac*Vdc*s+2*Kii*Rl*Vdc*s+Rl^2*s^2)/(-2*Rl*Vod^2*C^2*L*Rac^3*s^5-Rl^2*Vod^2*C^2*L*Rac^2*s^5-3*Rl^2*Vod^2*C^2*Rac^3*s^2*omega^2-2*Rl^3*Vod^2*C^2*Rac^2*s^2*omega^2-4*Rl*Vod^2*C*L*Rac^2*s^4-2*Rl^2*Vod^2*C*L*Rac*s^4-Rl^3*omega^2*C^4*Voq^2*Rac^4*s^4-Rl^3*omega^4*C^4*Voq^2*Rac^4*s^2-Rl^2*omega^2*C^3*Voq^2*Rac^4*s^3-2*Rl^3*omega^2*C^3*Voq^2*Rac^3*s^3-2*Rl^3*omega^2*C^2*Voq^2*Rac^2*s^2+2*omega^2*L^2*C^3*Vod^2*Rac^4*s^5-omega^4*L^3*C^4*Vod^2*Rac^4*s^5+omega^4*L^2*C^3*Vod^2*Rac^4*s^3-omega^6*L^3*C^4*Vod^2*Rac^4*s^3+4*omega^2*L^2*C^2*Vod^2*Rac^3*s^4-2*omega^4*L^3*C^3*Vod^2*Rac^3*s^4+omega^2*L^2*C*Vod^2*Rac^2*s^3-2*omega^4*L^3*C^2*Vod^2*Rac^2*s^3+Kpi^2*Vdc^2*Vod^2*C^2*Rac^3*s^4+2*Kpi^2*Vdc^2*Vod^2*s^3*C*Rac^2-omega^2*L^3*Voq^2*C^2*Rac^2*s^5-2*omega^4*L^3*Voq^2*C^2*Rac^2*s^3-2*omega^2*L^3*Voq^2*C*Rac*s^4+omega^2*L^2*Voq^2*C*Rac^2*s^3-Rac^4*Vod^2*C^2*Kpi*Vdc*s^4-Rac^4*Vod^2*C^2*Kii*Vdc*s^3-Rac^4*Vod^2*C^2*Rl*s^2*omega^2-Rac^3*Vod^2*Kpi*Vdc*s^3*C-Rac^3*Vod^2*C*Kii*Vdc*s^2+Rac^4*Vod^2*C^2*s^3*omega^2*L+2*Rac^3*Vod^2*s^2*omega^2*L*C-2*s^4*omega^4*L^3*C^3*Voq^2*Rac^3-s^2*omega^4*L^2*Vod^2*C^2*Rac^3+4*s^4*omega^2*Rac^3*Voq^2*C^2*L^2-s^2*omega^4*Rac^3*Voq^2*C^2*L^2+Kii^2*Vdc^2*Vod^2*C^2*Rac^3*s^2+Kii^2*Vdc^2*Vod^2*C^2*Rac^3*omega^2+2*Kii^2*Vdc^2*Vod^2*C*Rac^2*s+2*s^2*omega^2*Rac^3*Voq^2*C*L-s^2*omega^2*Rac^4*Voq^2*C^2*Rl+Kii*Vdc*Vod^2*s^2*L*Rac+Kii*Vdc*Vod^2*s^2*L*Rl+2*s*Kpi*Vdc^2*Vod^2*Kii*Rac+2*s*Kpi*Vdc^2*Vod^2*Kii*Rl+Kpi*Vdc*Vod^2*s^3*L*Rac+Kpi*Vdc*Vod^2*s^3*L*Rl+Kpi^2*Vdc^2*Voq^2*C^2*Rac^3*s^4+2*Kpi^2*Vdc^2*Voq^2*s^3*C*Rac^2-Rl^2*Voq^2*C^2*L*Rac^2*s^5-2*Rl*Voq^2*C^2*L*Rac^3*s^5-3*Rl^2*Voq^2*C^2*Rac^3*s^2*omega^2-2*Rl^2*Voq^2*C*L*Rac*s^4-4*Rl*Voq^2*C*L*Rac^2*s^4-Rac^4*Voq^2*C^2*Kpi*Vdc*s^4-Rac^4*Voq^2*C^2*Kii*Vdc*s^3-Rac^3*Voq^2*Kpi*Vdc*s^3*C-Rac^3*Voq^2*C*Kii*Vdc*s^2-omega^2*L^3*Vod^2*C^2*Rac^2*s^5-2*omega^2*L^3*Vod^2*C*Rac*s^4+2*omega^2*L^2*C^3*Voq^2*Rac^4*s^5-omega^4*L^3*C^4*Voq^2*Rac^4*s^5+omega^4*L^2*C^3*Voq^2*Rac^4*s^3-omega^6*L^3*C^4*Voq^2*Rac^4*s^3-Rl^3*omega^2*C^4*Vod^2*Rac^4*s^4-Rl^3*omega^4*C^4*Vod^2*Rac^4*s^2+2*Kii^2*Vdc^2*Voq^2*C*Rac*s*Rl+Kii^2*Vdc^2*omega^4*C^4*Vod^2*Rac^4*Rl-Rl^2*omega^2*C^4*Vod^2*Rac^4*L*s^5-Rl^2*omega^4*C^4*Vod^2*Rac^4*L*s^3+Kpi^2*Vdc^2*Voq^2*C^2*Rac^3*s^2*omega^2+Kpi*Vdc*Voq^2*C^2*L*Rac^3*s^5+2*Kpi*Vdc^2*Voq^2*C^2*Kii*Rac^3*s^3+2*Kpi*Vdc*Voq^2*C*L*Rac^2*s^4+4*Kpi*Vdc^2*Voq^2*C*Kii*Rac^2*s^2+Kii*Vdc*Voq^2*C^2*L*Rac^3*s^4+2*Kii*Vdc*Voq^2*C*L*Rac^2*s^3+Rac^4*Voq^2*C^2*s^3*omega^2*L-omega^2*L^3*Vod^2*s^3+2*omega^2*L*Voq^2*C^2*Kii*Rac^2*Vdc*s^2*Rl+2*omega^2*L*C^2*Vod^2*Rac^2*Kpi*Vdc*s^3*Rl+2*omega^2*L*C^2*Vod^2*Rac^2*Kii*Vdc*s^2*Rl+Kpi*Vdc*omega^2*C^4*Voq^2*Rac^4*L*s^5*Rl+Kpi*Vdc*omega^4*C^4*Voq^2*Rac^4*L*s^3*Rl+2*Kpi*Vdc^2*omega^2*C^4*Voq^2*Rac^4*Kii*s^3*Rl+4*Kpi*Vdc^2*omega^2*C^3*Voq^2*Rac^3*Kii*s^2*Rl+Rac^4*Voq^2*C^3*Kii*Vdc*s^3*omega^2*L+Rac^4*Voq^2*C^3*Kpi*Vdc*s^4*omega^2*L-Rl^3*Voq^2*s^2+Kii*Vdc*omega^2*C^4*Voq^2*Rac^4*L*s^4*Rl+Kii*Vdc*omega^4*C^4*Voq^2*Rac^4*L*s^2*Rl+2*Kpi^2*Vdc^2*Vod^2*C^2*Rac^2*s^2*omega^2*Rl+Kpi*Vdc*Vod^2*C^2*L*Rac^2*s^5*Rl+Rac^4*Vod^2*C^3*Kpi*Vdc*s^4*omega^2*L+Rac^4*Vod^2*C^3*Kii*Vdc*s^3*omega^2*L+3*Rac^3*Vod^2*Kpi*Vdc*s^3*C^2*omega^2*L+3*Rac^3*Vod^2*C^2*Kii*Vdc*s^2*omega^2*L-Rac^3*Vod^2*s^2-Rl^3*Vod^2*s^2-Rac^3*Voq^2*s^2+2*Kpi*Vdc^2*Vod^2*C^2*Kii*Rac^2*s^3*Rl+2*Kpi*Vdc*Vod^2*C*L*Rac*s^4*Rl+4*Kpi*Vdc^2*Vod^2*C*Kii*Rac*s^2*Rl+Kpi^2*Vdc^2*omega^2*C^4*Voq^2*Rac^4*s^4*Rl+Kpi^2*Vdc^2*omega^4*C^4*Voq^2*Rac^4*s^2*Rl+2*Kpi^2*Vdc^2*omega^2*C^3*Voq^2*Rac^3*s^3*Rl+3*Kpi*Vdc*omega^2*C^2*Voq^2*Rac^3*s^3*L+Kpi*Vdc*omega^2*C^3*Voq^2*Rac^4*s^3*Rl+2*Kpi^2*Vdc^2*omega^2*C^2*Voq^2*Rac^2*s^2*Rl+Kii*Vdc*Vod^2*C^2*L*Rac^2*s^4*Rl+2*Kii*Vdc*Vod^2*C*L*Rac*s^3*Rl+3*Kii*Vdc*omega^2*C^2*Voq^2*Rac^3*s^2*L+Kii*Vdc*omega^2*C^3*Voq^2*Rac^4*s^2*Rl+Kii^2*Vdc^2*omega^2*C^4*Voq^2*Rac^4*s^2*Rl+2*Kii^2*Vdc^2*omega^2*C^3*Voq^2*Rac^3*s*Rl+2*omega^2*L*Voq^2*C^3*Kpi*Rac^3*Vdc*s^4*Rl+2*s*Kpi*Vdc^2*Vod^2*C^2*Kii*Rac^3*omega^2+4*s*Kpi*Vdc^2*Vod^2*C^2*Kii*Rac^2*omega^2*Rl+2*s*Kpi*Vdc^2*omega^4*C^4*Voq^2*Rac^4*Kii*Rl+4*s*Kpi*Vdc^2*omega^2*C^2*Voq^2*Rac^2*Kii*Rl+2*omega^2*L*Voq^2*C^3*Kii*Rac^3*Vdc*s^3*Rl+2*omega^2*L*Voq^2*Kpi*Vdc*s^3*C^2*Rac^2*Rl+2*Rac^4*Vod^2*C^3*Rl*s^2*omega^4*L-Rac^2*Vod^2*Kpi*Vdc*s^3*C*Rl-Rac^2*Vod^2*C*Kii*Vdc*s^2*Rl-2*omega^2*L^2*Voq^2*s^3*Rl*C*Rac-s^4*omega^4*L^2*C^4*Voq^2*Rac^4*Rl-2*s^2*omega^4*L^2*Vod^2*C^2*Rac^2*Rl-s^2*omega^6*L^2*C^4*Voq^2*Rac^4*Rl-2*s^2*omega^4*L^2*C^2*Voq^2*Rac^2*Rl+2*s^4*omega^2*Rac^4*Voq^2*C^3*L*Rl+2*s^2*omega^4*Rac^4*Voq^2*C^3*L*Rl+Kii^2*Vdc^2*Vod^2*C^2*Rac^2*s^2*Rl+2*Kii^2*Vdc^2*Vod^2*C^2*Rac^2*omega^2*Rl+2*Kii^2*Vdc^2*Vod^2*C*Rac*s*Rl+Kii^2*Vdc^2*omega^4*C^4*Voq^2*Rac^4*Rl+2*Kii^2*Vdc^2*omega^2*C^2*Voq^2*Rac^2*Rl+2*Rac^2*Vod^2*s^2*Rl*omega^2*L*C+2*Rl*Vod^2*C^2*L*Rac^3*s^3*omega^2-2*Rl^2*Vod^2*C^2*L*Rac^2*s^3*omega^2-Rl^2*omega^2*C^4*Voq^2*Rac^4*L*s^5-Rl^2*omega^4*C^4*Voq^2*Rac^4*L*s^3+2*omega^2*L*C^3*Vod^2*Rac^4*Rl*s^4-2*omega^2*L*C^3*Vod^2*Rac^3*Rl^2*s^4-omega^4*L^2*C^4*Vod^2*Rac^4*Rl*s^4-omega^6*L^2*C^4*Vod^2*Rac^4*Rl*s^2+Kpi^2*Vdc^2*Vod^2*C^2*Rac^3*s^2*omega^2+Kpi*Vdc*Vod^2*C^2*L*Rac^3*s^5+2*Kpi*Vdc^2*Vod^2*C^2*Kii*Rac^3*s^3+2*Kpi*Vdc*Vod^2*C*L*Rac^2*s^4+4*Kpi*Vdc^2*Vod^2*C*Kii*Rac^2*s^2+Kii*Vdc*Vod^2*C^2*L*Rac^3*s^4+2*Kii*Vdc*Vod^2*C*L*Rac^2*s^3+2*omega^2*L*Voq^2*s^2*Rac^2*Rl*C+Kpi^2*Vdc^2*Voq^2*C^2*Rac^2*s^4*Rl+2*Kpi^2*Vdc^2*Voq^2*s^3*C*Rac*Rl-Rac^3*Voq^2*C^2*Kpi*Vdc*s^4*Rl-Rac^3*Voq^2*C^2*Kii*Vdc*s^3*Rl-Rac^2*Voq^2*Kpi*Vdc*s^3*C*Rl-Rac^2*Voq^2*C*Kii*Vdc*s^2*Rl+Kii^2*Vdc^2*Voq^2*C^2*Rac^2*s^2*Rl+Kpi^2*Vdc^2*Voq^2*s^2*Rl+Kpi^2*Vdc^2*Voq^2*s^2*Rac-Rl*Vod^2*C^2*L^2*Rac^2*s^4*omega^2-2*Rl^2*omega^2*C^3*Voq^2*Rac^3*s^4*L+2*Rl*omega^2*C^2*Voq^2*Rac^3*s^3*L-2*Rl^2*omega^2*C^2*Voq^2*Rac^2*s^3*L-2*omega^4*L^2*C^3*Vod^2*Rac^3*s^3*Rl-2*omega^2*L^2*C*Vod^2*Rac*s^3*Rl+Kpi^2*Vdc^2*Vod^2*C^2*Rac^2*s^4*Rl+2*Kpi^2*Vdc^2*Vod^2*s^3*C*Rac*Rl-2*omega^4*L^2*Voq^2*C^3*Rac^3*s^3*Rl-omega^2*L^2*Voq^2*C^2*Rac^2*s^4*Rl-Rac^3*Vod^2*C^2*Kpi*Vdc*s^4*Rl-Rac^3*Vod^2*C^2*Kii*Vdc*s^3*Rl-2*Rac*Vod^2*s^3*L*Rl-4*Rac^3*Vod^2*C*s^3*Rl-s^2*omega^2*L^2*Vod^2*Rac-s^2*omega^2*L^2*Vod^2*Rl+Kpi^2*Vdc^2*Vod^2*s^2*Rac+Kpi^2*Vdc^2*Vod^2*s^2*Rl-omega^2*L^2*Voq^2*s^2*Rac-omega^2*L^2*Voq^2*s^2*Rl-Rac^4*Vod^2*C^2*L*s^5-Rac^4*Vod^2*C^2*Rl*s^4-2*Rac^3*Vod^2*C^2*Rl^2*s^4-2*Rac^3*Vod^2*C*L*s^4-5*Rac^2*Vod^2*s^3*C*Rl^2-Rl^3*Vod^2*C^2*Rac^2*s^4-2*Rl^3*Vod^2*s^3*C*Rac-Rl^3*Voq^2*C^2*Rac^2*s^4-2*Rl^2*Voq^2*C^2*Rac^3*s^4-5*Rl^2*Voq^2*C*Rac^2*s^3-4*Rl*Voq^2*C*Rac^3*s^3-2*Rl^3*Voq^2*s^3*C*Rac-Rac^4*Voq^2*C^2*L*s^5-Rac^4*Voq^2*C^2*Rl*s^4-2*Rac^3*Voq^2*C*L*s^4-2*Rl*Voq^2*s^3*L*Rac-Rac^2*Vod^2*s^3*L-Rac^4*Vod^2*C*s^3-omega^2*L^3*Voq^2*s^3-3*Rac^2*Vod^2*s^2*Rl-3*Rl^2*Vod^2*s^2*Rac+Kii^2*Vdc^2*Vod^2*Rac+Kii^2*Vdc^2*Vod^2*Rl-Rl^2*Vod^2*s^3*L-3*Rl^2*Voq^2*s^2*Rac-3*Rac^2*Voq^2*s^2*Rl+Kii^2*Vdc^2*Voq^2*Rl+Kii^2*Vdc^2*Voq^2*Rac-Rl^2*Voq^2*s^3*L-Rac^4*Voq^2*C*s^3-Rac^2*Voq^2*s^3*L+2*s^4*omega^2*Kpi*Vdc*Vod^2*C^3*L*Rac^3*Rl+2*s*Kpi*Vdc^2*Voq^2*C^2*Kii*Rac^3*omega^2+2*s*Kpi*Vdc^2*omega^4*C^4*Vod^2*Rac^4*Kii*Rl+2*omega^2*L*Vod^2*C^3*Kii*Rac^3*Vdc*s^3*Rl+Kii*Vdc*Voq^2*C^2*L*Rac^2*s^4*Rl+2*Kii*Vdc*Voq^2*C*L*Rac*s^3*Rl+Kii*Vdc*omega^2*C^3*Vod^2*Rac^4*s^2*Rl+Kpi*Vdc*omega^2*C^4*Vod^2*Rac^4*L*s^5*Rl+Kpi*Vdc*omega^4*C^4*Vod^2*Rac^4*L*s^3*Rl+2*Kpi*Vdc^2*omega^2*C^4*Vod^2*Rac^4*Kii*s^3*Rl+4*Kpi*Vdc^2*omega^2*C^3*Vod^2*Rac^3*Kii*s^2*Rl+Kii*Vdc*omega^2*C^4*Vod^2*Rac^4*L*s^4*Rl+Kii*Vdc*omega^4*C^4*Vod^2*Rac^4*L*s^2*Rl+Kpi*Vdc*Voq^2*C^2*L*Rac^2*s^5*Rl+2*Kpi*Vdc^2*Voq^2*C^2*Kii*Rac^2*s^3*Rl+2*Kpi*Vdc*Voq^2*C*L*Rac*s^4*Rl+4*Kpi*Vdc^2*Voq^2*C*Kii*Rac*s^2*Rl+Kpi^2*Vdc^2*omega^2*C^4*Vod^2*Rac^4*s^4*Rl+Kpi^2*Vdc^2*omega^4*C^4*Vod^2*Rac^4*s^2*Rl+2*Kpi^2*Vdc^2*omega^2*C^3*Vod^2*Rac^3*s^3*Rl+Kpi*Vdc*omega^2*C^3*Vod^2*Rac^4*s^3*Rl+Kii^2*Vdc^2*omega^2*C^4*Vod^2*Rac^4*s^2*Rl+2*Kii^2*Vdc^2*omega^2*C^3*Vod^2*Rac^3*s*Rl-Rl^2*omega^2*C^3*Vod^2*Rac^4*s^3-2*Rl^3*omega^2*C^3*Vod^2*Rac^3*s^3+Kii^2*Vdc^2*Voq^2*C^2*Rac^3*s^2+Kii^2*Vdc^2*Voq^2*C^2*Rac^3*omega^2+2*Kii^2*Vdc^2*Voq^2*C*Rac^2*s+Kpi*Vdc*Voq^2*s^3*L*Rl+Kpi*Vdc*Voq^2*s^3*L*Rac+2*s*Kpi*Vdc^2*Voq^2*Kii*Rl+2*s*Kpi*Vdc^2*Voq^2*Kii*Rac+Kii*Vdc*Voq^2*s^2*L*Rl+Kii*Vdc*Voq^2*s^2*L*Rac);
% Rac^2*Vdc^2*(C^2*Kpi^2*Rac^2*Vdc^2*s^4+C^2*Kpi^2*Rac^2*Vdc^2*s^2*omega^2+2*C^2*Kpi*L*Rac^2*Vdc*s^5+2*C^2*Kpi*L*Rac^2*Vdc*s^3*omega^2+C^2*L^2*Rac^2*s^6+2*C^2*L^2*Rac^2*s^4*omega^2+C^2*L^2*Rac^2*s^2*omega^4+2*C^2*Kii*Kpi*Rac^2*Vdc^2*s^3+2*C^2*Kii*Kpi*Rac^2*Vdc^2*s*omega^2+2*C^2*Kii*L*Rac^2*Vdc*s^4+2*C^2*Kii*L*Rac^2*Vdc*s^2*omega^2+2*C^2*Kpi*Rac^2*Rl*Vdc*s^4+2*C^2*Kpi*Rac^2*Rl*Vdc*s^2*omega^2+2*C^2*L*Rac^2*Rl*s^5+2*C^2*L*Rac^2*Rl*s^3*omega^2+C^2*Kii^2*Rac^2*Vdc^2*s^2+C^2*Kii^2*Rac^2*Vdc^2*omega^2+2*C^2*Kii*Rac^2*Rl*Vdc*s^3+2*C^2*Kii*Rac^2*Rl*Vdc*s*omega^2+C^2*Rac^2*Rl^2*s^4+C^2*Rac^2*Rl^2*s^2*omega^2+2*C*Kpi^2*Rac*Vdc^2*s^3+4*C*Kpi*L*Rac*Vdc*s^4+2*C*L^2*Rac*s^5+2*C*L^2*Rac*s^3*omega^2+4*C*Kii*Kpi*Rac*Vdc^2*s^2+4*C*Kii*L*Rac*Vdc*s^3+2*C*Kpi*Rac^2*Vdc*s^3+4*C*Kpi*Rac*Rl*Vdc*s^3+2*C*L*Rac^2*s^4-2*C*L*Rac^2*s^2*omega^2+4*C*L*Rac*Rl*s^4+2*C*Kii^2*Rac*Vdc^2*s+2*C*Kii*Rac^2*Vdc*s^2+4*C*Kii*Rac*Rl*Vdc*s^2+2*C*Rac^2*Rl*s^3+2*C*Rac*Rl^2*s^3+Kpi^2*Vdc^2*s^2+2*Kpi*L*Vdc*s^3+L^2*s^4+L^2*s^2*omega^2+2*Kii*Kpi*Vdc^2*s+2*Kii*L*Vdc*s^2+2*Kpi*Rac*Vdc*s^2+2*Kpi*Rl*Vdc*s^2+2*Rac*Rl*s^2+2*L*Rac*s^3+2*L*Rl*s^3+Kii^2*Vdc^2+2*Kii*Rac*Vdc*s+2*Kii*Rl*Vdc*s+Rl^2*s^2+Rac^2*s^2)/(s*(Rac^3*Vod^2*s+Rl^3*Vod^2*s+Rl^3*Voq^2*s+Rac^3*Voq^2*s-2*Rac^3*Vod^2*s*omega^2*L*C+2*Rac^3*Vod^2*C^2*L*s^4*Rl-2*Rac^4*Vod^2*C^3*L^2*s^4*omega^2+4*Rac^2*Vod^2*C*L*s^3*Rl-4*Rac^3*Vod^2*C^2*L^2*s^3*omega^2+2*Rac*Vod^2*s*Vdc*Kpi*Rl-2*omega^2*s*Rac^3*Voq^2*C*L+omega^2*s*Rac^4*Voq^2*C^2*Rl-Rac^2*Vod^2*s^2*L^2*omega^2*C-Rac^4*Vod^2*C^2*s^2*omega^2*L+omega^2*L^3*Voq^2*C^2*Rac^2*s^4+2*omega^4*L^3*Voq^2*C^2*Rac^2*s^2+2*omega^2*L^3*Voq^2*C*Rac*s^3-omega^2*L^2*Voq^2*C*Rac^2*s^2+omega^2*L^2*Voq^2*s*Vdc*Kpi+Rac^4*Vod^2*C^2*Kpi*Vdc*s^3+Rac^4*Vod^2*C^2*Kii*Vdc*s^2+Rac^4*Vod^2*C^2*Kii*Vdc*omega^2+Rac^4*Vod^2*C^2*Rl*s*omega^2+3*Rac^3*Vod^2*C^2*Rl^2*s*omega^2+2*Rac^3*Vod^2*Kpi*Vdc*s^2*C+2*Rac^3*Vod^2*C*Kii*Vdc*s+Rl^2*Vod^2*C^2*L*Rac^2*s^4+2*Rl^3*Vod^2*C^2*Rac^2*s*omega^2+2*Rl^2*Vod^2*C*L*Rac*s^3+Rl^3*omega^2*C^4*Voq^2*Rac^4*s^3+Rl^3*omega^4*C^4*Voq^2*Rac^4*s+Rl^2*omega^2*C^3*Voq^2*Rac^4*s^2+2*Rl^3*omega^2*C^3*Voq^2*Rac^3*s^2+2*Rl^3*omega^2*C^2*Voq^2*Rac^2*s+omega^4*L^3*C^4*Vod^2*Rac^4*s^4-omega^4*L^2*C^3*Vod^2*Rac^4*s^2+omega^6*L^3*C^4*Vod^2*Rac^4*s^2+2*omega^4*L^3*C^3*Vod^2*Rac^3*s^3+2*omega^4*L^3*C^2*Vod^2*Rac^2*s^2-4*omega^2*s^3*Rac^3*Voq^2*C^2*L^2+omega^4*s*Rac^3*Voq^2*C^2*L^2+omega^4*s*L^2*Vod^2*C^2*Rac^3+2*omega^4*s^3*L^3*C^3*Voq^2*Rac^3+3*omega^2*s*Rl^2*C^2*Voq^2*Rac^3+Rl^2*Voq^2*C^2*L*Rac^2*s^4+2*Rl*Voq^2*C^2*L*Rac^3*s^4+2*Rl^2*Voq^2*C*L*Rac*s^3+4*Rl*Voq^2*C*L*Rac^2*s^3+Rac^4*Voq^2*C^2*Kpi*Vdc*s^3+Rac^4*Voq^2*C^2*Kii*Vdc*s^2+Rac^4*Voq^2*C^2*Kii*Vdc*omega^2+2*Rac^3*Voq^2*Kpi*Vdc*s^2*C+2*Rac^3*Voq^2*C*Kii*Vdc*s+omega^2*L^3*Vod^2*C^2*Rac^2*s^4+2*omega^2*L^3*Vod^2*C*Rac*s^3+omega^2*L^2*Vod^2*s*Vdc*Kpi-2*omega^2*L^2*C^3*Voq^2*Rac^4*s^4+omega^4*L^3*C^4*Voq^2*Rac^4*s^4-omega^4*L^2*C^3*Voq^2*Rac^4*s^2+omega^6*L^3*C^4*Voq^2*Rac^4*s^2+Rl^3*omega^2*C^4*Vod^2*Rac^4*s^3+Rl^3*omega^4*C^4*Vod^2*Rac^4*s+Rl^2*omega^2*C^3*Vod^2*Rac^4*s^2+2*Rl^3*omega^2*C^3*Vod^2*Rac^3*s^2+2*Rl*Voq^2*s*Vdc*Kpi*Rac-Rac^4*Voq^2*C^2*s^2*omega^2*L+omega^2*L^2*Voq^2*C^2*Kpi*Rac^2*Vdc*s^3+2*omega^4*L^2*Voq^2*C^2*Kpi*Rac^2*Vdc*s+omega^2*L^2*Voq^2*C^2*Kii*Rac^2*Vdc*s^2+2*omega^2*L^2*Voq^2*Kpi*Vdc*s^2*C*Rac+2*omega^2*L^2*Voq^2*C*Kii*Rac*Vdc*s+2*Rl*Vod^2*C^2*Kpi*Rac^3*Vdc*s*omega^2+2*Rl^2*Vod^2*C^2*Kpi*Rac^2*Vdc*s*omega^2+Rl^2*omega^2*C^4*Voq^2*Rac^4*Kpi*Vdc*s^3+Rl^2*omega^4*C^4*Voq^2*Rac^4*Kpi*Vdc*s+3*Rac^2*Vod^2*s*Rl+3*Rl^2*Vod^2*s*Rac+Rl^2*Vod^2*Kii*Vdc+Rl^2*Vod^2*s^2*L+Rac^2*Vod^2*Kii*Vdc+Rac^2*Vod^2*s^2*L+Rac^4*Vod^2*C*s^2+omega^2*L^3*Voq^2*s^2+3*Rl^2*Voq^2*s*Rac+3*Rac^2*Voq^2*s*Rl+omega^2*L^3*Vod^2*s^2+Rac^2*Voq^2*Kii*Vdc+Rac^2*Voq^2*s^2*L+Rac^4*Voq^2*C*s^2+Rl^2*Voq^2*Kii*Vdc+Rl^2*Voq^2*s^2*L-4*Rac^3*Voq^2*Kpi*Vdc*s^2*C^2*omega^2*L-4*Rac^3*Voq^2*C^2*Kii*Vdc*s*omega^2*L+2*Rl*Voq^2*C^2*Kpi*Rac^3*Vdc*s*omega^2+omega^2*L^2*Vod^2*C^2*Kpi*Rac^2*Vdc*s^3+omega^2*L^2*Vod^2*C^2*Kii*Rac^2*Vdc*s^2+2*omega^2*L^2*Vod^2*Kpi*Vdc*s^2*C*Rac+2*omega^2*L^2*Vod^2*C*Kii*Rac*Vdc*s+omega^4*L^2*C^4*Voq^2*Rac^4*Kpi*Vdc*s^3-2*omega^4*L*C^3*Voq^2*Rac^4*Kpi*Vdc*s+omega^6*L^2*C^4*Voq^2*Rac^4*Kpi*Vdc*s+omega^4*L^2*C^4*Voq^2*Rac^4*Kii*Vdc*s^2+2*omega^4*L^2*C^3*Voq^2*Rac^3*Kpi*Vdc*s^2+2*omega^4*L^2*C^3*Voq^2*Rac^3*Kii*Vdc*s-2*omega^2*L*C*Voq^2*Rac^2*s*Vdc*Kpi+Rl^2*omega^2*C^4*Vod^2*Rac^4*Kpi*Vdc*s^3+Rl^2*omega^4*C^4*Vod^2*Rac^4*Kpi*Vdc*s+Rl^2*omega^2*C^4*Vod^2*Rac^4*Kii*Vdc*s^2+2*Rl^2*omega^2*C^3*Vod^2*Rac^3*Kpi*Vdc*s^2+2*Rl^2*omega^2*C^3*Vod^2*Rac^3*Kii*Vdc*s-2*Rac^4*Voq^2*C^3*Kpi*Vdc*s^3*omega^2*L-2*Rac^4*Voq^2*C^3*Kii*Vdc*s^2*omega^2*L+Rl^2*omega^2*C^4*Voq^2*Rac^4*Kii*Vdc*s^2+2*Rl^2*omega^2*C^3*Voq^2*Rac^3*Kpi*Vdc*s^2+2*Rl^2*omega^2*C^3*Voq^2*Rac^3*Kii*Vdc*s+2*Rl^2*omega^2*C^2*Voq^2*Rac^2*s*Vdc*Kpi-2*omega^2*L*C^3*Vod^2*Rac^4*Kpi*Vdc*s^3+omega^4*L^2*C^4*Vod^2*Rac^4*Kpi*Vdc*s^3-2*omega^4*L*C^3*Vod^2*Rac^4*Kpi*Vdc*s+omega^6*L^2*C^4*Vod^2*Rac^4*Kpi*Vdc*s-2*omega^2*L*C^3*Vod^2*Rac^4*Kii*Vdc*s^2+omega^4*L^2*C^4*Vod^2*Rac^4*Kii*Vdc*s^2-4*omega^2*L*C^2*Vod^2*Rac^3*Kpi*Vdc*s^2+2*omega^4*L^2*C^3*Vod^2*Rac^3*Kpi*Vdc*s^2-4*omega^2*L*C^2*Vod^2*Rac^3*Kii*Vdc*s+2*omega^4*L^2*C^3*Vod^2*Rac^3*Kii*Vdc*s-2*omega^2*L*C*Vod^2*Rac^2*s*Vdc*Kpi+2*omega^4*L^2*C^2*Vod^2*Rac^2*s*Vdc*Kpi+2*omega^4*L^2*Voq^2*C^2*Kii*Rac^2*Vdc+omega^2*L^2*Voq^2*C^2*Rac^2*Rl*s^3+2*omega^2*L*Voq^2*C^3*Rac^3*Rl^2*s^3+2*omega^4*L^2*Voq^2*C^2*Rac^2*Rl*s+2*omega^2*L^2*Voq^2*s^2*C*Rac*Rl+2*omega^2*L*Voq^2*s^2*C^2*Rac^2*Rl^2+Rac^4*Vod^2*C^2*Kpi*Vdc*s*omega^2+2*Rl*Vod^2*C^2*Kpi*Rac^3*Vdc*s^3+Rl^2*Vod^2*C^2*Kpi*Rac^2*Vdc*s^3+2*Rl*Vod^2*C^3*L^2*Rac^3*s^2*omega^4+2*Rl*Vod^2*C^2*Kii*Rac^3*Vdc*s^2+Rl^2*Vod^2*C^2*Kii*Rac^2*Vdc*s^2+2*Rl*Vod^2*C^2*Kii*Rac^3*Vdc*omega^2+2*Rl^2*Vod^2*C^2*Kii*Rac^2*Vdc*omega^2+4*Rl*Vod^2*Kpi*Vdc*s^2*C*Rac^2+2*Rl^2*Vod^2*Kpi*Vdc*s^2*C*Rac+4*Rl*Vod^2*C*Kii*Rac^2*Vdc*s+2*Rl^2*Vod^2*C*Kii*Rac*Vdc*s+Rl^2*omega^2*C^4*Voq^2*Rac^4*L*s^4+2*Rl*omega^4*C^3*Voq^2*Rac^3*L^2*s^2+Rl^2*omega^4*C^4*Voq^2*Rac^4*L*s^2+Rl^2*omega^4*C^4*Voq^2*Rac^4*Kii*Vdc+2*Rl^2*omega^2*C^2*Voq^2*Rac^2*Kii*Vdc-2*omega^4*L*C^3*Vod^2*Rac^4*Kii*Vdc+omega^6*L^2*C^4*Vod^2*Rac^4*Kii*Vdc+omega^4*L^2*C^4*Vod^2*Rac^4*Rl*s^3-2*omega^4*L*C^3*Vod^2*Rac^4*Rl*s+omega^6*L^2*C^4*Vod^2*Rac^4*Rl*s-2*omega^2*L*C*Vod^2*Rac^2*Kii*Vdc+2*omega^4*L^2*C^2*Vod^2*Rac^2*Kii*Vdc-2*omega^2*L*C*Vod^2*Rac^2*s*Rl+2*omega^4*L^2*C^2*Vod^2*Rac^2*s*Rl+2*Rl*Vod^2*s^2*L^2*omega^2*C*Rac-2*omega^2*L*Voq^2*C^2*Rac^3*s^2*Rl+Rl*Vod^2*C^2*L^2*Rac^2*s^3*omega^2-2*omega^2*s^3*Rac^4*Voq^2*C^3*L*Rl-2*omega^4*s*Rac^4*Voq^2*C^3*L*Rl+omega^6*s*L^2*C^4*Voq^2*Rac^4*Rl+omega^4*s^3*L^2*C^4*Voq^2*Rac^4*Rl+2*Rac^3*Voq^2*C^2*Kpi*Vdc*s^3*Rl+2*Rac^3*Voq^2*C^2*Kii*Vdc*s^2*Rl+2*Rac^3*Voq^2*C^2*Kii*Vdc*omega^2*Rl-2*Rac^4*Voq^2*C^3*Kii*Vdc*omega^4*L+4*Rac^2*Voq^2*Kpi*Vdc*s^2*C*Rl+4*Rac^2*Voq^2*C*Kii*Vdc*s*Rl+Rl^2*Voq^2*C^2*Kpi*Rac^2*Vdc*s^3+Rl^2*Voq^2*C^2*Kii*Rac^2*Vdc*s^2+2*Rl^2*Voq^2*Kpi*Vdc*s^2*C*Rac+2*Rl^2*Voq^2*C*Kii*Rac*Vdc*s+Rac^4*Voq^2*C^2*Kpi*Vdc*s*omega^2+omega^6*L^2*C^4*Voq^2*Rac^4*Kii*Vdc-2*omega^2*L*C*Voq^2*Rac^2*Kii*Vdc+Rl^2*omega^2*C^4*Vod^2*Rac^4*L*s^4+Rl^2*omega^4*C^4*Vod^2*Rac^4*L*s^2+Rl^2*omega^4*C^4*Vod^2*Rac^4*Kii*Vdc+omega^2*L^2*Voq^2*Kii*Vdc+omega^2*L^2*Voq^2*s*Rac+omega^2*L^2*Voq^2*s*Rl+Rac^4*Vod^2*C^2*L*s^4+Rac^4*Vod^2*C^2*Rl*s^3+2*Rac^3*Vod^2*C^2*Rl^2*s^3+2*Rac^3*Vod^2*C*L*s^3+4*Rac^3*Vod^2*s^2*C*Rl+5*Rac^2*Vod^2*s^2*C*Rl^2+Rac^2*Vod^2*s*Vdc*Kpi+Rl^3*Vod^2*C^2*Rac^2*s^3+2*Rl^3*Vod^2*s^2*C*Rac+Rl^2*Vod^2*s*Vdc*Kpi+omega^2*s*L^2*Vod^2*Rac+omega^2*s*L^2*Vod^2*Rl+2*Rl*Vod^2*Kii*Vdc*Rac+2*Rl*Vod^2*s^2*L*Rac+Rl^3*Voq^2*C^2*Rac^2*s^3+2*Rl^2*Voq^2*C^2*Rac^3*s^3+5*Rl^2*Voq^2*C*Rac^2*s^2+4*Rl*Voq^2*C*Rac^3*s^2+2*Rl^3*Voq^2*s^2*C*Rac+Rl^2*Voq^2*s*Vdc*Kpi+Rac^4*Voq^2*C^2*L*s^4+Rac^4*Voq^2*C^2*Rl*s^3+2*Rac^3*Voq^2*C*L*s^3+Rac^2*Voq^2*s*Vdc*Kpi+omega^2*L^2*Vod^2*Kii*Vdc+2*Rac*Voq^2*Kii*Vdc*Rl+2*Rac*Voq^2*s^2*L*Rl+2*Rl^2*Vod^2*s^2*C^2*Rac^2*omega^2*L+2*Rl^2*Vod^2*C^3*Rac^3*s^3*omega^2*L-2*omega^2*L*Voq^2*s*Rac^2*Rl*C-2*Rac^4*Vod^2*C^3*Rl*s^3*omega^2*L-2*Rac^3*Vod^2*s^2*C^2*Rl*omega^2*L));

ZinwCdccl = 1/(1/Zincl+Cdc*s);

% load TFidqidqref_VSI.mat
load TFidqidqref_VSI_no_decoup.mat
% load TFidqidqref_VSI_no_decoup_delay.mat
tfididref_matlab = tfidqidqref_VSI(1,1);
tfidiqref_matlab = tfidqidqref_VSI(1,2);
tfiqidref_matlab = tfidqidqref_VSI(2,1);
tfiqiqref_matlab = tfidqidqref_VSI(2,2);

% load TFvdqidqref_VSI.mat
% load TFvdqidqref_VSI_no_decoup.mat
load TFvdqidqref_VSI_no_decoup_delay.mat
tfvdidref_matlab = tfvdqidqref_VSI(1,1);
tfvdiqref_matlab = tfvdqidqref_VSI(1,2);
tfvqidref_matlab = tfvdqidqref_VSI(2,1);
tfvqiqref_matlab = tfvdqidqref_VSI(2,2);

% load ZoutVSI_i_loop.mat
% load ZoutVSI_i_loop_no_decoup.mat
load ZoutVSI_i_loop_no_decoup_delay.mat
Zoutddicl_matlab = ZoutVSI(1,1);
Zoutdqicl_matlab = ZoutVSI(1,2);
Zoutqdicl_matlab = ZoutVSI(2,1);
Zoutqqicl_matlab = ZoutVSI(2,2);

% load ZinVSI_i_loop.mat
% load ZinVSI_i_loop_no_decoup.mat
load ZinVSI_i_loop_no_decoup_delay.mat
Zinicl_matlab = ZinVSI;
ZinwCdccl_matlab = 1/(1/Zinicl_matlab+Cdc*s);

% figure(1);
% bode(tfididref,P);
% hold on
% bode(tfididref_matlab,P);
% Bode_Darklines(3);
% hold off
% 
% figure(2);
% bode(tfidiqref,P);
% hold on
% bode(tfidiqref_matlab,P);
% Bode_Darklines(3);
% hold off
% 
% figure(3);
% bode(tfiqidref,P);
% hold on
% bode(tfiqidref_matlab,P);
% Bode_Darklines(3);
% hold off
% 
% figure(4);
% bode(tfiqiqref,P);
% hold on
% bode(tfiqiqref_matlab,P);
% Bode_Darklines(3);
% hold off
% 
% figure(5);
% bode(tfvdidref,P);
% hold on
% bode(tfvdidref_matlab,P);
% Bode_Darklines(3);
% hold off
% 
% figure(6);
% bode(tfvdiqref,P);
% hold on
% bode(tfvdiqref_matlab,P);
% Bode_Darklines(3);
% hold off
% 
% figure(7);
% bode(tfvqidref,P);
% hold on
% bode(tfvqidref_matlab,P);
% Bode_Darklines(3);
% hold off
% 
% figure(8);
% bode(tfvqiqref,P);
% hold on
% bode(tfvqiqref_matlab,P);
% Bode_Darklines(3);
% hold off
% 
figure(9);
bode(Zoutddcl,P);
hold on
bode(Zoutddicl_matlab,P);
Bode_Darklines(3);
hold off

figure(10);
bode(Zoutdqcl,P);
hold on
bode(Zoutdqicl_matlab,P);
Bode_Darklines(3);
hold off

figure(11);
bode(Zoutqdcl,P);
hold on
bode(Zoutqdicl_matlab,P);
Bode_Darklines(3);
hold off

figure(12);
bode(Zoutqqcl,P);
hold on
bode(Zoutqqicl_matlab,P);
Bode_Darklines(3);
hold off

figure(13);
bode(Zincl,P);
hold on
bode(Zinicl_matlab,P);
Bode_Darklines(3);
hold off

figure(14);
bode(ZinwCdccl,P);
hold on
bode(ZinwCdccl_matlab,P);
Bode_Darklines(3);
hold off


% Rac=15e10;
% tfvdidref = (Kpi*Vdc*s^2*C*Rac+C*L*Rac*s^3-C*L*Rac*s*omega^2+C*Kii*Rac*Vdc*s+s^2*C*Rac*Rl+s*Vdc*Kpi+s^2*L+Kii*Vdc+s*Rac+s*Rl)*Rac*Vdc*(s*Kpi+Kii)/(C^2*Kpi^2*Rac^2*Vdc^2*s^4+C^2*Kpi^2*Rac^2*Vdc^2*s^2*omega^2+2*C^2*Kpi*L*Rac^2*Vdc*s^5+2*C^2*Kpi*L*Rac^2*Vdc*s^3*omega^2+C^2*L^2*Rac^2*s^6+2*C^2*L^2*Rac^2*s^4*omega^2+C^2*L^2*Rac^2*s^2*omega^4+2*C^2*Kii*Kpi*Rac^2*Vdc^2*s^3+2*C^2*Kii*Kpi*Rac^2*Vdc^2*s*omega^2+2*C^2*Kii*L*Rac^2*Vdc*s^4+2*C^2*Kii*L*Rac^2*Vdc*s^2*omega^2+2*C^2*Kpi*Rac^2*Rl*Vdc*s^4+2*C^2*Kpi*Rac^2*Rl*Vdc*s^2*omega^2+2*C^2*L*Rac^2*Rl*s^5+2*C^2*L*Rac^2*Rl*s^3*omega^2+C^2*Kii^2*Rac^2*Vdc^2*s^2+C^2*Kii^2*Rac^2*Vdc^2*omega^2+2*C^2*Kii*Rac^2*Rl*Vdc*s^3+2*C^2*Kii*Rac^2*Rl*Vdc*s*omega^2+C^2*Rac^2*Rl^2*s^4+C^2*Rac^2*Rl^2*s^2*omega^2+2*C*Kpi^2*Rac*Vdc^2*s^3+4*C*Kpi*L*Rac*Vdc*s^4+2*C*L^2*Rac*s^5+2*C*L^2*Rac*s^3*omega^2+4*C*Kii*Kpi*Rac*Vdc^2*s^2+4*C*Kii*L*Rac*Vdc*s^3+2*C*Kpi*Rac^2*Vdc*s^3+4*C*Kpi*Rac*Rl*Vdc*s^3+2*C*L*Rac^2*s^4-2*C*L*Rac^2*s^2*omega^2+4*C*L*Rac*Rl*s^4+2*C*Kii^2*Rac*Vdc^2*s+2*C*Kii*Rac^2*Vdc*s^2+4*C*Kii*Rac*Rl*Vdc*s^2+2*C*Rac^2*Rl*s^3+2*C*Rac*Rl^2*s^3+Kpi^2*Vdc^2*s^2+2*Kpi*L*Vdc*s^3+L^2*s^4+L^2*s^2*omega^2+2*Kii*Kpi*Vdc^2*s+2*Kii*L*Vdc*s^2+2*Kpi*Rac*Vdc*s^2+2*Kpi*Rl*Vdc*s^2+2*Rac*Rl*s^2+2*L*Rac*s^3+2*L*Rl*s^3+Kii^2*Vdc^2+2*Kii*Rac*Vdc*s+2*Kii*Rl*Vdc*s+Rl^2*s^2+Rac^2*s^2);
% G_v_ll=tfvdidref*tfdelay;		  				% 1/(sLp+RLp)
% Rac=15;
% tfvdidref = (Kpi*Vdc*s^2*C*Rac+C*L*Rac*s^3-C*L*Rac*s*omega^2+C*Kii*Rac*Vdc*s+s^2*C*Rac*Rl+s*Vdc*Kpi+s^2*L+Kii*Vdc+s*Rac+s*Rl)*Rac*Vdc*(s*Kpi+Kii)/(C^2*Kpi^2*Rac^2*Vdc^2*s^4+C^2*Kpi^2*Rac^2*Vdc^2*s^2*omega^2+2*C^2*Kpi*L*Rac^2*Vdc*s^5+2*C^2*Kpi*L*Rac^2*Vdc*s^3*omega^2+C^2*L^2*Rac^2*s^6+2*C^2*L^2*Rac^2*s^4*omega^2+C^2*L^2*Rac^2*s^2*omega^4+2*C^2*Kii*Kpi*Rac^2*Vdc^2*s^3+2*C^2*Kii*Kpi*Rac^2*Vdc^2*s*omega^2+2*C^2*Kii*L*Rac^2*Vdc*s^4+2*C^2*Kii*L*Rac^2*Vdc*s^2*omega^2+2*C^2*Kpi*Rac^2*Rl*Vdc*s^4+2*C^2*Kpi*Rac^2*Rl*Vdc*s^2*omega^2+2*C^2*L*Rac^2*Rl*s^5+2*C^2*L*Rac^2*Rl*s^3*omega^2+C^2*Kii^2*Rac^2*Vdc^2*s^2+C^2*Kii^2*Rac^2*Vdc^2*omega^2+2*C^2*Kii*Rac^2*Rl*Vdc*s^3+2*C^2*Kii*Rac^2*Rl*Vdc*s*omega^2+C^2*Rac^2*Rl^2*s^4+C^2*Rac^2*Rl^2*s^2*omega^2+2*C*Kpi^2*Rac*Vdc^2*s^3+4*C*Kpi*L*Rac*Vdc*s^4+2*C*L^2*Rac*s^5+2*C*L^2*Rac*s^3*omega^2+4*C*Kii*Kpi*Rac*Vdc^2*s^2+4*C*Kii*L*Rac*Vdc*s^3+2*C*Kpi*Rac^2*Vdc*s^3+4*C*Kpi*Rac*Rl*Vdc*s^3+2*C*L*Rac^2*s^4-2*C*L*Rac^2*s^2*omega^2+4*C*L*Rac*Rl*s^4+2*C*Kii^2*Rac*Vdc^2*s+2*C*Kii*Rac^2*Vdc*s^2+4*C*Kii*Rac*Rl*Vdc*s^2+2*C*Rac^2*Rl*s^3+2*C*Rac*Rl^2*s^3+Kpi^2*Vdc^2*s^2+2*Kpi*L*Vdc*s^3+L^2*s^4+L^2*s^2*omega^2+2*Kii*Kpi*Vdc^2*s+2*Kii*L*Vdc*s^2+2*Kpi*Rac*Vdc*s^2+2*Kpi*Rl*Vdc*s^2+2*Rac*Rl*s^2+2*L*Rac*s^3+2*L*Rl*s^3+Kii^2*Vdc^2+2*Kii*Rac*Vdc*s+2*Kii*Rl*Vdc*s+Rl^2*s^2+Rac^2*s^2);
% G_v_nl=tfvdidref*tfdelay;
% Rac=5;
% tfvdidref = (Kpi*Vdc*s^2*C*Rac+C*L*Rac*s^3-C*L*Rac*s*omega^2+C*Kii*Rac*Vdc*s+s^2*C*Rac*Rl+s*Vdc*Kpi+s^2*L+Kii*Vdc+s*Rac+s*Rl)*Rac*Vdc*(s*Kpi+Kii)/(C^2*Kpi^2*Rac^2*Vdc^2*s^4+C^2*Kpi^2*Rac^2*Vdc^2*s^2*omega^2+2*C^2*Kpi*L*Rac^2*Vdc*s^5+2*C^2*Kpi*L*Rac^2*Vdc*s^3*omega^2+C^2*L^2*Rac^2*s^6+2*C^2*L^2*Rac^2*s^4*omega^2+C^2*L^2*Rac^2*s^2*omega^4+2*C^2*Kii*Kpi*Rac^2*Vdc^2*s^3+2*C^2*Kii*Kpi*Rac^2*Vdc^2*s*omega^2+2*C^2*Kii*L*Rac^2*Vdc*s^4+2*C^2*Kii*L*Rac^2*Vdc*s^2*omega^2+2*C^2*Kpi*Rac^2*Rl*Vdc*s^4+2*C^2*Kpi*Rac^2*Rl*Vdc*s^2*omega^2+2*C^2*L*Rac^2*Rl*s^5+2*C^2*L*Rac^2*Rl*s^3*omega^2+C^2*Kii^2*Rac^2*Vdc^2*s^2+C^2*Kii^2*Rac^2*Vdc^2*omega^2+2*C^2*Kii*Rac^2*Rl*Vdc*s^3+2*C^2*Kii*Rac^2*Rl*Vdc*s*omega^2+C^2*Rac^2*Rl^2*s^4+C^2*Rac^2*Rl^2*s^2*omega^2+2*C*Kpi^2*Rac*Vdc^2*s^3+4*C*Kpi*L*Rac*Vdc*s^4+2*C*L^2*Rac*s^5+2*C*L^2*Rac*s^3*omega^2+4*C*Kii*Kpi*Rac*Vdc^2*s^2+4*C*Kii*L*Rac*Vdc*s^3+2*C*Kpi*Rac^2*Vdc*s^3+4*C*Kpi*Rac*Rl*Vdc*s^3+2*C*L*Rac^2*s^4-2*C*L*Rac^2*s^2*omega^2+4*C*L*Rac*Rl*s^4+2*C*Kii^2*Rac*Vdc^2*s+2*C*Kii*Rac^2*Vdc*s^2+4*C*Kii*Rac*Rl*Vdc*s^2+2*C*Rac^2*Rl*s^3+2*C*Rac*Rl^2*s^3+Kpi^2*Vdc^2*s^2+2*Kpi*L*Vdc*s^3+L^2*s^4+L^2*s^2*omega^2+2*Kii*Kpi*Vdc^2*s+2*Kii*L*Vdc*s^2+2*Kpi*Rac*Vdc*s^2+2*Kpi*Rl*Vdc*s^2+2*Rac*Rl*s^2+2*L*Rac*s^3+2*L*Rl*s^3+Kii^2*Vdc^2+2*Kii*Rac*Vdc*s+2*Kii*Rl*Vdc*s+Rl^2*s^2+Rac^2*s^2);
% G_v_hl=tfvdidref*tfdelay;
% figure(15)
% bode(G_v_ll,G_v_nl,G_v_hl)
% Bode_Darklines(3);
% 
% %% designed new controller from sisotool
% fprintf('New current controller design:\n')
% kpi=94.221*0.00022
% kii=94.221
% kpv=147.22*0.00016
% kiv=147.22
% R=5;
% L=1.2e-3; RL=50e-3;Rl=50e-3; 
% C=22e-6;
% Cdc=0.15e-3;
% Rcdc=1e-6;
% Ldc=0.99e-3;
% Rsdc=1e-6;
% 
% f=400; 
% w=2*pi*f; 
% fsw=20e+3; 
% Tsw=1/fsw; 
% Tdelay = 0.5*Tsw;
% Tdelay1= 0.5*Tsw;
% wsw=2*pi*fsw;
% DT=1e-6;
% 
% J=[0 -1; 1 0];
% 
% Vdc=270; Vdcref=270;
% Vse=57.5;
% Vsm=Vse*sqrt(2);
% Vsdq=[sqrt(3/2)*Vsm; 0];
% Vdref=Vsdq(1); Vqref=Vsdq(2);
% omega = 400*2*pi;
% Rac =  R;
% Vod = Vdref;
% Voq = Vqref;
% % Rl=RL;
% Dd = -(omega*L*Voq-Rac*Vod-Rl*Vod+Rl*omega*C*Voq*Rac+omega^2*L*C*Vod*Rac)/(Rac*Vdc)
% Dq = -(-Rl*Voq-Rac*Voq-omega*L*Vod+omega^2*L*C*Voq*Rac-Rl*omega*C*Vod*Rac)/(Rac*Vdc)
% Idref = -(-Vod+omega*C*Voq*Rac)/Rac
% Iqref = (Voq+omega*C*Vod*Rac)/Rac