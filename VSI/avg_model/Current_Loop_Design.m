clc
close all
clear all
%% current loop design
% system parameters
%% power stage parameters for average model
R=15;                           % load resistor
L=1.10e-3;%0.64e-3;%1.1e-3;%(0.75879e-3)/(2/3);           % output filter inductance
RL=150e-3;%150e-3;%                                       % filter inductor ESR
C=41.308e-6;%20.4e-6;%31.4e-6;%;                          % output filter capacitance
Vdcref=270;
Vdc=Vdcref;
Vse=57.5;
Vsm=Vse*sqrt(2);
Vsdq=[sqrt(3/2)*Vsm; 0];
Vdref=Vsdq(1);
% Vd reference step factor
Vdstep=1.2;
Vqref=Vsdq(2);
fline=400; 
w=2*pi*fline;
omega=w;
fsw=20e+3;
Tsw=1/fsw;
wsw=2*pi*fsw;
Tdelay=1.5*Tsw;
s=tf([1,0],[0,1]);
% omega=400*2*pi;
tfdelay = (1-Tdelay*s)/(1+Tdelay*s);
Rl=RL;
% id to dd transfer function:
Rac=1.5e3;
tfidd_ln = (L*C^2*Rac^2*s^3+2*s^2*C*Rac*L+s*L+L*s*C^2*Rac^2*omega^2+C^2*Rac^2*Rl*s^2+2*s*C*Rac*Rl+s*C*Rac^2+C^2*Rac^2*omega^2*Rl+Rac+Rl)*Vdc/(L^2*C^2*Rac^2*s^4+2*L^2*C*Rac*s^3+s^2*L^2+2*L^2*s^2*C^2*Rac^2*omega^2+2*L^2*C*Rac*omega^2*s+omega^2*L^2+L^2*C^2*Rac^2*omega^4+2*L*C^2*Rac^2*Rl*s^3+2*L*s^2*C*Rac^2+4*L*C*Rac*Rl*s^2+2*L*s*C^2*Rac^2*omega^2*Rl+2*L*s*Rac+2*s*L*Rl-2*L*omega^2*C*Rac^2+C^2*Rac^2*Rl^2*s^2+2*s*C*Rac^2*Rl+2*C*Rac*Rl^2*s+Rl^2+Rac^2+2*Rac*Rl+C^2*Rac^2*omega^2*Rl^2);
G_i_ll=tfidd_ln*tfdelay;		  				% 1/(sLp+RLp)
Rac=15;
tfidd_ln = (L*C^2*Rac^2*s^3+2*s^2*C*Rac*L+s*L+L*s*C^2*Rac^2*omega^2+C^2*Rac^2*Rl*s^2+2*s*C*Rac*Rl+s*C*Rac^2+C^2*Rac^2*omega^2*Rl+Rac+Rl)*Vdc/(L^2*C^2*Rac^2*s^4+2*L^2*C*Rac*s^3+s^2*L^2+2*L^2*s^2*C^2*Rac^2*omega^2+2*L^2*C*Rac*omega^2*s+omega^2*L^2+L^2*C^2*Rac^2*omega^4+2*L*C^2*Rac^2*Rl*s^3+2*L*s^2*C*Rac^2+4*L*C*Rac*Rl*s^2+2*L*s*C^2*Rac^2*omega^2*Rl+2*L*s*Rac+2*s*L*Rl-2*L*omega^2*C*Rac^2+C^2*Rac^2*Rl^2*s^2+2*s*C*Rac^2*Rl+2*C*Rac*Rl^2*s+Rl^2+Rac^2+2*Rac*Rl+C^2*Rac^2*omega^2*Rl^2);
G_i_nl=tfidd_ln*tfdelay;
Rac=5;
tfidd_ln = (L*C^2*Rac^2*s^3+2*s^2*C*Rac*L+s*L+L*s*C^2*Rac^2*omega^2+C^2*Rac^2*Rl*s^2+2*s*C*Rac*Rl+s*C*Rac^2+C^2*Rac^2*omega^2*Rl+Rac+Rl)*Vdc/(L^2*C^2*Rac^2*s^4+2*L^2*C*Rac*s^3+s^2*L^2+2*L^2*s^2*C^2*Rac^2*omega^2+2*L^2*C*Rac*omega^2*s+omega^2*L^2+L^2*C^2*Rac^2*omega^4+2*L*C^2*Rac^2*Rl*s^3+2*L*s^2*C*Rac^2+4*L*C*Rac*Rl*s^2+2*L*s*C^2*Rac^2*omega^2*Rl+2*L*s*Rac+2*s*L*Rl-2*L*omega^2*C*Rac^2+C^2*Rac^2*Rl^2*s^2+2*s*C*Rac^2*Rl+2*C*Rac*Rl^2*s+Rl^2+Rac^2+2*Rac*Rl+C^2*Rac^2*omega^2*Rl^2);
G_i_hl=tfidd_ln*tfdelay;
% controller design
figure(15)
bode(G_i_ll,G_i_nl,G_i_hl)
Bode_Darklines(3);
grid on

wc=2*pi*fsw/35;                                         % desired bandwidth
mfi=30/180*pi;                                          % desired phase margin
% derivation of PI parameters
[mag,phase]=bode(G_i_ll,wc);
Tz=1/wc*tan(mfi-pi/2-phase/180*pi);                     % zero Tz of PI control
Ti=abs(mag*(1+1i*wc*Tz)/wc);                            % time constant Ti of PI control
% check of the design
GH_ll=G_i_ll*tf([Tz  1],[Ti  0]);					% open loop gain
WW_ll=feedback(GH_ll,1,-1);						% transfer function bw i_ref and i
% figure(1)
% step(WW_ll)										% step response
% grid on
% figure(2)
% margin(GH_ll)
% grid on

kpin=Tz/Ti;
kiin=1/Ti;
kpi=2*pi*800*L/Vdc*0.5; 
kii=2*pi*800*(RL)/Vdc*10;
% kpi=kpin; 
% kii=kiin;



tf_comp=tf([1 2*pi*2e3],[1 2*pi*1]);

GH_ll=G_i_ll*tf([kpi  kii],[1  0]);					% open loop gain
GH_ll_comp=GH_ll*tf_comp;
WW_ll=feedback(GH_ll,1,-1);						% transfer function bw i_ref and i
WW_ll_comp=feedback(GH_ll_comp,1,-1)
figure(3)
step(WW_ll)										% step response
hold on
% step(WW_ll_comp)
grid on
figure(4)
margin(GH_ll)
hold on
% margin(GH_ll_comp)
grid on
%% normal load loop gain
GH_nl=G_i_nl*tf([Tz  1],[Ti  0]);

WW_nl=feedback(GH_nl,1,-1);						% transfer function bw i_ref and i

% figure(5)
% step(WW_nl)										% step response
% grid on
% figure(6)
% margin(GH_nl)
% grid on

GH_nl=G_i_nl*tf([kpi  kii],[1  0]);
GH_nl_comp=GH_nl*tf_comp;
WW_nl=feedback(GH_nl,1,-1);						% transfer function bw i_ref and i
WW_nl_comp=feedback(GH_nl_comp,1,-1)
figure(7)
step(WW_nl)										% step response
hold on
% step(WW_nl_comp)
grid on
figure(8)
margin(GH_nl)
hold on
% margin(GH_nl_comp)
grid on

%% high load loop gain
GH_hl=G_i_hl*tf([Tz  1],[Ti  0]);
WW_hl=feedback(GH_hl,1,-1);						% transfer function bw i_ref and i
% figure(9)
% step(WW_hl)	
% grid on% step response
% figure(10)
% margin(GH_hl)
% grid on

GH_hl=G_i_hl*tf([kpi  kii],[1  0]);
WW_hl=feedback(GH_hl,1,-1);						% transfer function bw i_ref and i
% figure(11)
% step(WW_hl)	
% grid on% step response
% figure(12)
% margin(GH_hl)
% grid on

% kpi=kpin;
% kii=kiin;



Ra=16.08;
Rb=16.3;
Rc=15.46;
L=(0.75879e-3)/(2/3);%1.2e-3;
RL=80e-3; 
C=21.308e-6;
RC=0.054;

Vdc=270;
Vse=57.5;
Vsm=Vse*sqrt(2);
Vsdq=[sqrt(3/2)*Vsm; 0];
Vdref=Vsdq(1);
% Vd reference step factor
Vdstep=1.2;
Vqref=Vsdq(2);
Iin_limit = 15*sqrt(2);
fline=400; 
w=2*pi*fline;
fsw=20e+3;
Tsw=1/fsw;
wsw=2*pi*fsw;
%parameters for IGBT 6MBP30RH060-50
Vfs= 0.5;
Vfd= 0.7;
% suppose R IGBT and R Diode are the same
Rsd= 0.05;
Rss= Rsd;
DT = 80;
DeadTime = 1e-6;
DT_en = 0;
omega = 400*2*pi;
Rl = RL;
Rac =  R;
Vod = Vdref;
Voq = Vqref;
% Rl = RL;
Dd = -(omega*L*Voq-Rac*Vod-Rl*Vod+Rl*omega*C*Voq*Rac+omega^2*L*C*Vod*Rac)/(Rac*Vdc);
Dq = -(-Rl*Voq-Rac*Voq-omega*L*Vod+omega^2*L*C*Voq*Rac-Rl*omega*C*Vod*Rac)/(Rac*Vdc);
Idref = -(-Vod+omega*C*Voq*Rac)/Rac;
Iqref = (Voq+omega*C*Vod*Rac)/Rac;
% initial for PWM generation
PWM_Cycle = 4000;  %time 12.5ns clock
PWM_mode = 0;
% fv=100; fi=1000;
% 
% % controller 1
% kpi=2*pi*fi*L;%/Vdc; 
% kii=2*pi*fi*(RL)%;/Vdc;
% kpv=2*pi*fv*C;
% kiv=2*pi*fv/R;

% controller 2
% kpi=0.0279;%2*pi*fi*L/Vdc; 
% kii=1.1636;%2*pi*fi*(RL)/Vdc;
% kpv=0.0138*3;%2*pi*fv*C;
% kiv=41.8879*3;%2*pi*fv/R;


%----parameter for DQ component filtering
% kpi=94.221*0.00022*270
% kii=94.221*270
% kpv=147.22*0.00016
% kiv=147.22

%% initial for controller parameters
% Time_stamp = 0;
% control_para=[Time_stamp R C L RL Vfs Rss Vfd Rsd fline Tsw DeadTime DT_en Iin_limit ...
%     Vdref Vdstep Vqref Vdc kpv kiv kpi kii PWM_mode];

fprintf('\ndone\n');



