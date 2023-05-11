clc
clear all;
% close all;

% bdclose all
% set_param(0,'CharacterEncoding','ISO-8859-1')
Tstart=0.015;
Rs=0.000001;
Vdcref=270;
alpha=1e10;         %% 0.5
Cdc=150e-6; 
R=90*alpha/(1+alpha);
RCdc=0.070;
P=Vdcref^2/R;
Vse=57.5;
Vse=Vse+P/3/Vse*Rs*1;
Vsm=Vse*sqrt(2);
L=500e-6; 
RL=100e-3;

f=400; w=2*pi*f;
fsw=20e+3; Tsw=1/fsw; wsw=2*pi*fsw;
%parameters for IGBT
Vfs= 0.5;
Vfd=0.7;
Rsd=0.05;
DT=1e-6;
DeadTime=1e-6;
PWM_Cycle=4000;
%% PLL
fpll = 1;       %% 400
damp_pll = 0.707;
a1=2*damp_pll*fpll*6.28318530717959*0.57735026918963/Vse;
a2=-(2*damp_pll*fpll*6.28318530717959-fpll*fpll*39.47841760435743*Tsw)*0.57735026918963/Vse;

DEF_pll=800;							%control loop for PLL (Hz)
DEF_pll_damp=0.707;					%damping factor for the PLL controller
DEF_Vin=57.5;						%input phase to neutral rms voltage
DEF_Tsw=0.00005;						%switching period
FWPI_a1=2*DEF_pll_damp*DEF_pll*6.28318530717959*0.57735026918963/DEF_Vin;					%digital PI parameters for frequency loop
FWPI_a2=-(2*DEF_pll_damp*DEF_pll*6.28318530717959-DEF_pll*DEF_pll*39.47841760435743*DEF_Tsw)*0.57735026918963/DEF_Vin;

tf_pll_z=tf([FWPI_a1 FWPI_a2],[1 -1],DEF_Tsw);
tf_pll=d2c(tf_pll_z)


fv=100; fi=1000;
% kpv=2*pi*fv*Cdc; 
% kiv=2*pi*fv/R;
% kpi=2*pi*fi*L/Vdcref;
% kii=2*pi*fi*RL/Vdcref;

Lboost_con=500e-6;
RLboost_con=100e-3;
Cdc_con=150e-6;
Rdc_con=90;
kpv=2*pi*fv*Cdc_con; kiv=2*pi*fv/Rdc_con;
kpi=2*pi*fi*Lboost_con/Vdcref; kii=2*pi*fi*RLboost_con/Vdcref;

Vsdq=[sqrt(3/2)*Vsm; 0];
Vsd=Vsdq(1)
Vsq=Vsdq(2)



fprintf('done!\n')


% %% model linearization:
% model = 'AFE_avg_ol_pll';
% 
% %% Linearize the model
% sys = linearize(model,0.5);
% 
% %% Plot the resulting linearization.
% % bode(sys)
% H=tf(sys);
% 
% Vavg1=[H(1,1) H(1,2);
%     H(2,1) H(2,2);];
% 
% Ilavg1=[H(3,1) H(3,2);
%     H(4,1) H(4,2);];
% 
% Isavg1=[H(5,1) H(5,2);
%     H(6,1) H(6,2);];
% Zsavg1=Vavg1/Isavg1; Zlavg1=Vavg1/Ilavg1;
% figure
% bode(Zlavg1)
% Bode_Darklines(3)
%% analitical derivation