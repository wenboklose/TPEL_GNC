clc
clear all;
close all

Bode_O=bodeoptions;
Bode_O.XLabel.FontSize=14;
Bode_O.YLabel.FontSize=14;
Bode_O.TickLabel.FontSize=14;
Bode_O.Title.FontSize=14;
Bode_O.title.String=' ';
Bode_O.Grid='on';
Bode_O.XLim={[1 1e4]};
Bode_O.XLimMode={'manual'};
Bode_O.FreqUnits='Hz';
Bode_O.PhaseWrapping='off';

% bdclose all
% set_param(0,'CharacterEncoding','ISO-8859-1')
Tstart=0.015;
Rs=1e-6;
Vdcref=270;
alpha=1e10;         %% 0.5
Cdc=150e-6; 
R=90;
RCdc=0.070;
P=Vdcref^2/R;
Vse=57.5;
Vse=Vse+P/3/Vse*Rs*1;
Vsm=Vse*sqrt(2);
Vsm_o=Vse*sqrt(2);
L=510e-6; 
RL=150e-3;

f=400; w=2*pi*f;
fsw=20e+3; Tsw=1/fsw; wsw=2*pi*fsw;
% %parameters for IGBT
% Vfs= 0.5;
% Vfd=0.7;
% Rsd=0.05;
% DT=1e-6;
% DeadTime=1e-6;
% PWM_Cycle=4000;
%% PLL
DEF_pll=10;							%control loop for PLL (Hz)
DEF_pll_damp=0.707;					%damping factor for the PLL controller
DEF_Vin=57.5;						%input phase to neutral rms voltage
DEF_Tsw=0.00005;						%switching period
FWPI_a1=2*DEF_pll_damp*DEF_pll*6.28318530717959*0.57735026918963/DEF_Vin;					%digital PI parameters for frequency loop
FWPI_a2=-(2*DEF_pll_damp*DEF_pll*6.28318530717959-DEF_pll*DEF_pll*39.47841760435743*DEF_Tsw)*0.57735026918963/DEF_Vin;
tf_pll_z=tf([FWPI_a1 FWPI_a2],[1 -1],DEF_Tsw);
tf_pll=d2c(tf_pll_z);

fv=100; fi=1000;
Lboost_con=500e-6;
RLboost_con=100e-3;
Cdc_con=150e-6;
Rdc_con=90;
kpv=2*pi*fv*Cdc_con; kiv=2*pi*fv/Rdc_con;
kpi=2*pi*fi*Lboost_con/Vdcref; kii=2*pi*fi*RLboost_con/Vdcref;

Vsdq=[sqrt(3/2)*Vsm; 0];
Vsd=Vsdq(1);
Vsq=Vsdq(2);
Dd=0.3657;
Dq=-0.0382;
Id=8.2109;
Iq=0.0227;
fprintf('initialization is done!\n')

%% Zin open loop without pll model linearization:

sys = linearize('AFE_avg_Zol',0.5);

H=ss(sys);
Vavg1=[H(1,1) H(1,2);
    H(2,1) H(2,2);];
Ilavg1=[H(3,1) H(3,2);
    H(4,1) H(4,2);];
Zin_ol_avg_sim=Vavg1/Ilavg1;

%% Zin open loop with pll model linearization:

sys = linearize('AFE_avg_Zol_pll',0.5);

H=ss(sys);
Vavg1=[H(1,1) H(1,2);
    H(2,1) H(2,2);];
Ilavg1=[H(3,1) H(3,2);
    H(4,1) H(4,2);];
Zin_ol_pll_avg_sim=Vavg1/Ilavg1;

%% Zin closed loop with pll model linearization:

[sys op] = linearize('AFE_avg_vl_pll',0.5);

H=ss(sys);
Vavg1=[H(1,1) H(1,2);
    H(2,1) H(2,2);];
Ilavg1=[H(3,1) H(3,2);
    H(4,1) H(4,2);];
Zin_vl_pll_avg_sim=Vavg1/Ilavg1;

% %% Zin closed loop without pll model linearization:
% 
% sys = linearize('AFE_avg_vl',0.5);
% 
% H=ss(sys);
% Vavg1=[H(1,1) H(1,2);
%     H(2,1) H(2,2);];
% Ilavg1=[H(3,1) H(3,2);
%     H(4,1) H(4,2);];
% Zin_vl_avg_sim=Vavg1/Ilavg1;

%% Calculation based on model
s=tf([1 0],[0 1]);
Zdc=R/(R*Cdc*s+1);
GL=1*L*s+RL;
omega=w;
Vdc=270.1678;
Zin_cal=[Zdc*Dd^2+GL Zdc*Dd*Dq-1*omega*L; Zdc*Dd*Dq+1*omega*L Zdc*Dq^2+GL];
Yin_cal=1/Zin_cal;
Dt=[Dd Dq;0 0];
It=[Id Iq;0 0];
Gve_cal=Zdc*Dt*Yin_cal;
den1=(-Zdc*Dd^2-GL)*(-Zdc*Dq^2-GL)-(-1*omega*L-Zdc*Dd*Dq)*(1*omega*L-Zdc*Dd*Dq);
num11=(-Zdc*Id*Dd-Vdc)*(-Zdc*Dq^2-GL)-(-Zdc*Id*Dq)*(1*omega*L-Zdc*Dd*Dq);
num12=(-Zdc*Dd*Iq)*(-Zdc*Dq^2-GL)-(-Zdc*Dq*Iq-Vdc)*(1*omega*L-Zdc*Dd*Dq);
den2=-den1;
num21=(-Zdc*Id*Dd-Vdc)*(-1*omega*L-Zdc*Dd*Dq)-(-Zdc*Id*Dq)*(-Zdc*Dd^2-GL);
num22=(-Zdc*Dd*Iq)*(-1*omega*L-Zdc*Dd*Dq)-(-Zdc*Dq*Iq-Vdc)*(-Zdc*Dd^2-GL);
Gid_cal=[-num11/den1 -num12/den1;-num21/den2 -num22/den2];
Gvd_cal=Zdc*(Dt*Gid_cal+1*It);

% Vse=57.5;
% Vsm=Vse*sqrt(2);
% Vsdq=[sqrt(3/2)*Vsm; 0];
% Vsd=Vsdq(1)
% Vsq=Vsdq(2)
E0=Vsd;
Gpll=tf_pll/(s+E0*tf_pll);
Gipll=[0 Iq*Gpll;0 -Id*Gpll];
Gdpll=[0 -Dq*Gpll;0 +Dd*Gpll];

% Zin_ol_pll=1/(Gppll+Yin_ol*Gspll);
Zin_pll_cal=1/(1/Zin_cal+Gid_cal*Gdpll);

%% Calculatoin with current and voltage control loop
Gdei = [0 1/Vdcref*Lboost_con*w; -1/Vdcref*Lboost_con*w 0];
Gci = -[kpi+kii/s 0; 0 kpi+kii/s];
Ki = [1 0; 0 1];
Gcv = [kpv+kiv/s 0; 0 0];
Kv = [1 0; 0 0];
Gdel = [1 0; 0 1];
I = [1 0; 0 1];
A = Gid_cal/(I+Gdel*(-Gdei+Gci)*Gcv*Kv*Gvd_cal);
B = Gdel*(-Gdei+Gci)*Ki;
C = Gdpll-Gdel*(-Gdei+Gci)*(Gipll+Gcv*Kv*Gve_cal);
Yin_vil_pll = (I-A*B)\(A*C+Yin_cal);

%% calculation and simulation comparison
figure(1)
bode(Zin_ol_avg_sim,Zin_cal,Bode_O)
legend('Zin\_avg\_sim','Zin\_cal')
Bode_Darklines(3)

figure(2)
bode(Zin_ol_pll_avg_sim,Zin_pll_cal,Bode_O)
legend('Zin\_pll\_avg\_sim','Zin\_pll\_cal')
Bode_Darklines(3)

figure(3)
bode(Zin_vl_pll_avg_sim,Bode_O)
legend('Zin\_pll\_avg\_sim','Zin\_pll\_cal')
Bode_Darklines(3)

figure(4)
bode(1/Zin_vl_pll_avg_sim,Bode_O)
Bode_Darklines(3)

figure(5)
nyquist(1/Zin_vl_pll_avg_sim)
Bode_Darklines(3)

% Yin_vl_avg_sim=1/Zin_vl_avg_sim;
% 
% save('Yin_vl','Yin_vl_avg_sim');