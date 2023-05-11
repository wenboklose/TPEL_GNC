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
L=500e-6; 
RL=100e-3;
I = [1 0; 0 1];
f=400; w=2*pi*f;
omega=w;
fsw=20e+3; Tsw=1/fsw; wsw=2*pi*fsw;
s=tf([1 0],[0 1]);
%% PLL
DEF_pll=100;							%control loop for PLL (Hz)
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

Id=8.1579;
Iq=0.0227;

Hidd = (-.3353*s^4-0.356e5*s^3+0.243e10*s^2+0.3791e15*s+0.1094e20)/(s^4+0.2257e6*s^3+0.1939e11*s^2+0.7515e15*s+0.1108e20);
Hidq = (0.1007e9*s^2+0.2217e14*s+0.9166e18)/(s^4+0.2257e6*s^3+0.1939e11*s^2+0.7515e15*s+0.1108e20);
Hi = [Hidd 0*Hidq;-0*Hidq Hidd];
Hvdd = (0.2769e-1*s^4-5135*s^3+0.3403e9*s^2+0.1703e14*s+0.3006e18)/(s^4+0.7407e5*s^3+0.249e10*s^2+0.4142e14*s+0.3073e18);
Hvdq = (-0.1806e8*s^2+0.2702e13*s+0.599e17)/(s^4+0.7407e5*s^3+0.249e10*s^2+0.4142e14*s+0.3073e18);
Hv = [Hvdd Hvdq; -Hvdq Hvdd];
Gdeldd = exp((-1/fsw)*s)*cos(omega*(1/fsw));
Gdeldq = exp((-1/fsw)*s)*sin(omega*(1/fsw));
Gsvmdd = (1-.25*s/fsw)/(1+.25*s/fsw);%exp(-0.5/fsw*s);
Gsvmdq = -(1-.25*s/fsw)*0/(1+.25*s/fsw);
Gsvm = [Gsvmdd 0*Gsvmdq;-0*Gsvmdq Gsvmdd];
Gdel = Gsvm;%[Gdeldd 1*Gdeldq;-1*Gdeldq Gdeldd];
Ki = I;%Hi;%[1 0; 0 1];

fprintf('initialization is done!\n')

%% Zin model linearization:

% sys = linearize('AFE_avg_il',0.5);
sys = linearize('AFE_avg_il_Ki_SVM',0.5);
%% Zin model linearization
H=ss(sys);
Vavg1=[H(1,1) H(1,2);
    H(2,1) H(2,2);];
Ilavg1=[H(3,1) H(3,2);
    H(4,1) H(4,2);];
Zin_il_avg_sim=Vavg1/Ilavg1;
Dd = Dd(length(Dd)-10);
Dq = Dq(length(Dq)-10);
Vdc = Vdc(length(Vdc)-10);

%% Calculation based on model
E0=Vsd;
Gpll=tf_pll/(s+E0*tf_pll);
Gipll=[0 Iq*Gpll;0 -Id*Gpll];
Gdpll=[0 -Dq*Gpll;0 +Dd*Gpll];
Zdc=R/(R*Cdc*s+1);
GL=1*L*s+RL;
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

% Zin_ol_pll=1/(Gppll+Yin_ol*Gspll);
Zin_pll_cal=1/(1/Zin_cal+Gid_cal*Gdpll);
Ypll=Gid_cal*Gdpll;
Yin_pll_cal=Yin_cal+Ypll;

%% Calculatoin with current control loop
Gdei = [0 1/Vdcref*Lboost_con*w; -1/Vdcref*Lboost_con*w 0];
Gci = -[kpi+kii/s 0; 0 kpi+kii/s];

% Gdel = [1 0; 0 1];

Ypll_il = (-Gid_cal*(Gdel*(-Gdei + Gci)*Gipll-Gdpll));
T = (I+Gid_cal*Gdel*(-Gdei + Gci)*Ki);
Yin_il_avg_cal = T\Yin_cal;
Yin_il_pll_avg_cal = T\(Ypll_il+Yin_cal);
Zin_il_pll_avg_cal = (Ypll_il+Yin_cal)\T;
Ypll_il_T = T\Ypll_il;
Yin_cal_T = T\Yin_cal;

figure(2)
bode(1/Zin_il_avg_sim,Yin_il_avg_cal,Bode_O)
legend('Yin\_il\_avg\_sim','Yin\_il\_avg\_cal')
Bode_Darklines(3)
figure(21)
bode(Zin_il_avg_sim,1/Yin_il_avg_cal,Bode_O)
legend('Zin\_il\_avg\_sim','Zin\_il\_avg\_cal')
Bode_Darklines(3)
%% Zin_pll model linearization:

% sys = linearize('AFE_avg_il_pll',0.5);
sys = linearize('AFE_avg_il_pll_Ki_SVM',0.5);
%% Zin model linearization
H=ss(sys);
Vavg1=[H(1,1) H(1,2);
    H(2,1) H(2,2);];
Ilavg1=[H(3,1) H(3,2);
    H(4,1) H(4,2);];
Zin_il_pll_avg_sim=Vavg1/Ilavg1;

Dd = Dd(length(Dd)-10);
Dq = Dq(length(Dq)-10);
Vdc = Vdc(length(Vdc)-10);
E0=Vsd;
Gpll=tf_pll/(s+E0*tf_pll);
Gipll=[0 Iq*Gpll;0 -Id*Gpll];
Gdpll=[0 -Dq*Gpll;0 +Dd*Gpll];
%% Calculation based on model
Zdc=R/(R*Cdc*s+1);
GL=1*L*s+RL;
omega=w;
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

% Zin_ol_pll=1/(Gppll+Yin_ol*Gspll);
Zin_pll_cal=1/(1/Zin_cal+Gid_cal*Gdpll);
Ypll=Gid_cal*Gdpll;
Yin_pll_cal=Yin_cal+Ypll;

%% Calculatoin with current control loop
I = [1 0; 0 1];
Gdei = 0*[0 1/Vdcref*Lboost_con*w; -1/Vdcref*Lboost_con*w 0];
Gci = -[kpi+kii/s 0; 0 kpi+kii/s];

% Gdel = [1 0; 0 1];

Ypll_il = (-Gid_cal*(Gdel*(-Gdei + Gci)*Gipll-Gdpll));
T = (I+Gid_cal*Gdel*(-Gdei + Gci)*Ki);
Yin_il_avg_cal = T\Yin_cal;
Yin_il_pll_avg_cal = T\(Ypll_il+Yin_cal);
Zin_il_pll_avg_cal = T/(Ypll_il+Yin_cal);
Zin_il_pll_avg_cal = (Ypll_il+Yin_cal)\T;
Ypll_il_T = T\Ypll_il;
Yin_cal_T = T\Yin_cal;
%% calculation and simulation comparison



figure(3)
bode(1/Zin_il_pll_avg_sim,Yin_il_pll_avg_cal,Bode_O)
legend('Yin\_il\_avg\_pll\_sim','Yin\_il\_avg\_pll\_cal')
Bode_Darklines(3)

figure(31)
bode(Zin_il_pll_avg_sim,Zin_il_pll_avg_cal,Bode_O)
legend('Zin\_il\_avg\_pll\_sim','Zin\_il\_avg\_pll\_cal')
Bode_Darklines(3)