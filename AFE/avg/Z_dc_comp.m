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
DEF_pll=300;							%control loop for PLL (Hz)
DEF_pll_damp=0.707;                     %damping factor for the PLL controller
DEF_Vin=57.5;                           %input phase to neutral rms voltage
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

sys = linearize('AFE_avg_vl_pll_Zdc',0.5);

H=ss(sys);
Zo_dc_avg_sim_100=H(1)/(3);



fv=50; fi=1000;
Lboost_con=500e-6;
RLboost_con=100e-3;
Cdc_con=150e-6;
Rdc_con=90;
kpv=2*pi*fv*Cdc_con; kiv=2*pi*fv/Rdc_con;
kpi=2*pi*fi*Lboost_con/Vdcref; kii=2*pi*fi*RLboost_con/Vdcref;

sys = linearize('AFE_avg_vl_pll_Zdc',0.5);

H=ss(sys);
Zo_dc_avg_sim_50=H(1)/(3);


fv=50; fi=500;
Lboost_con=500e-6;
RLboost_con=100e-3;
Cdc_con=150e-6;
Rdc_con=90;
kpv=2*pi*fv*Cdc_con; kiv=2*pi*fv/Rdc_con;
kpi=2*pi*fi*Lboost_con/Vdcref; kii=2*pi*fi*RLboost_con/Vdcref;

sys = linearize('AFE_avg_vl_pll_Zdc',0.5);

H=ss(sys);
Zo_dc_avg_sim_50_500=H(1)/(3);


Cdc=35e-6; 
fv=50; fi=1000;
Lboost_con=500e-6;
RLboost_con=100e-3;
Cdc_con=35e-6;
Rdc_con=90;
kpv=2*pi*fv*Cdc_con; kiv=2*pi*fv/Rdc_con;
kpi=2*pi*fi*Lboost_con/Vdcref; kii=2*pi*fi*RLboost_con/Vdcref;

sys = linearize('AFE_avg_vl_pll_Zdc',0.5);

H=ss(sys);
Zo_dc_avg_sim_50_500_50=H(1)/(3);

%% calculation and simulation comparison
figure(1)
bode(Zo_dc_avg_sim_100,Zo_dc_avg_sim_50,Zo_dc_avg_sim_50_500,...
    Zo_dc_avg_sim_50_500_50,Bode_O)
legend('Zo_dc_avg_sim_100','Zo_dc_avg_sim_50','Zo_dc_avg_sim_50_500',...
    'Zo_dc_avg_sim_50_500_50')
Bode_Darklines(3)
