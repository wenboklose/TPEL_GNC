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
Rs=1e-1;
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

Id=8.1579;
Iq=0.0257;

%% signal conditioning filter
C1=75e-12;C2=39e-12;R1=15e3;R2=15e3;
tf_filter = tf([0 1],[C1*C2*R1*R2 C2*(R1+R2) 1]);
tf_filter_dq = tf(JF_DQFromABC(tf_filter,omega));
tf_filter_dq = [tf_filter_dq(1,1) 1*tf_filter_dq(1,2); -1*tf_filter_dq(1,2) tf_filter_dq(2,2)];

Gdeldd = exp((-1/fsw)*s)*cos(omega*(1/fsw));
Gdeldq = exp((-1/fsw)*s)*sin(omega*(1/fsw));
Tdelay = 1.5/fsw;
Gsvmdd = (1-.5*Tdelay*s)/(1+.5*Tdelay*s);%exp(-0.5/fsw*s);
Gsvmdq = -(1-.5*Tdelay*s)*0/(1+.5*Tdelay*s);
Gsvm = [Gsvmdd 0*Gsvmdq;-0*Gsvmdq Gsvmdd];

fprintf('initialization is done!\n')

% %% Zin model linearization:
% 
% % sys = linearize('AFE_avg_il',0.5);
% sys = linearize('AFE_avg_il_Ki_SVM',0.5);
% %% Zin model linearization
% H=ss(sys);
% Vavg1=[H(1,1) H(1,2);
%     H(2,1) H(2,2);];
% Ilavg1=[H(3,1) H(3,2);
%     H(4,1) H(4,2);];
% Zin_il_avg_sim=Vavg1/Ilavg1;
% Dd = Dd(length(Dd)-10);
% Dq = Dq(length(Dq)-10);
% Vdc = Vdc(length(Vdc)-10);
% 
% %% Calculation based on model
% E0=Vsd;
% Gpll=tf_pll/(s+E0*tf_pll);
% Gipll=[0 Iq*Gpll;0 -Id*Gpll];
% Gdpll=[0 -Dq*Gpll;0 +Dd*Gpll];
% Zdc=R/(R*Cdc*s+1);
% GL=1*L*s+RL;
% Zin_cal=[Zdc*Dd^2+GL Zdc*Dd*Dq-1*omega*L; Zdc*Dd*Dq+1*omega*L Zdc*Dq^2+GL];
% Yin_cal=1/Zin_cal;
% Dt=[Dd Dq;0 0];
% It=[Id Iq;0 0];
% Gve_cal=Zdc*Dt*Yin_cal;
% den1=(-Zdc*Dd^2-GL)*(-Zdc*Dq^2-GL)-(-1*omega*L-Zdc*Dd*Dq)*(1*omega*L-Zdc*Dd*Dq);
% num11=(-Zdc*Id*Dd-Vdc)*(-Zdc*Dq^2-GL)-(-Zdc*Id*Dq)*(1*omega*L-Zdc*Dd*Dq);
% num12=(-Zdc*Dd*Iq)*(-Zdc*Dq^2-GL)-(-Zdc*Dq*Iq-Vdc)*(1*omega*L-Zdc*Dd*Dq);
% den2=-den1;
% num21=(-Zdc*Id*Dd-Vdc)*(-1*omega*L-Zdc*Dd*Dq)-(-Zdc*Id*Dq)*(-Zdc*Dd^2-GL);
% num22=(-Zdc*Dd*Iq)*(-1*omega*L-Zdc*Dd*Dq)-(-Zdc*Dq*Iq-Vdc)*(-Zdc*Dd^2-GL);
% Gid_cal=[-num11/den1 -num12/den1;-num21/den2 -num22/den2];
% Gvd_cal=Zdc*(Dt*Gid_cal+1*It);
% 
% % Zin_ol_pll=1/(Gppll+Yin_ol*Gspll);
% Zin_pll_cal=1/(1/Zin_cal+Gid_cal*Gdpll);
% Ypll=Gid_cal*Gdpll;
% Yin_pll_cal=Yin_cal+Ypll;
% 
% %% Calculatoin with current control loop
% Gdei = [0 1/Vdcref*Lboost_con*w; -1/Vdcref*Lboost_con*w 0];
% Gci = -[kpi+kii/s 0; 0 kpi+kii/s];
% 
% % Gdel = [1 0; 0 1];
% 
% Ypll_il = (-Gid_cal*(Gdel*(-Gdei + Gci)*Gipll-Gdpll));
% T = (I+Gid_cal*Gdel*(-Gdei + Gci)*Ki);
% Yin_il_avg_cal = T\Yin_cal;
% Zin_il_avg_cal = Yin_cal\T;
% Yin_cal_T = T\Yin_cal;
% 
% figure(2)
% bode(1/Zin_il_avg_sim,Yin_il_avg_cal,Bode_O)
% legend('Yin\_il\_avg\_sim','Yin\_il\_avg\_cal')
% Bode_Darklines(3)
% figure(21)
% bode(Zin_il_avg_sim,Zin_il_avg_cal,Bode_O)
% legend('Zin\_il\_avg\_sim','Zin\_il\_avg\_cal')
% Bode_Darklines(3)
%% Zin_pll model linearization:

% sys = linearize('AFE_avg_il_pll',0.5);
I = [tf([0 1],[0 1]) tf([0 0],[0 1]); tf([0 0],[0 1]) tf([0 1],[0 1])];
Gdel = Gsvm;
Ki = tf_filter_dq;
Kv = tf_filter_dq;%Hv;
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
Id = Id(length(Id)-10);
Iq = 1*Iq(length(Iq)-10);
Vsd = Vsd(length(Vsd)-10);
Vsq = Vsq(length(Vsq)-10);
Vdc = Vdc(length(Vdc)-10);


%% Calculation based on model of frd
n= 1e3;
w = logspace(-2,8,n);
w = 2*pi*w;

E0=Vsd;
Gpll=tf_pll/(s+E0*tf_pll);
Gipll=[tf([0 0],[0 1]) Iq*Gpll;tf([0 0],[0 1]) -Id*Gpll];

Gdel_frd = [frd(freqresp(Gdel(1,1),w),w) frd(freqresp(Gdel(1,2),w),w);...
    frd(freqresp(Gdel(2,1),w),w) frd(freqresp(Gdel(2,2),w),w)];

Ki_frd = [frd(freqresp(Ki(1,1),w),w) frd(freqresp(Ki(1,2),w),w);...
    frd(freqresp(Ki(2,1),w),w) frd(freqresp(Ki(2,2),w),w)];

Kv_frd = [frd(freqresp(Kv(1,1),w),w) frd(freqresp(Kv(1,2),w),w);...
    frd(freqresp(Kv(2,1),w),w) frd(freqresp(Kv(2,2),w),w)];

Gipll_frd = [frd(freqresp(Gipll(1,1),w),w) frd(freqresp(Gipll(1,2),w),w);...
    frd(freqresp(Gipll(2,1),w),w) frd(freqresp(Gipll(2,2),w),w)];

Gdpll=[tf([0 0],[0 1]) -Dq*Gpll;tf([0 0],[0 1]) +Dd*Gpll];
Gdpll_frd = [frd(freqresp(Gdpll(1,1),w),w) frd(freqresp(Gdpll(1,2),w),w);...
    frd(freqresp(Gdpll(2,1),w),w) frd(freqresp(Gdpll(2,2),w),w)];
Zdc=R/(R*Cdc*s+1);
GL=1*L*s+RL;
Zin_cal=[Zdc*Dd^2+GL Zdc*Dd*Dq-1*omega*L; Zdc*Dd*Dq+1*omega*L Zdc*Dq^2+GL];
Yin_cal=1/Zin_cal;
Yin_cal_frd = [frd(freqresp(Yin_cal(1,1),w),w) frd(freqresp(Yin_cal(1,2),w),w);...
    frd(freqresp(Yin_cal(2,1),w),w) frd(freqresp(Yin_cal(2,2),w),w)];

Dt=[Dd Dq;0 0];
It=[Id Iq;0 0];
Gve_cal=Zdc*Dt*Yin_cal;
Gve_cal_frd = [frd(freqresp(Gve_cal(1,1),w),w) frd(freqresp(Gve_cal(1,2),w),w);...
    frd(freqresp(Gve_cal(2,1),w),w) frd(freqresp(Gve_cal(2,2),w),w)];
den1=(-Zdc*Dd^2-GL)*(-Zdc*Dq^2-GL)-(-1*omega*L-Zdc*Dd*Dq)*(1*omega*L-Zdc*Dd*Dq);
num11=(-Zdc*Id*Dd-Vdc)*(-Zdc*Dq^2-GL)-(-Zdc*Id*Dq)*(1*omega*L-Zdc*Dd*Dq);
num12=(-Zdc*Dd*Iq)*(-Zdc*Dq^2-GL)-(-Zdc*Dq*Iq-Vdc)*(1*omega*L-Zdc*Dd*Dq);
den2=-den1;
num21=(-Zdc*Id*Dd-Vdc)*(-1*omega*L-Zdc*Dd*Dq)-(-Zdc*Id*Dq)*(-Zdc*Dd^2-GL);
num22=(-Zdc*Dd*Iq)*(-1*omega*L-Zdc*Dd*Dq)-(-Zdc*Dq*Iq-Vdc)*(-Zdc*Dd^2-GL);
Gid_cal=[-num11/den1 -num12/den1;-num21/den2 -num22/den2];
Gid_cal_frd = [frd(freqresp(Gid_cal(1,1),w),w) frd(freqresp(Gid_cal(1,2),w),w);...
    frd(freqresp(Gid_cal(2,1),w),w) frd(freqresp(Gid_cal(2,2),w),w)];
Gvd_cal=Zdc*(Dt*Gid_cal+1*It);
Gvd_cal_frd = [frd(freqresp(Gvd_cal(1,1),w),w) frd(freqresp(Gvd_cal(1,2),w),w);...
    frd(freqresp(Gvd_cal(2,1),w),w) frd(freqresp(Gvd_cal(2,2),w),w)];

% Ypll_frd=Gid_cal_frd*Gdpll_frd;
% Yin_pll_cal_frd=Yin_cal_frd+Ypll_frd;

%% Calculatoin with current control loop
I_frd = [frd(freqresp(I(1,1),w),w) frd(freqresp(I(1,2),w),w);...
    frd(freqresp(I(2,1),w),w) frd(freqresp(I(2,2),w),w)];
Gdei = [tf([0 0],[0 1]) tf([0 1/Vdcref*Lboost_con*omega],[0 1]);...
    -tf([0 1/Vdcref*Lboost_con*omega],[0 1]) tf([0 0],[0 1])];
Gdei = [tf([0 0],[0 1]) tf([0 0],[0 1]);...
    tf([0 0],[0 1]) tf([0 0],[0 1])];
Gdei_frd = [frd(freqresp(Gdei(1,1),w),w) frd(freqresp(Gdei(1,2),w),w);...
    frd(freqresp(Gdei(2,1),w),w) frd(freqresp(Gdei(2,2),w),w)];
% Gci = -[kpi+kii/s tf([0 0],[0 1]); tf([0 0],[0 1]) kpi+kii/s];
Gci = -[kpi+kii*s/(s^2+omega^2) kii*omega/(s^2+omega^2); -kii*omega/(s^2+omega^2) kpi+kii*s/(s^2+omega^2)];
Gci_frd = [frd(freqresp(Gci(1,1),w),w) frd(freqresp(Gci(1,2),w),w);...
    frd(freqresp(Gci(2,1),w),w) frd(freqresp(Gci(2,2),w),w)];

Ypll_il_frd = (-Gid_cal_frd*Gdel_frd*((-Gdei_frd + Gci_frd)*Gipll_frd*Kv_frd-Gdpll_frd*Kv_frd));
T_frd = (I_frd+Gid_cal_frd*Gdel_frd*(-Gdei_frd + Gci_frd)*Ki_frd);
Yin_il_pll_avg_cal_frd = T_frd\(Ypll_il_frd+Yin_cal_frd);
Zin_il_pll_avg_cal_frd = (Ypll_il_frd+Yin_cal_frd)\T_frd;
%% calculation and simulation comparison

load ZAFE_pll_il_100_1000.mat
Zin_il_pll_sw_sim = Z.ZDQLRAW;

figure(3)
bode(1/Zin_il_pll_avg_sim,Yin_il_pll_avg_cal_frd,Bode_O)
legend('Yin\_il\_avg\_pll\_sim','Yin\_il\_avg\_pll\_cal')
Bode_Darklines(3)

figure(4)
bode(Yin_il_pll_avg_cal_frd,1/Zin_il_pll_sw_sim,Bode_O)
legend('Yin\_il\_PLL\_cal','Yin\_il\_PLL\_sim')
Bode_Darklines(3)

figure(31)
bode(Zin_il_pll_avg_sim,Zin_il_pll_avg_cal_frd,Bode_O)
legend('Zin\_il\_avg\_pll\_sim','Zin\_il\_avg\_pll\_cal')
Bode_Darklines(3)

% 
% figure(4)
% bode(Gvd_avg_sim,Gvd_cal(1,:),Bode_O)
% legend('Gvd\_avg\_sim','Gvd\_cal')
% Bode_Darklines(3)
% 
% figure(5)
% bode(Zin_pll_avg_sim,Zin_pll_cal,Bode_O)
% legend('Zin\_pll\_avg\_sim','Zin\_pll\_cal')
% Bode_Darklines(3)
% Zin_pll_1_avg_sim=Zin_pll_avg_sim;
% Zin_pll_1_cal=Zin_pll_cal;
% save('ZAFE_in_pll_800.mat','Zin_pll_1_avg_sim','Zin_pll_1_cal');

% figure(6)
% bode(Ypll_il,Yin_cal,Yin_il_pll_avg_cal,Bode_O)
% legend('Ypll\_il','Yin','Yin\_il\_pll')
% Bode_Darklines(3)
% 
% figure(7)
% bode(Ypll_il,Yin_cal,Ypll_il+Yin_cal,Bode_O)
% Bode_Darklines(3)
% legend('Ypll\_il','Yin','sum\_Ypll\_il\_Yin')
% 
% figure(8)
% bode(1/Ypll_il,1/Yin_cal,1/(Ypll_il+Yin_cal),Bode_O)
% Bode_Darklines(3)
% legend('Zpll\_il','Zin','p\_Zpll\_il\_Yin')
% 
% figure(9)
% bode(1/Yin_il_avg_cal,1/Yin_cal,Zin_il_pll_avg_cal,Bode_O)
% Bode_Darklines(3)
% legend('Zpll\_il','Zin','p\_Zpll\_il\_Yin')
% 
% figure(10)
% bode(1/Yin_cal,1/Yin_cal_T,1/Ypll_il_T,Zin_il_pll_avg_cal,Bode_O)
% Bode_Darklines(3)
% legend('Zpll\_il','Zin','p\_Ypll\_il\_Yin')