clc
clear all;
close all
%% Bode plot options
Bode_O=bodeoptions;
Bode_O.XLabel.FontSize=18;
Bode_O.YLabel.FontSize=18;
Bode_O.TickLabel.FontSize=18;
Bode_O.Title.FontSize=18;
Bode_O.title.String=' ';
Bode_O.Grid='on';
Bode_O.XLim={[1 1e4]};
Bode_O.XLimMode={'manual'};
Bode_O.FreqUnits='Hz';
Bode_O.PhaseWrapping='off';

%% power stage parameter initialization
Tstart=0.015;
Rs=1e-1;
Vdcref=270;
alpha=1e10;         %% 0.5
Cdc=105e-6; 
R=96;
RCdc=0.049;
P=Vdcref^2/R;
Vse=57.5;
Vse=Vse+P/3/Vse*Rs*1;
Vsm=Vse*sqrt(2);
Vsm_o=Vse*sqrt(2);
L=470e-6;
RL=110e-3;
I = [1 0; 0 1];
f=400;
omega=2*pi*f;
fsw=20e+3; Tsw=1/fsw; wsw=2*pi*fsw;
Vsdq=[sqrt(3/2)*Vsm; 0];
Vsd_s=Vsdq(1);
Vsq_s=Vsdq(2);
Id_ref=7.6307;
Iq_ref=0;
s=tf([1 0],[0 1]);

f_pll = [50 100 200];
f_v = [50 153 100 220 260 300]; % 70 is fixed for VSI cases
ki_i = [2.5 5 20];

for i=1:1:3%length(f_v)
    
    %% PLL parameters calculation
DEF_pll=f_pll(2);                     		%control loop for PLL (Hz)
DEF_pll_damp=0.707;                     %damping factor for the PLL controller
DEF_Vin=57.5;                           %input phase to neutral rms voltage
DEF_Tsw=0.00005;						%switching period
FWPI_a1=2*DEF_pll_damp*DEF_pll*6.28318530717959*0.57735026918963/DEF_Vin;					%digital PI parameters for frequency loop
FWPI_a2=-(2*DEF_pll_damp*DEF_pll*6.28318530717959-DEF_pll*DEF_pll*39.47841760435743*DEF_Tsw)*0.57735026918963/DEF_Vin;
tf_pll_z=tf([FWPI_a1 FWPI_a2],[1 -1],DEF_Tsw);
tf_pll=d2c(tf_pll_z);

%% controller parameters calculation (need three sets of parameters)
fv=f_v(6); fi=1000;
Lboost_con=500e-6;
RLboost_con=100e-3;
Cdc_con=105e-6;
Rdc_con=96;
% kpv=2*pi*fv*Cdc_con;%4.4396*0.06;%2*pi*fv*Cdc_con; %0.3770;%
% kiv=2.5*2*pi*fv/Rdc_con;%4.4396;%2*pi*fv/Rdc_con; %27.9253;%
% kpi=ki_i(2)*0.00068;%2*pi*fi*Lboost_con/Vdcref; %0.0233;% 
% kii=ki_i(2);%2*pi*fi*RLboost_con/Vdcref; %4.6542;%

kpv=2*pi*fv*Cdc_con;%4.4396*0.06;%2*pi*fv*Cdc_con; %0.3770;%
kiv=2*pi*fv/Rdc_con;%4.4396;%2*pi*fv/Rdc_con; %27.9253;%
kpi=2*pi*fi*Lboost_con/Vdcref; %0.0233;% 
kii=ki_i(i)*2*pi*fi*RLboost_con/Vdcref; %4.6542;%


%% signal conditioning filter
C1=75e-12;C2=39e-12;R1=15e3;R2=15e3;
tf_filter = tf([0 1],[C1*C2*R1*R2 C2*(R1+R2) 1]);
tf_filter_dq = tf(JF_DQFromABC(tf_filter,omega));
tf_filter_dq = [tf_filter_dq(1,1) 0*tf_filter_dq(1,2); 0*tf_filter_dq(2,1) tf_filter_dq(2,2)];

Gdeldd = exp((-1/fsw)*s)*cos(-omega*(1.5/fsw));
Gdeldq = exp((-1/fsw)*s)*sin(-omega*(1.5/fsw));
Tdelay = 1.5/fsw;
Gsvmdd = (1-.5*Tdelay*s)/(1+.5*Tdelay*s);%exp(-0.5/fsw*s);
Gsvmdq = -(1-.5*Tdelay*s)*1/(1+.5*Tdelay*s);
Gsvm = [Gsvmdd 0*Gsvmdq;-0*Gsvmdq Gsvmdd];
% Gsvm = tf(JF_DQFromABC(Gsvmdd,omega));
I = [tf([0 1],[0 1]) tf([0 0],[0 1]); tf([0 0],[0 1]) tf([0 1],[0 1])];
Gdel = Gsvm;%I;
Ki = tf_filter_dq;
Kv = tf_filter_dq;
fprintf('initialization is done!\n')

%% get the current loop gain:
%% Tid model linearization:
model = 'AFE_avg_il_pll_Ki_SVM_LA_Ti';
% load_system(model);
% open_system(model);
io=getlinio(model);
%% Linearize the model
Tid = -linearize(model,0.3,io);
%% tfvdcidref model linearization:
% model = 'AFE_avg_il_pll_Ki_SVM_LA_tfcli';
% % load_system(model);
% % open_system(model);
% io=getlinio(model);
% %% Linearize the model
% tfvdcidref = -linearize(model,0.3,io);

figHandle=figure(1);
set(figHandle, 'Position', [50, 10, 1000, 800]);
bode(Tid,Bode_O)
Bode_Darklines(3)
hold on
%% Tvdc model linearization:
model = 'AFE_avg_vl_pll_Ki_SVM_LA_Tv';
% load_system(model);
% open_system(model);
io=getlinio(model);
%% Linearize the model
Tvdc = -linearize(model,0.6,io);
figHandle=figure(2);
set(figHandle,'Position',[50, 10, 1000, 800]);
bode(Tvdc,Bode_O)
Bode_Darklines(3)
hold on

Ti(i)=tf(Tid);
kpv
kiv
kpi
kii
end