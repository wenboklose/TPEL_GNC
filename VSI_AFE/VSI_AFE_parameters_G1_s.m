% close all;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% AFE parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Vse=57.5;
Vsm=Vse*sqrt(2);
%% boost inductor parameters
Lboost_a    = 474e-6;           % inductor No. 2
RLboost_a   = 0.080;
Lboost_b    = 471e-6;           % inductor No. 1
RLboost_b   = 0.084;
Lboost_c    = 463e-6;           % inductor No. 3
RLboost_c   = 0.110;
%% dc link cap parameters
Cdc_afe     = 98.8e-6; 
RCdc_afe    = 0.049;
LCdc_afe    = 0.65e-6;
%% load resistor
Rdc         = 90;
%% dead time and PWM parameters
DeadTime =0.5e-6;
PWM_Cycle=4000;
%% system parameters
f=400; w=2*pi*f;
fsw=20e+3; Tsw=1/fsw; wsw=2*pi*fsw;
m=0.4;
%parameters for IGBT 6MBP30RH060-50
Vfs=1.0;
Vfd=0.75;
Rsd=0.05;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% VSI parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% output filter inductor parameters
La  = 962e-6;       % inductor No. 1
RLa = 0.123;        
Lb  = 985e-6;       % inductor No. 2
RLb = 0.132;
Lc  = 972e-6;       % inductor No. 3
RLc = 0.089;
%% output filter cap parameters
Ca  = 32e-6;        % cap No. 1
LCa = 0.61e-6;
RCa = 0.051;
Cb  = 31.5e-6;      % cap No. 2
LCb = 0.61e-6;
RCb = 0.054;
Cc  = 31.8e-6;      % cap No. 3
LCc = 0.63e-6;
RCc = 0.052;
%% dc link cap parameters
Cdc_vsi     = 147.3e-6;
LCdc_vsi    = 0.57e-6;
RCdc_vsi    = 0.046;
%% load resistor
R=16;
%% input dc voltage
Vdc=270/1;
%% signal conditioning filter
C1=75e-12;C2=39e-12;R1=15e3;R2=15e3;
tf_filter = tf([0 1],[C1*C2*R1*R2 C2*(R1+R2) 1]);

fprintf('\ndone\n')