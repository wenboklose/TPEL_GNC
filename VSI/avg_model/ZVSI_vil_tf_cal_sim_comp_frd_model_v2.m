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
%% for phase to netural derivation
s=tf([1,0],[0,1]);
I = [tf([0 1],[0 1]) tf([0 0],[0 1]); tf([0 0],[0 1]) tf([0 1],[0 1])];
L=1.0e-3;
RL=0.13;
Vdc=270;
Vod=99.6;
Voq=0;
R=1200;
C=1*31.5e-6;
Cdc=150e-6;
f=400;
omega=2*pi*f;
w=2*pi*f;
fsw=20e+3; Tsw=1/fsw; wsw=2*pi*fsw;

%% PWM delay
Gdeldd = exp((-1/fsw)*s)*cos(omega*(1/fsw));
Gdeldq = exp((-1/fsw)*s)*sin(omega*(1/fsw));
Tdelay = 1.5/fsw;
Gsvmdd = (1-.5*Tdelay*s)/(1+.5*Tdelay*s);%exp(-0.5/fsw*s);
Gsvmdq = -(1-.5*Tdelay*s)*0/(1+.5*Tdelay*s);
Gsvm = tf(JF_DQFromABC(Gsvmdd,omega));
Gsvm = [Gsvmdd 0*Gsvmdq;-0*Gsvmdq Gsvmdd];
Gdel = Gsvm;%I;
%% signal conditioning filter
C1=75e-12;C2=39e-12;R1=18e3;R2=18e3;
tf_filter = tf([0 1],[C1*C2*R1*R2 C2*(R1+R2) 1]);
tf_filter_dq = tf(JF_DQFromABC(tf_filter,omega));
tf_filter_dq = [tf_filter_dq(1,1) 0*tf_filter_dq(1,2); -0*tf_filter_dq(1,2) tf_filter_dq(2,2)];

Ki = tf_filter_dq;
Kv = tf_filter_dq;%Hv;

Vdc=270; Vdcref=270;
Vse=57.5;
Vsm=Vse*sqrt(2);
Vsdq=[sqrt(3/2)*Vsm; 0];
Vdref=Vsdq(1); Vqref=Vsdq(2);
omega = f*2*pi;
Rac =  R;
Vod = Vdref;
Voq = Vqref;
Rl=RL;
Dd = 0.29997;%-(omega*L*Voq-Rac*Vod-Rl*Vod+Rl*omega*C*Voq*Rac+omega^2*L*C*Vod*Rac)/(Rac*Vdc);
Dq = 0.08;%-(-Rl*Voq-Rac*Voq-omega*L*Vod+omega^2*L*C*Voq*Rac-Rl*omega*C*Vod*Rac)/(Rac*Vdc);
Idref = 8.2994;%-(-Vod+omega*C*Voq*Rac)/Rac;
Iqref = 7.8846;%(Voq+omega*C*Vod*Rac)/Rac;

f_v= [80 150 250 200];  % 200 is fix for AFE cases
fi=4000;

for i=1:1:length(f_v)
    fv=f_v(i)
    kpi=2*pi*fi*L/Vdcref*0.5
    kii=2*pi*fi*(Rl)/Vdcref*10
    kpv=2*pi*fv*C
    kiv=2*pi*fv/R%130;%41.8879;%2*pi*fv/R;

    %% get the current loop gain:
    %% Tid model linearization:
    model = 'VSI_v_loop_Z_ac_Kv_Ki_SVM_Ti'; 
    io=getlinio(model);
    %% Linearize the model
    Tid = -linearize(model,5,io);
    

    figHandle=figure(1);
    set(figHandle, 'Position', [50, 10, 1000, 800]);
    bode(Tid,Bode_O)
    Bode_Darklines(3)

    %% Tvd model linearization:
    model = 'VSI_v_loop_Z_ac_Kv_Ki_SVM_Tv';
    % load_system(model);
    % open_system(model);
    io=getlinio(model);
    %% Linearize the model
    Tvd = -linearize(model,25,io);
    figHandle=figure(2);
    set(figHandle,'Position',[50, 10, 1000, 800]);
    bode(Tvd,Bode_O)
    Bode_Darklines(3)
    hold on

    kpi
    kii
    kpv
    kiv

end
% figHandle=figure(3);
% set(figHandle,'Position',[50, 10, 1000, 800]);
% bode(tf_filter,Bode_O)
