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
RL=0.15;
Vdc=270;
Vod=99.6;
Voq=0;
R=15;
C=31.3e-6;
Cdc=150e-6;
f=400;
omega=2*pi*f;
w=2*pi*f;
fsw=20e+3; Tsw=1/fsw; wsw=2*pi*fsw;

%% on board signal condition filter
%% signal conditioning filter
C1=75e-12;C2=39e-12;R1=15e3;R2=15e3;
tf_filter = tf([0 1],[C1*C2*R1*R2 C2*(R1+R2) 1]);
tf_filter_dq = tf(JF_DQFromABC(tf_filter,omega));

%% PWM delay
Gdeldd = exp((-1/fsw)*s)*cos(omega*(1/fsw));
Gdeldq = exp((-1/fsw)*s)*sin(omega*(1/fsw));
Tdelay = 1.5/fsw;
Gsvmdd = (1-.5*Tdelay*s)/(1+.5*Tdelay*s);%exp(-0.5/fsw*s);
Gsvmdq = -(1-.5*Tdelay*s)*0/(1+.5*Tdelay*s);
Gsvm = [Gsvmdd 1*Gsvmdq;-1*Gsvmdq Gsvmdd];
Gdel = Gsvm;%I;
%% signal conditioning filter
C1=75e-12;C2=39e-12;R1=15e3;R2=15e3;
tf_filter = tf([0 1],[C1*C2*R1*R2 C2*(R1+R2) 1]);
tf_filter_dq = tf(JF_DQFromABC(tf_filter,omega));
tf_filter_dq = [tf_filter_dq(1,1) 1*tf_filter_dq(1,2); -1*tf_filter_dq(1,2) tf_filter_dq(2,2)];

Ki = tf_filter_dq;
Kv = tf_filter_dq;%Hv;

Vdc=270; Vdcref=270;
Vse=57.5;
Vsm=Vse*sqrt(2);
Vsdq=[sqrt(3/2)*Vsm; 0];
Vdref=Vsdq(1); Vqref=Vsdq(2);
omega = 400*2*pi;
Rac =  R;
Vod = Vdref;
Voq = Vqref;
Rl=RL;
Dd = -(omega*L*Voq-Rac*Vod-Rl*Vod+Rl*omega*C*Voq*Rac+omega^2*L*C*Vod*Rac)/(Rac*Vdc);
Dq = -(-Rl*Voq-Rac*Voq-omega*L*Vod+omega^2*L*C*Voq*Rac-Rl*omega*C*Vod*Rac)/(Rac*Vdc);
Idref = -(-Vod+omega*C*Voq*Rac)/Rac
Iqref = (Voq+omega*C*Vod*Rac)/Rac

fv=100; fi=1000;

kpi=(46.8*8.3e-5)/Vdcref;%0.0116;%2*pi*fi*L/Vdcref*0.5;
kii=46.8/Vdcref;%34.9066;%2*pi*fi*(Rl)/Vdcref*10;
kpv=130*2.2e-5;%0.0066;%2*pi*fv*C;
kiv=130;%41.8879;%2*pi*fv/R;

%% Zo_VSI_sim model linearization:
model = 'VSI_i_loop_Z_ac';

sys = linearize('VSI_i_loop_Z_ac_Ki_SVM',2);
%% Zin model linearization
H=ss(sys);
Vavg1=[H(1,1) H(1,2);
    H(2,1) H(2,2);];

Ilavg1=[H(3,1) H(3,2);
    H(4,1) H(4,2);];

Isavg1=[H(5,1) H(5,2);
    H(6,1) H(6,2);];
Zo_il_avg_sim=Vavg1/Isavg1;
% figure(1)
% bode(Zo_il_avg_sim,Bode_O)
% legend('Zo\_avg\_sim','Zo\_cal')
% Bode_Darklines(3)

Dd = Dd(length(Dd)-10);
Dq = Dq(length(Dq)-10);

% %% Zo_VSI_sim model linearization:
% model = 'VSI_o_loop_Z_ac';
% 
% sys = linearize('VSI_o_loop_Z_ac',0.5);
% %% Zin model linearization
% H=ss(sys);
% Vavg1=[H(1,1) H(1,2);
%     H(2,1) H(2,2);];
% 
% Ilavg1=[H(3,1) H(3,2);
%     H(4,1) H(4,2);];
% 
% Isavg1=[H(5,1) H(5,2);
%     H(6,1) H(6,2);];
% Zin_avg_sim=Vavg1/Isavg1;

% % w=2*pi*f;
% %% Gii_VSI_sim model linearization:
% model = 'VSI_o_loop_Gii';
% 
% sys = linearize('VSI_o_loop_Gii',0.5);
% %% Gii_VSI model linearization
% H=ss(sys);
% Vavg1=[H(1,1) H(1,2);
%     H(2,1) H(2,2);];
% 
% Ilavg1=[H(3,1) H(3,2);
%     H(4,1) H(4,2);];
% 
% Isavg1=[H(5,1) H(5,2);
%     H(6,1) H(6,2);];
% Gii_avg_sim=Vavg1/Isavg1;
% %% Gid model linearization:
% model = 'VSI_o_loop_Gid';
% 
% %% Create the linearization I/O as specified in AFE_avg_Gid
% ios(4) = linio('VSI_o_loop_Gid/ilq',1,'out');
% ios(3) = linio('VSI_o_loop_Gid/ild',1,'out');
% ios(2) = linio('VSI_o_loop_Gid/Dq',1,'in');
% ios(1) = linio('VSI_o_loop_Gid/Dd',1,'in');
% 
% %% Linearize the model
% Gid_avg_sim = linearize(model,0.5,ios);
% 
% %% Gvd model linearization:
% model = 'VSI_o_loop_Gvd';
% 
% %% Create the linearization I/O as specified in AFE_avg_Gid
% ios(4) = linio('VSI_o_loop_Gvd/vq',1,'out');
% ios(3) = linio('VSI_o_loop_Gvd/vd',1,'out');
% ios(2) = linio('VSI_o_loop_Gvd/Dq',1,'in');
% ios(1) = linio('VSI_o_loop_Gvd/Dd',1,'in');
% 
% %% Linearize the model
% Gvd_avg_sim = linearize(model,0.5,ios);
% 
% %% Gig_VSI_sim model linearization:
% sys = linearize('VSI_o_loop_Gig',0.5);
% %% Gig_VSI model linearization
% H=ss(sys);
% 
% Vavg=H(1,1);
% 
% Ilavg=H(2,1);
% 
% Isavg=H(3,1);
% Zsavg=Isavg/Vavg; Zlavg=Ilavg/Vavg;
% 
% Gig_avg_sim = [Zlavg 0;Zsavg 0];
% 
% %% Gig_VSI_sim model linearization:
% sys = linearize('VSI_o_loop_Gvg',0.5);
% %% Gii_VSI model linearization
% H=ss(sys);
% 
% Vavg=H(1,1);
% 
% Ilavg=H(2,1);
% 
% Isavg=H(3,1);
% Zsavg=Isavg/Vavg; Zlavg=Ilavg/Vavg;
% 
% Gvg_avg_sim = [Zlavg 0;Zsavg 0];
% 
% %% Yin_VSI_sim model linearization:
% sys = linearize('VSI_o_loop_Z_dc',0.5);
% %% Gii_VSI model linearization
% H=ss(sys);
% 
% Vavg=H(1,1);
% 
% Ilavg=H(2,1);
% 
% Isavg=H(3,1);
% Yin_avg_sim=Ilavg/Vavg;

%% Zo_VSI_cal
C=C/3;
GL=3*L*s+3*Rl;
den1=(-3^2*omega^2*L^2*C*s-(GL)^2*C*s-GL)^2+(3^2*omega^2*L^2*omega*C-3*omega*L+(GL)^2*omega*C)^2;
num11=(3^2*omega^2*L^2+(GL)^2)*(-(3^2*omega^2*L^2*C*s+(GL)^2*C*s+GL));
num12=-(3^2*omega^2*L^2+(GL)^2)*(3^2*omega^2*L^2*omega*C-3*omega*L+(GL)^2*omega*C);
num21=-num12;
num22=num11;
Zo_cal=[num11/den1 num12/den1;num21/den1 num22/den1];
w = 1:2:1e4;
w = 2*pi*w;
Zo_cal_frd = -1/3*[frd(freqresp(Zo_cal(1,1),w),w) frd(freqresp(Zo_cal(1,2),w),w);...
    frd(freqresp(Zo_cal(2,1),w),w) frd(freqresp(Zo_cal(2,2),w),w)];

%% Gii_VSI_cal
% C=C/3;
% GL=3*L*s+3*Rl;
den1=(1-omega*C*3*omega*L+C*s*GL)^2+(omega*C*GL+C*s*3*omega*L)^2;
num11=1-omega*C*3*omega*L+C*s*GL;
num12=omega*C*GL+C*s*3*omega*L;
num21=-num12;
num22=num11;
Gii_cal=[num11/den1 num12/den1;num21/den1 num22/den1];
Gii_cal_frd = -[frd(freqresp(Gii_cal(1,1),w),w) frd(freqresp(Gii_cal(1,2),w),w);...
    frd(freqresp(Gii_cal(2,1),w),w) frd(freqresp(Gii_cal(2,2),w),w)];

%% Gid_VSI_cal
% C=C/3;
% GL=3*L*s+3*Rl;
den1=[((GL*C*s+1)*(C^2*s^2+omega^2*C^2)-omega^2*C^2)^2+(-3*omega*L*C*s*(C^2*s^2+omega^2*C^2)+omega*C*C*s)^2]/[C*s*(C^2*s^2+omega^2*C^2)*Vdc];
num11=(GL*C*s+1)*(C^2*s^2+omega^2*C^2)-omega^2*C^2;
num12=-(-3*omega*L*C*s*(C^2*s^2+omega^2*C^2)+omega*C*C*s);
num21=-num12;
num22=num11;
Gid_cal=[num11/den1 num12/den1;num21/den1 num22/den1];
Gid_cal_frd = 3*[frd(freqresp(Gid_cal(1,1),w),w) frd(freqresp(Gid_cal(1,2),w),w);...
    frd(freqresp(Gid_cal(2,1),w),w) frd(freqresp(Gid_cal(2,2),w),w)];

%% Gvd_VSI_cal
% C=C/3;
% GL=3*L*s+3*Rl;
den1=[(GL*C*s-3*omega*L*omega*C+1)^2+(GL*omega*C+3*omega*L*C*s)^2]/Vdc;
num11=GL*C*s-3*omega*L*omega*C+1;
num12=(GL*omega*C+3*omega*L*C*s);
num21=-num12;
num22=num11;
Gvd_cal=[num11/den1 num12/den1;num21/den1 num22/den1];
Gvd_cal_frd = 1*[frd(freqresp(Gvd_cal(1,1),w),w) frd(freqresp(Gvd_cal(1,2),w),w);...
    frd(freqresp(Gvd_cal(2,1),w),w) frd(freqresp(Gvd_cal(2,2),w),w)];

%% Gig_VSI_cal
% C=C/3;
% GL=3*L*s+3*Rl;
den1=((GL*C*s+1)*(C^2*s^2+omega^2*C^2)-omega^2*C^2)^2+(3*omega*L*C*s*(C^2*s^2+omega^2*C^2)-omega*C^2*s)^2;
num11=((GL*C*s+1)*(C^2*s^2+omega^2*C^2)-omega^2*C^2)*Dd*C*s*(C^2*s^2+omega^2*C^2)+(3*omega*L*C*s*(C^2*s^2+omega^2*C^2)-omega*C^2*s)*Dq*C*s*(C^2*s^2+omega^2*C^2);
num21=((GL*C*s+1)*(C^2*s^2+omega^2*C^2)-omega^2*C^2)*Dq*C*s*(C^2*s^2+omega^2*C^2)-(3*omega*L*C*s*(C^2*s^2+omega^2*C^2)-omega*C^2*s)*Dd*C*s*(C^2*s^2+omega^2*C^2);
num12=0;
num22=0;

Gig_cal=[num11/den1 num12/den1;num21/den1 num22/den1];
Gig_cal_frd = 3*[frd(freqresp(Gig_cal(1,1),w),w) frd(freqresp(Gig_cal(1,2),w),w);...
    frd(freqresp(Gig_cal(2,1),w),w) frd(freqresp(Gig_cal(2,2),w),w)];

%% Gvg_VSI_cal
% C=C/3;
% GL=3*L*s+3*Rl;
den1=(GL*C*s-3*omega^2*L*C+1)^2+(GL*omega*C+3*omega*L*C*s)^2;
num11=(GL*C*s-3*omega^2*L*C+1)*Dd+(GL*omega*C+3*omega*L*C*s)*Dq;
num21=-(GL*omega*C+3*omega*L*C*s)*Dd+(GL*C*s-3*omega^2*L*C+1)*Dq;
num12=0;
num22=0;

Gvg_cal=[num11/den1 num12/den1;num21/den1 num22/den1];
Gvg_cal_frd = [frd(freqresp(Gvg_cal(1,1),w),w) frd(freqresp(Gvg_cal(1,2),w),w);...
    frd(freqresp(Gvg_cal(2,1),w),w) frd(freqresp(Gvg_cal(2,2),w),w)];

%% Yin_VSI_cal
% C=C/3;
% GL=3*L*s+3*Rl;
den1=((GL*C*s+1)*(C^2*s^2+omega^2*C^2)-omega^2*C^2)^2+(3*omega*L*C*s*(C^2*s^2+omega^2*C^2)-omega*C^2*s)^2;
num11=((GL*C*s+1)*(C^2*s^2+omega^2*C^2)-omega^2*C^2)*Dd*C*s*(C^2*s^2+omega^2*C^2)+(3*omega*L*C*s*(C^2*s^2+omega^2*C^2)-omega*C^2*s)*Dq*C*s*(C^2*s^2+omega^2*C^2);
num21=((GL*C*s+1)*(C^2*s^2+omega^2*C^2)-omega^2*C^2)*Dq*C*s*(C^2*s^2+omega^2*C^2)-(3*omega*L*C*s*(C^2*s^2+omega^2*C^2)-omega*C^2*s)*Dd*C*s*(C^2*s^2+omega^2*C^2);
num12=0;
num22=0;

Yin_cal = 3*(Dd*num11/den1 + Dq*num21/den1)+1/(1/(Cdc*s)+0.01);
Yin_cal_frd = frd(freqresp(Yin_cal,w),w);

%% Calculatoin with current control loop


Gdel_frd = [frd(freqresp(Gdel(1,1),w),w) frd(freqresp(Gdel(1,2),w),w);...
    frd(freqresp(Gdel(2,1),w),w) frd(freqresp(Gdel(2,2),w),w)];

Ki_frd = [frd(freqresp(Ki(1,1),w),w) frd(freqresp(Ki(1,2),w),w);...
    frd(freqresp(Ki(2,1),w),w) frd(freqresp(Ki(2,2),w),w)];

Kv_frd = [frd(freqresp(Kv(1,1),w),w) frd(freqresp(Kv(1,2),w),w);...
    frd(freqresp(Kv(2,1),w),w) frd(freqresp(Kv(2,2),w),w)];


I_frd = [frd(freqresp(I(1,1),w),w) frd(freqresp(I(1,2),w),w);...
    frd(freqresp(I(2,1),w),w) frd(freqresp(I(2,2),w),w)];
Gdei = [tf([0 0],[0 1]) tf([0 1/Vdcref*L*omega],[0 1]);...
    -tf([0 1/Vdcref*L*omega],[0 1]) tf([0 0],[0 1])];
% Gdei = [tf([0 0],[0 1]) tf([0 0],[0 1]);...
%     tf([0 0],[0 1]) tf([0 0],[0 1])];
Gdei_frd = [frd(freqresp(Gdei(1,1),w),w) frd(freqresp(Gdei(1,2),w),w);...
    frd(freqresp(Gdei(2,1),w),w) frd(freqresp(Gdei(2,2),w),w)];
Gci = [kpi+kii/s tf([0 0],[0 1]); tf([0 0],[0 1]) kpi+kii/s];
Gci_frd = [frd(freqresp(Gci(1,1),w),w) frd(freqresp(Gci(1,2),w),w);...
    frd(freqresp(Gci(2,1),w),w) frd(freqresp(Gci(2,2),w),w)];

%% calculation of VSI output impedance with current loop control
Zo_il_cal_frd=Gvd_cal_frd*((I_frd+Gdel_frd*(-Gdei_frd+Gci_frd)*Ki_frd*Gid_cal_frd)\Gdel_frd*(Gdei_frd-Gci_frd)*Ki_frd*Gii_cal_frd)+Zo_cal_frd;
% calculation and simulation comparison
figure(1)
bode(Zo_il_avg_sim,Zo_il_cal_frd,Zo_cal_frd,Bode_O)
legend('Zo\_il\_avg\_sim','Zo\_il\_cal','Zo\_cal')
Bode_Darklines(3)

% figure(2)
% bode(-Gii_avg_sim,Gii_cal_frd,Bode_O)
% legend('Gii\_avg\_sim','Gii\_cal')
% Bode_Darklines(3)
% 
% figure(3)
% bode(Gid_avg_sim,Gid_cal_frd,Bode_O)
% legend('Gid\_avg\_sim','Gid\_cal')
% Bode_Darklines(3)
% 
% figure(4)
% bode(Gvd_avg_sim,Gvd_cal_frd/3,Bode_O)
% legend('Gvd\_avg\_sim','Gvd\_cal')
% Bode_Darklines(3)
% 
% figure(5)
% bode(Gig_avg_sim,Gig_cal_frd,Bode_O)
% legend('Gig\_avg\_sim','Gig\_cal')
% Bode_Darklines(3)
% 
% figure(6)
% bode(Gvg_avg_sim,Gvg_cal_frd,Bode_O)
% legend('Gvg\_avg\_sim','Gvg\_cal')
% Bode_Darklines(3)
% 
% figure(7)
% bode(1/Yin_avg_sim,1/Yin_cal_frd,Bode_O)
% legend('Yin\_avg\_sim','Yin\_cal')
% Bode_Darklines(3)
