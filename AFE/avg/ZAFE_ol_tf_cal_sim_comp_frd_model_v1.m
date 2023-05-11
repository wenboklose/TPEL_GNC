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
Bode_O.XLim={[40 1e4]};
Bode_O.XLimMode={'manual'};
Bode_O.FreqUnits='Hz';
Bode_O.PhaseWrapping='off';

% bdclose all
% set_param(0,'CharacterEncoding','ISO-8859-1')
Tstart=0.015;
Rs=0e-6;
Vdcref=270;
alpha=1e10;         %% 0.5
Cdc=100e-6; 
R=90;
RCdc=0.070;
P=Vdcref^2/R;
Vse=57.5;
Vse=Vse+P/3/Vse*Rs*1;
Vsm=Vse*sqrt(2);
Vsm_o=Vse*sqrt(2);
L=500e-6; 
RL=500e-3;

f=400; w=2*pi*f;
omega=2*pi*f;
fsw=20e+3; Tsw=1/fsw; wsw=2*pi*fsw;

f_pll = [10 100 200];

for i=1:1:1
    %% PLL
    DEF_pll=f_pll(2);							%control loop for PLL (Hz)
    DEF_pll_damp=0.707;					%damping factor for the PLL controller
    DEF_Vin=57.5;						%input phase to neutral rms voltage
    DEF_Tsw=0.00005;						%switching period
    FWPI_a1=2*DEF_pll_damp*DEF_pll*6.28318530717959*0.57735026918963/DEF_Vin;					%digital PI parameters for frequency loop
    FWPI_a2=-(2*DEF_pll_damp*DEF_pll*6.28318530717959-DEF_pll*DEF_pll*39.47841760435743*DEF_Tsw)*0.57735026918963/DEF_Vin;
    tf_pll_z=tf([FWPI_a1 FWPI_a2],[1 -1],DEF_Tsw);
    tf_pll=d2c(tf_pll_z);

    Vsdq=[sqrt(3/2)*Vsm; 0];
    Vsd=Vsdq(1);
    Vsq=Vsdq(2);
    Dd=0.3657;%72.4/270;%
    Dq=-0.0382;%-5.2/270;%
    Id=8.2109;
    Iq=0.0227;
    fprintf('initialization is done!\n')

    I = [tf([0 1],[0 1]) tf([0 0],[0 1]); tf([0 0],[0 1]) tf([0 1],[0 1])];
    s=tf([1 0],[0 1]);
    %% on board signal condition filter
    %% signal conditioning filter
    C1=75e-12;C2=39e-12;R1=15e3;R2=15e3;
    tf_filter = tf([0 1],[C1*C2*R1*R2 C2*(R1+R2) 1]);
    tf_filter_dq = tf(JF_DQFromABC(tf_filter,omega));
    tf_filter_dq = [tf_filter_dq(1,1) tf([0 0],[0 1]); tf([0 0],[0 1]) tf_filter_dq(2,2)];
    %% PWM delay
    Gdeldd = exp((-1/fsw)*s)*cos(omega*(1/fsw));
    Gdeldq = exp((-1/fsw)*s)*sin(omega*(1/fsw));
    Tdelay = 1.5/fsw;
    Gsvmdd = (1-.5*Tdelay*s)/(1+.5*Tdelay*s);%exp(-0.5/fsw*s);
    Gsvmdq = -(1-.5*Tdelay*s)*0/(1+.5*Tdelay*s);
    Gsvm = [Gsvmdd 0*Gsvmdq;-0*Gsvmdq Gsvmdd];
    Gdel = I;%Gsvm;%I;

    %% Zin_pll model linearization:
    model = 'AFE_avg_Zol_pll';

    sys = linearize('AFE_avg_Zol_pll',0.2);
    %% Zin model linearization
    H=ss(sys);
    Vavg1=[H(1,1) H(1,2);
        H(2,1) H(2,2);];
    Ilavg1=[H(3,1) H(3,2);
        H(4,1) H(4,2);];
    Zin_pll_avg_sim=Vavg1/Ilavg1;

    Dd = Dd_rec(length(Dd_rec)-10);
    Dq = Dq_rec(length(Dq_rec)-10);
    Id = Id_rec(length(Id_rec)-10);
    Iq = 1*Iq_rec(length(Iq_rec)-10);
    Vsd = Vsd_rec(length(Vsd_rec)-10);
    Vsq = Vsq_rec(length(Vsq_rec)-10);
    Vdc = Vdc_rec(length(Vdc_rec)-10);

    %% Calculation based on model of frd
    w = 1:1:1e4;
    w = 2*pi*w;

    E0=Vsd;
    Gpll=tf_pll/(s+E0*tf_pll);
    Gipll=[tf([0 0],[0 1]) Iq*Gpll;tf([0 0],[0 1]) -Id*Gpll];

    Gdel_frd = [frd(freqresp(Gdel(1,1),w),w) frd(freqresp(Gdel(1,2),w),w);...
        frd(freqresp(Gdel(2,1),w),w) frd(freqresp(Gdel(2,2),w),w)];

    tf_filter_dq_frd = [frd(freqresp(tf_filter_dq(1,1),w),w) frd(freqresp(tf_filter_dq(1,2),w),w);...
        frd(freqresp(tf_filter_dq(2,1),w),w) frd(freqresp(tf_filter_dq(2,2),w),w)];

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

    Ypll_frd=Gid_cal_frd*Gdel_frd*Gdpll_frd*tf_filter_dq_frd;
    Yin_pll_cal_frd=Yin_cal_frd+Ypll_frd;

%     figure(5)
%     bode(Yin_pll_cal_frd,Bode_O)
%     legend('Yin\_pll\_cal')
%     Bode_Darklines(3)
%     hold on

    figure(6)
    bode(Zin_pll_avg_sim,Bode_O)
%     legend('Zin\_ol\_PLL\_sim')
    Bode_Darklines(3)
    hold on

%     figure(7)
%     bode(Ypll_frd,Bode_O)
%     Bode_Darklines(3)
%     hold on

end
bode(1/Yin_pll_cal_frd,Bode_O)
Bode_Darklines(3)
figure(7)
bode(Yin_pll_cal_frd,Bode_O)
Bode_Darklines(3)
