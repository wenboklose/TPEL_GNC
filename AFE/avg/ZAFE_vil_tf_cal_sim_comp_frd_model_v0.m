clc
clear all;
close all
%% Bode plot options
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
Bode_O.PhaseWrapping='on';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
linewidth=3; fontsize=18;
linestyle1='b'; linestyle2='g';
linestyle3='r'; linestyle4='k';
linestyle5='m'; linestyle6='c';

linestyle=[linestyle1; linestyle2; linestyle3; linestyle4; linestyle5; linestyle6];

%% power stage parameter initialization
Tstart=0.015;
Rs=1e-1;
Vdcref=270;
alpha=1e10;         %% 0.5
Cdc=105e-6; 
R=93;
RCdc=0.049;
P=Vdcref^2/R;
Vse=57.5;
Vse=Vse+P/3/Vse*Rs*1;
Vsm=Vse*sqrt(2);
Vsm_o=Vse*sqrt(2);
L=500e-6;
RL=80e-3;
I = [1 0; 0 1];
f=400;
omega=2*pi*f;
fsw=20e+3; Tsw=1/fsw; wsw=2*pi*fsw;
Vsdq=[sqrt(3/2)*Vsm; 0];
Vsd_s=Vsdq(1);
Vsq_s=Vsdq(2);
Id_ref=7.6907;
Iq_reff=[-5 -2 0 2 5];
s=tf([1 0],[0 1]);

f_pll = [1 50 200 800];
f_v = [30 35 70 105 135 150];
ki_i = [20];

for i=1:1:length(ki_i)
    Iq_ref = 0*Iq_reff(3);
    %% PLL parameters calculation
    DEF_pll=f_pll(3);                     	%control loop for PLL (Hz)
    DEF_pll_damp=0.707;                     %damping factor for the PLL controller
    DEF_Vin=57.5;                           %input phase to neutral rms voltage
    DEF_Tsw=0.00005;						%switching period
    FWPI_a1=2*DEF_pll_damp*DEF_pll*6.28318530717959*0.57735026918963/DEF_Vin;					%digital PI parameters for frequency loop
    FWPI_a2=-(2*DEF_pll_damp*DEF_pll*6.28318530717959-DEF_pll*DEF_pll*39.47841760435743*DEF_Tsw)*0.57735026918963/DEF_Vin;
    tf_pll_z=tf([FWPI_a1 FWPI_a2],[1 -1],DEF_Tsw);
    tf_pll=d2c(tf_pll_z)
%     kp=2*49*.01;
%     ki=12*pi*2*49*.01*f_pll(i);
%     tf_pll=tf([kp ki],[1 0])
    %% controller parameters calculation (need three sets of parameters)
    fv=70;%f_v(5);%100; 
    fi=1000;
    Lboost_con=500e-6;
    RLboost_con=100e-3;
    Cdc_con=100e-6;
    Rdc_con=90;
%     kpv=2*pi*fv*Cdc_con;%4.4396*0.06;%2*pi*fv*Cdc_con; %0.3770;%
%     kiv=2.5*2*pi*fv/Rdc_con;%4.4396;%2*pi*fv/Rdc_con; %27.9253;%
%     kpi=ki_i(2)*0.00068;%2*pi*fi*Lboost_con/Vdcref; %0.0233;% 
%     kii=ki_i(2);%2*pi*fi*RLboost_con/Vdcref; %4.6542;%
    
    kpv=2*pi*fv*Cdc_con; %0.3770;%
    kiv=2*pi*fv/Rdc_con; %27.9253;%
    kpi=2*pi*fi*Lboost_con/Vdcref; %0.0233;% 
    kii=ki_i(1)*2*pi*fi*RLboost_con/Vdcref; %4.6542;%

    %% signal conditioning filter
    C1=75e-12;C2=39e-12;R1=15e3;R2=15e3;
    tf_filter = tf([0 1],[C1*C2*R1*R2 C2*(R1+R2) 1]);
    tf_filter_dq = tf(JF_DQFromABC(tf_filter,omega));
    tf_filter_dq = [tf_filter_dq(1,1) 0*tf_filter_dq(1,2); -0*tf_filter_dq(2,1) tf_filter_dq(2,2)];

    Gdeldd = exp((-1/fsw)*s)*cos(-omega*(1.5/fsw));
    Gdeldq = exp((-1/fsw)*s)*sin(-omega*(1.5/fsw));
    Tdelay = 1.5/fsw;
    Gsvmdd = (1-.5*Tdelay*s)/(1+.5*Tdelay*s);%exp(-0.5/fsw*s);
    Gsvm = tf(JF_DQFromABC(Gsvmdd,omega));
%     Gsvmdq = -(1-.5*Tdelay*s)*1/(1+.5*Tdelay*s);
%     Gsvm = [Gsvmdd 1*Gsvmdq;-1*Gsvmdq Gsvmdd];
    Gsvm = [Gdeldd 0*Gdeldq;-0*Gdeldq Gdeldd];
    I = [tf([0 1],[0 1]) tf([0 0],[0 1]); tf([0 0],[0 1]) tf([0 1],[0 1])];
    Gdel = Gsvm;%I;
    Ki = tf_filter_dq;
    Kv = tf_filter_dq;
    fprintf('initialization is done!\n')

    %% dq average model simulation for verification
    sys = linearize('AFE_avg_vl_pll_Ki_SVM',0.2);
    %% Zin model linearization
    H=ss(sys);
    Vavg1=[H(1,1) H(1,2);
        H(2,1) H(2,2);];
    Ilavg1=[H(3,1) H(3,2);
        H(4,1) H(4,2);];
    Zin_vl_pll_avg_sim=Vavg1/Ilavg1;

    Dd = Dd_rec(length(Dd_rec)-10);
    Dq = Dq_rec(length(Dq_rec)-10);
    Id = Id_rec(length(Id_rec)-10);
    Iq = 1*Iq_rec(length(Iq_rec)-10);
    Vsd = Vsd_rec(length(Vsd_rec)-10);
    Vsq = Vsq_rec(length(Vsq_rec)-10);
    Vdc = Vdc_rec(length(Vdc_rec)-10);


    %% Calculation based on model of frd
    n= 1e3;
    w = logspace(-2,4,n);
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

    Gdpll=1*[tf([0 0],[0 1]) -Dq*Gpll;tf([0 0],[0 1]) +Dd*Gpll];
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
        tf([0 -1/Vdcref*Lboost_con*omega],[0 1]) tf([0 0],[0 1])];
    Gdei = [tf([0 0],[0 1]) tf([0 0],[0 1]);...
        tf([0 0],[0 1]) tf([0 0],[0 1])];
    Gdei_frd = 1*[frd(freqresp(Gdei(1,1),w),w) frd(freqresp(Gdei(1,2),w),w);...
        frd(freqresp(Gdei(2,1),w),w) frd(freqresp(Gdei(2,2),w),w)];

    % Gci = -[kpi+kii/s tf([0 0],[0 1]); tf([0 0],[0 1]) kpi+kii/s];
    Gci = -[kpi+kii/s 0; 0 kpi+kii/s];
    Gci_frd = [frd(freqresp(Gci(1,1),w),w) frd(freqresp(Gci(1,2),w),w);...
        frd(freqresp(Gci(2,1),w),w) frd(freqresp(Gci(2,2),w),w)];
    Gcv = [kpv+kiv/s tf([0 0],[0 1]); tf([0 0],[0 1]) tf([0 0],[0 1])];
    Gcv_frd = [frd(freqresp(Gcv(1,1),w),w) frd(freqresp(Gcv(1,2),w),w);...
        frd(freqresp(Gcv(2,1),w),w) frd(freqresp(Gcv(2,2),w),w)];
    A = Gid_cal_frd/(I_frd+Gdel_frd*Gci_frd*Gcv_frd*Kv_frd*Gvd_cal_frd);
%     A_cal = Gid_cal/(I+Gdel*Gci*Gcv*Kv*Gvd_cal);
    B = Gdel_frd*(-Gdei_frd+Gci_frd)*Ki_frd;
%     B_cal = Gdel*(-Gdei+Gci)*Ki;
    C = Gdel_frd*(Gdpll_frd*Kv_frd-(-Gdei_frd+Gci_frd)*Gipll_frd*Kv_frd-...
        Gci_frd*Gcv_frd*Kv_frd*Gve_cal_frd);
%     C_cal = Gdel*(Gdpll*Kv-(-Gdei+Gci)*Gipll*Kv-...
%         Gci*Gcv*Kv*Gve_cal);
    C1= Gdel_frd*(-Gci_frd*Gcv_frd*Kv_frd*Gve_cal_frd);
    Yin_vil_pll_cal_frd = (I_frd+A*B)\(A*C+Yin_cal_frd);
    Zin_vil_pll_cal_frd = (A*C+Yin_cal_frd)\(I_frd+A*B);%1/Yin_vil_pll_cal_frd;%
%     Zin_vil_pll_cal = (A_cal*C_cal+Yin_cal)\(I+A_cal*B_cal);%1/Yin_vil_pll_cal_frd;%
    Zin_vil_cal_frd = (A*C1+Yin_cal_frd)\(I_frd+A*B);
    
    Ypll_il_frd = (-Gid_cal_frd*(Gdel_frd*(-Gdei_frd + Gci_frd)*Gipll_frd-Gdpll_frd));
    T_frd = (I_frd+Gid_cal_frd*Gdel_frd*(-Gdei_frd + Gci_frd)*Ki_frd);
    Yin_il_pll_avg_cal_frd = T_frd\(Ypll_il_frd+Yin_cal_frd);
    Zin_il_pll_avg_cal_frd = (Ypll_il_frd+Yin_cal_frd)\T_frd;
    
    Zqq_pll = Vsd^2/P+Id*Gpll;
  
    figure(6)
    bode(Zin_vl_pll_avg_sim,Zin_vil_pll_cal_frd,Bode_O)
    hold on
    Bode_Darklines(3)
    
    figure(7)
    bode(tf_pll,Bode_O)
    hold on
    Bode_Darklines(3)
    
    figure(8)
    bode(1/Zin_vl_pll_avg_sim,Yin_il_pll_avg_cal_frd,Bode_O)
    hold on
    Bode_Darklines(3)
end
save_file='Z_AFE_fi_1k_fv_70_pll_1_avg.mat';
save(save_file,'Zin_vl_pll_avg_sim')