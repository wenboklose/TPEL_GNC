clc
clear all;
% close all

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

nsamples=1000;
K_norm=1;L_norm=1e-3;C_norm=32e-6;C1_norm=75e-12;C2_norm=39e-12;R1_norm=15e3;R2_norm=15e3;Vdc_norm=270;R_norm=15;R_rc_norm=1e3;C_rc_norm=51e-12;
Nom=[L_norm C_norm K_norm C1_norm C2_norm R1_norm R2_norm Vdc_norm R_norm R_rc_norm C_rc_norm];
Nk=nsamples;
Nc=size(Nom,2);
randn('state',sum(100*clock));
%
Tn=zeros(Nc,Nk); % Reserve space; decreases run time
%
% Create tolerance array T;
% 
Tr=0.1;Td=0.1;Te=0.05;Tf=0.001;Tg=0.1;
T=[-Tr*1 -Te*1 -Td*1 -Te*1 -Te*1 -Tf*1 -Tf*1 -Tg*1 -Tg*1 -Tf -Te;Tr*1 Te*1 Td*1 Te*1 Te*1 Tf*1 Tf*1 Tg*1 Tg*1 Tf Te];
% T below for fig 20 and fig 21
%T=[ -Tr -Tr -Tr -Tr -Te -Te;0.03 Tr 0.03 Tr Te Te];
% Convert to tolerance multipliers
for k=1:Nk
   for w=1:Nc
      Tn(w,k)=Nom(w)*(((T(2,w)-T(1,w))/6)*(randn+3)+T(1,w)+1);
   end % end w loop
end % end k loop
L_rec=Tn(1,:);
C_rec=Tn(2,:);
K_rec=Tn(3,:);
C1_rec=Tn(4,:);
C2_rec=Tn(5,:);
R1_rec=Tn(6,:);
R2_rec=Tn(7,:);
Vdc_rec=Tn(8,:);
R_rec=Tn(9,:);
R_rc_rec=Tn(10,:);
C_rc_rec=Tn(11,:);

for i=1:1:nsamples
    
    %% for phase to netural derivation
        s=tf([1,0],[0,1]);
        I = [tf([0 1],[0 1]) tf([0 0],[0 1]); tf([0 0],[0 1]) tf([0 1],[0 1])];
        L=L_rec(i);%1.0e-3;
        RL=0.13;
        Vdc=Vdc_rec(i);%270; 
        Vdcref=270;
        Vse=57.5;
        R=R_rec(i);%15;
        C=C_rec(i);%31.5e-6;
        Cdc=150e-6;
        f=400;
        omega=2*pi*f;
        w=2*pi*f;
        fsw=20e+3; Tsw=1/fsw; wsw=2*pi*fsw;
        %% PWM delay
        Gdeldd = exp((-1/fsw)*s)*cos(omega*(1/fsw));
        Gdeldq = exp((-1/fsw)*s)*sin(omega*(1/fsw));
        Tdelay = 1/fsw;
        Gsvmdd = (1-.5*Tdelay*s)/(1+.5*Tdelay*s);%exp(-0.5/fsw*s);
        Gsvmdq = -(1-.5*Tdelay*s)*0/(1+.5*Tdelay*s);
        Gsvm = [Gsvmdd 1*Gsvmdq;-1*Gsvmdq Gsvmdd];
        Gdel = Gsvm;%I;
        %% signal conditioning filter
        C1=C1_rec(i);C2=C2_rec(i);R1=R1_rec(i);R2=R2_rec(i);
        K=K_rec(i);
        R_rc=R_rc_rec(i);C_rc=C_rc_rec(i);
%         C1=75e-12;C2=39e-12;R1=15e3;R2=15e3;
        tf_filter = K*tf([0 1],[C1*C2*R1*R2 C2*(R1+R2) 1])*tf([0 1],[R_rc*C_rc 1]);
        tf_filter_dq_s = tf(JF_DQFromABC(tf_filter,omega));
        tf_filter_dq = [tf_filter_dq_s(1,1) 1*tf_filter_dq_s(1,2); -1*tf_filter_dq_s(1,2) tf_filter_dq_s(2,2)];

        Ki = tf_filter_dq;
        Kv = tf_filter_dq;%Hv;

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
        Idref = -(-Vod+omega*C*Voq*Rac)/Rac;
        Iqref = (Voq+omega*C*Vod*Rac)/Rac;

        fv=180; fi=1000;

        kpi=(46.8*8.3e-5)/0.5;%0.0116;%2*pi*fi*L/Vdcref*0.5;
        kii=46.8/0.5;%34.9066;%2*pi*fi*(Rl)/Vdcref*10;
        kpv=2*pi*fv*29.925e-6;%C;%130*2.2e-5;%0.0066;%2*pi*fv*C;
        kiv=2*pi*fv/15;%130;%41.8879;%2*pi*fv/R;

        %% Zo_VSI_cal
        C=C/3;
        GL=3*L*s+3*Rl;
        den1=(-3^2*omega^2*L^2*C*s-(GL)^2*C*s-GL)^2+(3^2*omega^2*L^2*omega*C-3*omega*L+(GL)^2*omega*C)^2;
        num11=(3^2*omega^2*L^2+(GL)^2)*(-(3^2*omega^2*L^2*C*s+(GL)^2*C*s+GL));
        num12=-(3^2*omega^2*L^2+(GL)^2)*(3^2*omega^2*L^2*omega*C-3*omega*L+(GL)^2*omega*C);
        num21=-num12;
        num22=num11;
        Zo_cal=[num11/den1 num12/den1;num21/den1 num22/den1];
        w = 1:5:1e4;
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
        Gid_cal_frd_rec{i}=Gid_cal_frd;
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
        Gvd_cal_frd_rec{i}=Gvd_cal_frd;
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
        Gig_cal_frd_rec{i}=Gig_cal_frd;
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
        Gvg_cal_frd_rec{i}=Gvg_cal_frd;
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
        Gcv = [kpv+kiv/s tf([0 0],[0 1]); tf([0 0],[0 1]) kpv+kiv/s];
        Gcv_frd = [frd(freqresp(Gcv(1,1),w),w) frd(freqresp(Gcv(1,2),w),w);...
            frd(freqresp(Gcv(2,1),w),w) frd(freqresp(Gcv(2,2),w),w)];
        %% calculation of VSI output impedance with current loop control
        Tv_frd=I_frd+Gvd_cal_frd*((I_frd+Gdel_frd*(Gci_frd-Gdei_frd)*Ki_frd*Gid_cal_frd)\Gdel_frd*Gci_frd*Gcv_frd*Kv_frd);

        Zo_il_cal_frd = Gvd_cal_frd*((I_frd+Gdel_frd*(-Gdei_frd+Gci_frd)*Ki_frd*Gid_cal_frd)\Gdel_frd*(Gdei_frd-Gci_frd)*Ki_frd*Gii_cal_frd)+Zo_cal_frd;

        Zo_vil_cal_frd{i}=Tv_frd\Zo_il_cal_frd;
        Zo_cal_frd_open{i}=Zo_cal_frd;
end

%% calculate the standard deviation of Zovil_cal_frd 
% since Zovil_cal_frd are complex numbers, there are two ways
% to separate them, 1. real and imaginary 2. magnitude and angle
% for the following the first method is used

% for i=1:1:nsamples
%     a_dd(:,i)=real(squeeze(Zo_vil_cal_frd{i}(1,1).ResponseData));
%     b_dd(:,i)=imag(squeeze(Zo_vil_cal_frd{i}(1,1).ResponseData));
%     a_dq(:,i)=real(squeeze(Zo_vil_cal_frd{i}(1,2).ResponseData));
%     b_dq(:,i)=imag(squeeze(Zo_vil_cal_frd{i}(1,2).ResponseData));
%     a_qd(:,i)=real(squeeze(Zo_vil_cal_frd{i}(2,1).ResponseData));
%     b_qd(:,i)=imag(squeeze(Zo_vil_cal_frd{i}(2,1).ResponseData));
%     a_qq(:,i)=real(squeeze(Zo_vil_cal_frd{i}(2,2).ResponseData));
%     b_qq(:,i)=imag(squeeze(Zo_vil_cal_frd{i}(2,2).ResponseData));
% end
% for i=1:1:length(a_dd(:,1))
%     a_dd_mean(i)=mean(a_dd(i,:));
%     a_dd_std(i)=std(a_dd(i,:));
%     a_dd_95_upper(i)=a_dd_mean(i)+1.95996*a_dd_std(i);
%     a_dd_95_lower(i)=a_dd_mean(i)-1.95996*a_dd_std(i);
%     b_dd_mean(i)=mean(b_dd(i,:));
%     b_dd_std(i)=std(b_dd(i,:));
%     b_dd_95_upper(i)=b_dd_mean(i)+1.95996*b_dd_std(i);
%     b_dd_95_lower(i)=b_dd_mean(i)-1.95996*b_dd_std(i);
%     
%     a_dq_mean(i)=mean(a_dq(i,:));
%     a_dq_std(i)=std(a_dq(i,:));
%     a_dq_95_upper(i)=a_dq_mean(i)+1.95996*a_dq_std(i);
%     a_dq_95_lower(i)=a_dq_mean(i)-1.95996*a_dq_std(i);
%     b_dq_mean(i)=mean(b_dq(i,:));
%     b_dq_std(i)=std(b_dq(i,:));
%     b_dq_95_upper(i)=b_dq_mean(i)+1.95996*b_dq_std(i);
%     b_dq_95_lower(i)=b_dq_mean(i)-1.95996*b_dq_std(i);
%     a_qd_mean(i)=mean(a_qd(i,:));
%     a_qd_std(i)=std(a_qd(i,:));
%     a_qd_95_upper(i)=a_qd_mean(i)+1.95996*a_qd_std(i);
%     a_qd_95_lower(i)=a_qd_mean(i)-1.95996*a_qd_std(i);
%     b_qd_mean(i)=mean(b_qd(i,:));
%     b_qd_std(i)=std(b_qd(i,:));
%     b_qd_95_upper(i)=b_qd_mean(i)+1.95996*b_qd_std(i);
%     b_qd_95_lower(i)=b_qd_mean(i)-1.95996*b_qd_std(i);
%     a_qq_mean(i)=mean(a_qq(i,:));
%     a_qq_std(i)=std(a_qq(i,:));
%     a_qq_95_upper(i)=a_qq_mean(i)+1.95996*a_qq_std(i);
%     a_qq_95_lower(i)=a_qq_mean(i)-1.95996*a_qq_std(i);
%     b_qq_mean(i)=mean(b_qq(i,:));
%     b_qq_std(i)=std(b_qq(i,:));
%     b_qq_95_upper(i)=b_qq_mean(i)+1.95996*b_qq_std(i);
%     b_qq_95_lower(i)=b_qq_mean(i)-1.95996*b_qq_std(i);
% end
% Zo_vil_cal_frd_dd_mean=frd(complex(a_dd_mean,b_dd_mean),w);
% Zo_vil_cal_frd_dd_95_upper=frd(complex(a_dd_95_upper,b_dd_95_upper),w);
% Zo_vil_cal_frd_dd_95_lower=frd(complex(a_dd_95_lower,b_dd_95_lower),w);
% Zo_vil_cal_frd_dq_mean=frd(complex(a_dq_mean,b_dq_mean),w);
% Zo_vil_cal_frd_dq_95_upper=frd(complex(a_dq_95_upper,b_dq_95_upper),w);
% Zo_vil_cal_frd_dq_95_lower=frd(complex(a_dq_95_lower,b_dq_95_lower),w);
% Zo_vil_cal_frd_qd_mean=frd(complex(a_qd_mean,b_qd_mean),w);
% Zo_vil_cal_frd_qd_95_upper=frd(complex(a_qd_95_upper,b_qd_95_upper),w);
% Zo_vil_cal_frd_qd_95_lower=frd(complex(a_qd_95_lower,b_qd_95_lower),w);
% Zo_vil_cal_frd_qq_mean=frd(complex(a_qq_mean,b_qq_mean),w);
% Zo_vil_cal_frd_qq_95_upper=frd(complex(a_qq_95_upper,b_qq_95_upper),w);
% Zo_vil_cal_frd_qq_95_lower=frd(complex(a_qq_95_lower,b_qq_95_lower),w);

for i=1:1:nsamples
    a_dd(:,i)=abs(squeeze(Zo_vil_cal_frd{i}(1,1).ResponseData));
    b_dd(:,i)=angle(squeeze(Zo_vil_cal_frd{i}(1,1).ResponseData));
    a_dq(:,i)=abs(squeeze(Zo_vil_cal_frd{i}(1,2).ResponseData));
    b_dq(:,i)=angle(squeeze(Zo_vil_cal_frd{i}(1,2).ResponseData));
    a_qd(:,i)=abs(squeeze(Zo_vil_cal_frd{i}(2,1).ResponseData));
    b_qd(:,i)=angle(squeeze(Zo_vil_cal_frd{i}(2,1).ResponseData));
    a_qq(:,i)=abs(squeeze(Zo_vil_cal_frd{i}(2,2).ResponseData));
    b_qq(:,i)=angle(squeeze(Zo_vil_cal_frd{i}(2,2).ResponseData));
end
for i=1:1:length(a_dd(:,1))
    a_dd_mean(i)=mean(a_dd(i,:));
    a_dd_std(i)=std(a_dd(i,:));
    a_dd_max(i)=max(a_dd(i,:));
    a_dd_min(i)=min(a_dd(i,:));
    a_dd_95_upper(i)=a_dd_mean(i)+1.95996*a_dd_std(i);
    a_dd_95_lower(i)=a_dd_mean(i)-1.95996*a_dd_std(i);
    b_dd_mean(i)=mean(b_dd(i,:));
    b_dd_std(i)=std(b_dd(i,:));
    b_dd_max(i)=max(b_dd(i,:));
    b_dd_min(i)=min(b_dd(i,:));
    b_dd_95_upper(i)=b_dd_mean(i)+1.95996*b_dd_std(i);
    b_dd_95_lower(i)=b_dd_mean(i)-1.95996*b_dd_std(i);
    
    a_dq_mean(i)=mean(a_dq(i,:));
    a_dq_std(i)=std(a_dq(i,:));
    a_dq_95_upper(i)=a_dq_mean(i)+1.95996*a_dq_std(i);
    a_dq_95_lower(i)=a_dq_mean(i)-1.95996*a_dq_std(i);
    b_dq_mean(i)=mean(b_dq(i,:));
    b_dq_std(i)=std(b_dq(i,:));
    b_dq_95_upper(i)=b_dq_mean(i)+1.95996*b_dq_std(i);
    b_dq_95_lower(i)=b_dq_mean(i)-1.95996*b_dq_std(i);
    a_qd_mean(i)=mean(a_qd(i,:));
    a_qd_std(i)=std(a_qd(i,:));
    a_qd_95_upper(i)=a_qd_mean(i)+1.95996*a_qd_std(i);
    a_qd_95_lower(i)=a_qd_mean(i)-1.95996*a_qd_std(i);
    b_qd_mean(i)=mean(b_qd(i,:));
    b_qd_std(i)=std(b_qd(i,:));
    b_qd_95_upper(i)=b_qd_mean(i)+1.95996*b_qd_std(i);
    b_qd_95_lower(i)=b_qd_mean(i)-1.95996*b_qd_std(i);
    a_qq_mean(i)=mean(a_qq(i,:));
    a_qq_std(i)=std(a_qq(i,:));
    a_qq_95_upper(i)=a_qq_mean(i)+1.95996*a_qq_std(i);
    a_qq_95_lower(i)=a_qq_mean(i)-1.95996*a_qq_std(i);
    b_qq_mean(i)=mean(b_qq(i,:));
    b_qq_std(i)=std(b_qq(i,:));
    b_qq_95_upper(i)=b_qq_mean(i)+1.95996*b_qq_std(i);
    b_qq_95_lower(i)=b_qq_mean(i)-1.95996*b_qq_std(i);
end
Zo_vil_cal_frd_mean_mag=a_dd_mean;
Zo_vil_cal_frd_95_upper_mag=a_dd_95_upper;
Zo_vil_cal_frd_95_lower_mag=a_dd_95_lower;
Zo_vil_cal_frd_mean_ang=b_dd_mean;
Zo_vil_cal_frd_95_upper_ang=b_dd_95_upper;
Zo_vil_cal_frd_95_lower_ang=b_dd_95_lower;
figure(1)
semilogx(w/2/pi,20*log10(Zo_vil_cal_frd_mean_mag),'b');
hold on
semilogx(w/2/pi,20*log10(Zo_vil_cal_frd_95_upper_mag),'r');
semilogx(w/2/pi,20*log10(Zo_vil_cal_frd_95_lower_mag),'g');
figure(2)
semilogx(w/2/pi,Zo_vil_cal_frd_mean_ang*180/pi,'b');
hold on
semilogx(w/2/pi,Zo_vil_cal_frd_95_upper_ang*180/pi,'r');
semilogx(w/2/pi,Zo_vil_cal_frd_95_lower_ang*180/pi,'g');
figure(51)
semilogx(w/2/pi,20*log10(Zo_vil_cal_frd_mean_mag),'b');
hold on
semilogx(w/2/pi,20*log10(a_dd_max),'r');
semilogx(w/2/pi,20*log10(a_dd_min),'g');
% semilogx(w/2/pi,20*log10(a_dd_std),'k');
figure(41)
semilogx(w/2/pi,Zo_vil_cal_frd_mean_ang*180/pi,'b');
hold on
semilogx(w/2/pi,b_dd_max*180/pi,'r');
semilogx(w/2/pi,b_dd_min*180/pi,'g');
% semilogx(w/2/pi,b_dd_std*180/pi,'k');

for i=1:1:nsamples
    a_dd(:,i)=abs(squeeze(Zo_cal_frd_open{i}(1,1).ResponseData));
    b_dd(:,i)=angle(squeeze(Zo_cal_frd_open{i}(1,1).ResponseData));
    a_dq(:,i)=abs(squeeze(Zo_cal_frd_open{i}(1,2).ResponseData));
    b_dq(:,i)=angle(squeeze(Zo_cal_frd_open{i}(1,2).ResponseData));
    a_qd(:,i)=abs(squeeze(Zo_cal_frd_open{i}(2,1).ResponseData));
    b_qd(:,i)=angle(squeeze(Zo_cal_frd_open{i}(2,1).ResponseData));
    a_qq(:,i)=abs(squeeze(Zo_cal_frd_open{i}(2,2).ResponseData));
    b_qq(:,i)=angle(squeeze(Zo_cal_frd_open{i}(2,2).ResponseData));
end
for i=1:1:length(a_dd(:,1))
    a_dd_mean(i)=mean(a_dd(i,:));
    a_dd_std(i)=std(a_dd(i,:));
    a_dd_max(i)=max(a_dd(i,:));
    a_dd_min(i)=min(a_dd(i,:));
    a_dd_95_upper(i)=a_dd_mean(i)+1.95996*a_dd_std(i);
    a_dd_95_lower(i)=a_dd_mean(i)-1.95996*a_dd_std(i);
    b_dd_mean(i)=mean(b_dd(i,:));
    b_dd_std(i)=std(b_dd(i,:));
    b_dd_max(i)=max(b_dd(i,:));
    b_dd_min(i)=min(b_dd(i,:));
    b_dd_95_upper(i)=b_dd_mean(i)+1.95996*b_dd_std(i);
    b_dd_95_lower(i)=b_dd_mean(i)-1.95996*b_dd_std(i);
    
    a_dq_mean(i)=mean(a_dq(i,:));
    a_dq_std(i)=std(a_dq(i,:));
    a_dq_95_upper(i)=a_dq_mean(i)+1.95996*a_dq_std(i);
    a_dq_95_lower(i)=a_dq_mean(i)-1.95996*a_dq_std(i);
    b_dq_mean(i)=mean(b_dq(i,:));
    b_dq_std(i)=std(b_dq(i,:));
    b_dq_95_upper(i)=b_dq_mean(i)+1.95996*b_dq_std(i);
    b_dq_95_lower(i)=b_dq_mean(i)-1.95996*b_dq_std(i);
    a_qd_mean(i)=mean(a_qd(i,:));
    a_qd_std(i)=std(a_qd(i,:));
    a_qd_95_upper(i)=a_qd_mean(i)+1.95996*a_qd_std(i);
    a_qd_95_lower(i)=a_qd_mean(i)-1.95996*a_qd_std(i);
    b_qd_mean(i)=mean(b_qd(i,:));
    b_qd_std(i)=std(b_qd(i,:));
    b_qd_95_upper(i)=b_qd_mean(i)+1.95996*b_qd_std(i);
    b_qd_95_lower(i)=b_qd_mean(i)-1.95996*b_qd_std(i);
    a_qq_mean(i)=mean(a_qq(i,:));
    a_qq_std(i)=std(a_qq(i,:));
    a_qq_95_upper(i)=a_qq_mean(i)+1.95996*a_qq_std(i);
    a_qq_95_lower(i)=a_qq_mean(i)-1.95996*a_qq_std(i);
    b_qq_mean(i)=mean(b_qq(i,:));
    b_qq_std(i)=std(b_qq(i,:));
    b_qq_95_upper(i)=b_qq_mean(i)+1.95996*b_qq_std(i);
    b_qq_95_lower(i)=b_qq_mean(i)-1.95996*b_qq_std(i);
end
Zo_cal_frd_open_mean_mag=a_dd_mean;
Zo_cal_frd_open_95_upper_mag=a_dd_95_upper;
Zo_cal_frd_open_95_lower_mag=a_dd_95_lower;
Zo_cal_frd_open_mean_ang=b_dd_mean;
Zo_cal_frd_open_95_upper_ang=b_dd_95_upper;
Zo_cal_frd_open_95_lower_ang=b_dd_95_lower;
figure(3)
semilogx(w/2/pi,20*log10(Zo_cal_frd_open_mean_mag),'b');
hold on
semilogx(w/2/pi,20*log10(Zo_cal_frd_open_95_upper_mag),'r');
semilogx(w/2/pi,20*log10(Zo_cal_frd_open_95_lower_mag),'g');
semilogx(w/2/pi,20*log10(a_dd),'k');
figure(4)
semilogx(w/2/pi,Zo_cal_frd_open_mean_ang*180/pi,'b');
hold on
semilogx(w/2/pi,Zo_cal_frd_open_95_upper_ang*180/pi,'r');
semilogx(w/2/pi,Zo_cal_frd_open_95_lower_ang*180/pi,'g');
semilogx(w/2/pi,b_dd*180/pi,'k');

figure(5)
semilogx(w/2/pi,20*log10(Zo_cal_frd_open_mean_mag),'b');
hold on
semilogx(w/2/pi,20*log10(a_dd_max),'r');
semilogx(w/2/pi,20*log10(a_dd_min),'g');
% semilogx(w/2/pi,20*log10(a_dd_std),'k');
figure(4)
semilogx(w/2/pi,Zo_cal_frd_open_mean_ang*180/pi,'b');
hold on
semilogx(w/2/pi,b_dd_max*180/pi,'r');
semilogx(w/2/pi,b_dd_min*180/pi,'g');
% semilogx(w/2/pi,b_dd_std*180/pi,'k');