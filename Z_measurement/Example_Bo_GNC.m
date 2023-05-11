%% verification of Bar-on example:
clc
close all
clear all

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
Bode_O.PhaseWrapping='on';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
linewidth=1; fontsize=18;
linestyle1='b'; linestyle2='g';
linestyle3='r'; linestyle4='k';
linestyle5='m'; linestyle6='c';

linestyle=[linestyle1; linestyle2; linestyle3; linestyle4; linestyle5; linestyle6];

% define GNC plot frequency range
k_freq_r = .5;%1;%0.6;

load ('Z_vsi_100p_afe_fi_1k_fv70_vsi_fi_4k_fv_250.mat')
Z_VSI_o_400_fi_4k_fv_250_test = Z.ZDQ/3;%Z.ZDQRAW/3;%
load ('Z_vsi_100p_afe_fi_1k_fv70_vsi_fi_4k_fv_200.mat')
Z_VSI_o_400_fi_4k_fv_200_test = Z.ZDQ/3;%Z.ZDQRAW/3;%
load ('Z_vsi_100p_afe_fi_1k_fv70_vsi_fi_4k_fv_150.mat')
Z_VSI_o_400_fi_4k_fv_150_test = Z.ZDQ/3;%Z.ZDQRAW/3;%
load ('Z_vsi_100p_afe_fi_1k_fv70_vsi_fi_4k_fv_110.mat')
Z_VSI_o_400_fi_4k_fv_110_test = Z.ZDQ/3;%Z.ZDQRAW/3;%
load ('Z_vsi_100p_afe_fi_1k_fv70_vsi_fi_4k_fv_100.mat')
Z_VSI_o_400_fi_4k_fv_100_test = Z.ZDQ/3;%Z.ZDQRAW/3;%
load ('Z_vsi_100p_afe_fi_1k_fv70_vsi_fi_4k_fv_80.mat')
Z_VSI_o_400_fi_4k_fv_80_test = Z.ZDQ/3;%Z.ZDQRAW/3;%

load ('Z_afe_100p_vsi_fi_4k_fv200_afe_fi_1k_fv_70.mat')
Z_AFE_400_fi_1000_fv_70_pll_1_test = Z.ZDQRAW/3;

%% GNC plot using measured exprimental results
% to use the exprimental data, first task to do is to reconstruct the test
% data using the dame frequency list. although all the test datas are
% collected using 50 points between 40 Hz and 10 kHz, the frequency value
% of different test may varies for a little bit
freq_test = Z_VSI_o_400_fi_4k_fv_80_test(1,1).Frequency;
Z_VSI_o_400_fi_4k_fv_80_test_1 = Z_VSI_o_400_fi_4k_fv_80_test;
% reconstruct every test data:
Z_dd_resp = Z_VSI_o_400_fi_4k_fv_100_test(1,1).ResponseData;
Z_dd = frd(Z_dd_resp,freq_test);
Z_dq_resp = Z_VSI_o_400_fi_4k_fv_100_test(1,2).ResponseData;
Z_dq = frd(Z_dq_resp,freq_test);
Z_qd_resp = Z_VSI_o_400_fi_4k_fv_100_test(2,1).ResponseData;
Z_qd = frd(Z_qd_resp,freq_test);
Z_qq_resp = Z_VSI_o_400_fi_4k_fv_100_test(2,2).ResponseData;
Z_qq = frd(Z_qq_resp,freq_test);
Z_VSI_o_400_fi_4k_fv_100_test_1 = [Z_dd Z_dq; Z_qd Z_qq];

Z_dd_resp = Z_VSI_o_400_fi_4k_fv_110_test(1,1).ResponseData;
Z_dd = frd(Z_dd_resp,freq_test);
Z_dq_resp = Z_VSI_o_400_fi_4k_fv_110_test(1,2).ResponseData;
Z_dq = frd(Z_dq_resp,freq_test);
Z_qd_resp = Z_VSI_o_400_fi_4k_fv_110_test(2,1).ResponseData;
Z_qd = frd(Z_qd_resp,freq_test);
Z_qq_resp = Z_VSI_o_400_fi_4k_fv_110_test(2,2).ResponseData;
Z_qq = frd(Z_qq_resp,freq_test);
Z_VSI_o_400_fi_4k_fv_110_test_1 = [Z_dd Z_dq; Z_qd Z_qq];

Z_dd_resp = Z_VSI_o_400_fi_4k_fv_150_test(1,1).ResponseData;
Z_dd = frd(Z_dd_resp,freq_test);
Z_dq_resp = Z_VSI_o_400_fi_4k_fv_150_test(1,2).ResponseData;
Z_dq = frd(Z_dq_resp,freq_test);
Z_qd_resp = Z_VSI_o_400_fi_4k_fv_150_test(2,1).ResponseData;
Z_qd = frd(Z_qd_resp,freq_test);
Z_qq_resp = Z_VSI_o_400_fi_4k_fv_150_test(2,2).ResponseData;
Z_qq = frd(Z_qq_resp,freq_test);
Z_VSI_o_400_fi_4k_fv_150_test_1 = [Z_dd Z_dq; Z_qd Z_qq];


Z_dd_resp = Z_VSI_o_400_fi_4k_fv_200_test(1,1).ResponseData;
Z_dd = frd(Z_dd_resp,freq_test);
Z_dq_resp = Z_VSI_o_400_fi_4k_fv_200_test(1,2).ResponseData;
Z_dq = frd(Z_dq_resp,freq_test);
Z_qd_resp = Z_VSI_o_400_fi_4k_fv_200_test(2,1).ResponseData;
Z_qd = frd(Z_qd_resp,freq_test);
Z_qq_resp = Z_VSI_o_400_fi_4k_fv_200_test(2,2).ResponseData;
Z_qq = frd(Z_qq_resp,freq_test);
Z_VSI_o_400_fi_4k_fv_200_test_1 = [Z_dd Z_dq; Z_qd Z_qq];


Z_dd_resp = Z_VSI_o_400_fi_4k_fv_250_test(1,1).ResponseData;
Z_dd = frd(Z_dd_resp,freq_test);
Z_dq_resp = Z_VSI_o_400_fi_4k_fv_250_test(1,2).ResponseData;
Z_dq = frd(Z_dq_resp,freq_test);
Z_qd_resp = Z_VSI_o_400_fi_4k_fv_250_test(2,1).ResponseData;
Z_qd = frd(Z_qd_resp,freq_test);
Z_qq_resp = Z_VSI_o_400_fi_4k_fv_250_test(2,2).ResponseData;
Z_qq = frd(Z_qq_resp,freq_test);
Z_VSI_o_400_fi_4k_fv_250_test_1 = [Z_dd Z_dq; Z_qd Z_qq];


Z_dd_resp = Z_AFE_400_fi_1000_fv_70_pll_1_test(1,1).ResponseData;
Z_dd = frd(Z_dd_resp,freq_test);
Z_dq_resp = Z_AFE_400_fi_1000_fv_70_pll_1_test(1,2).ResponseData;
Z_dq = frd(Z_dq_resp,freq_test);
Z_qd_resp = Z_AFE_400_fi_1000_fv_70_pll_1_test(2,1).ResponseData;
Z_qd = frd(Z_qd_resp,freq_test);
Z_qq_resp = Z_AFE_400_fi_1000_fv_70_pll_1_test(2,2).ResponseData;
Z_qq = frd(Z_qq_resp,freq_test);
Z_AFE_400_fi_1000_fv_70_pll_1_test = [Z_dd Z_dq; Z_qd Z_qq];
Yin_vil_pll_cal_frd = 1/Z_AFE_400_fi_1000_fv_70_pll_1_test;           % afe is the same

freq_test = Z_VSI_o_400_fi_4k_fv_80_test(1,1).Frequency;

k_d = +5;
k_g = 1.1;

% six case VSI
Zo_vsi_vil = Z_VSI_o_400_fi_4k_fv_250_test_1;                           % then the unstable vsi case
Zo_vsi_vil_s = [Zo_vsi_vil(1,1) 0; 0 Zo_vsi_vil(2,2)];
Yin_vil_pll_cal_frd_s = [Yin_vil_pll_cal_frd(1,1) 0; 0 Yin_vil_pll_cal_frd(2,2)];

% RR is the full GNC return ration
RR=Zo_vsi_vil*Yin_vil_pll_cal_frd;
% rr is the return ratio defined in Rolando's paper
rr=Zo_vsi_vil_s*Yin_vil_pll_cal_frd_s;
RRdd=rr(1,1);
RRdq=rr(1,2);
RRqd=rr(2,1);
RRqq=rr(2,2);
% calculate ACindex defined by Rolando
ACindex=(RRdd-RRqq)*(RRdd-RRqq)/(4*RRqd*RRdq);

% Eigenvalues and Sorting
RR_res=[RR(1,1).ResponseData RR(1,2).ResponseData;RR(2,1).ResponseData RR(2,2).ResponseData];
I = [tf([0 1],[0 1]), 0; 0 tf([0 1],[0 1])];
Resp_I = freqresp(I,freq_test);
M_res=RR_res./(Resp_I+RR_res);
% [RR_res]=sortloci(RR_res);
for k=1:length(freq_test)*k_freq_r
%         Lddeigenvalues(:,k)=eig(Lddresp(:,:,k));
%         Lqqeigenvalues(:,k)=eig(Lqqresp(:,:,k));
    Leigenvalues(:,k)=eig(RR_res(:,:,k));
    Singular_L(:,k) = svd(RR_res(:,:,k));
end
[Leigenvalues]=sortloci(Leigenvalues);
% plot singular value and abs of eigenvalue
figure(1)
semilogx(freq_test(1:length(freq_test)*k_freq_r),Singular_L(1,:),'b*');
hold on
semilogx(freq_test(1:length(freq_test)*k_freq_r),Singular_L(2,:),'b*');
semilogx(freq_test(1:length(freq_test)*k_freq_r),abs(Leigenvalues(1,:)),'k');
semilogx(freq_test(1:length(freq_test)*k_freq_r),abs(Leigenvalues(2,:)),'k');
grid on


% plot the full GNC Characteristic Loci
figHandle=figure(2);
set(figHandle, 'Position', [50, 10, 1000, 800]);
plot(Leigenvalues(1,:),linestyle1,'LineWidth',linewidth)
hold on

% plot(conj(Leigenvalues(1,1:length(freq_test)*k_freq_r)),linestyle1,'LineWidth',linewidth)
plot(Leigenvalues(2,:),linestyle2,'LineWidth',linewidth)
% plot(conj(Leigenvalues(2,1:length(freq_test)*k_freq_r)),linestyle2,'LineWidth',linewidth)

legend({'{\it\lambda}_{1}','{\it\lambda}_{1}_{mirror}','{\it\lambda}_{2}','{\it\lambda}_{2}_{mirror}'},'Fontsize',fontsize,'FontWeight','bold')
plot(-1,0,'r+','LineWidth',linewidth)
% hold off
title('Characteristic Loci of Ldd','Fontsize',fontsize,'FontWeight','bold')
grid
set(gca,'FontSize',fontsize);
ezplot('x^2+y^2=1')

Mag=abs(Leigenvalues(2,:));
Phase=angle(Leigenvalues(2,:))*180/pi;
figure(3)
subplot(2,1,1)
semilogx(freq_test(1:length(freq_test)*k_freq_r),20*log10(Mag),'b*');
grid on
subplot(2,1,2)
semilogx(freq_test(1:length(freq_test)*k_freq_r),Phase,'r*');
hold off
grid on


angle_eigs_L=angle(Leigenvalues);
Mag_eigs_L=abs(Leigenvalues);

% for phase margin, find at which frequency Mag of L' eig reach unity
% circle

index_1=find(Mag_eigs_L(1,:)>0.99&Mag_eigs_L(1,:)<1.01);
index_2=find(Mag_eigs_L(2,:)>0.99&Mag_eigs_L(2,:)<1.01);
PM_1=(-angle_eigs_L(1,index_1)+pi)*180/pi;
PM_2=(-angle_eigs_L(2,index_2)+pi)*180/pi;

if isempty(PM_1)
    [PM Index]= min(PM_2);
elseif isempty(PM_2)
    [PM Index] = min(PM_1);
elseif min(PM_1)>=min(PM_2)
    [PM Index] = min(PM_2);
else
    [PM Index] = min(PM_1);
end
Phase_Margin=PM
% [eig_L]=sortloci(eig_L);
figure(4)
plot(Leigenvalues(1,1:length(Leigenvalues)),'b')
hold on
plot(Leigenvalues(2,1:length(Leigenvalues)),'r')
grid on
ezplot('x^2+y^2=1')

U=[exp(1i*PM*pi/180) 0;0 exp(1i*PM*pi/180)];
for k=1:length(freq_test)*k_freq_r
   eig_L_1(:,k) = eig(RR_res(:,:,k)*(U));
end
[eig_L_1]=sortloci(eig_L_1);
% figure(41)
plot(eig_L_1(1,1:length(eig_L_1)),'b*')
% hold on
plot(eig_L_1(2,1:length(eig_L_1)),'r*')
% grid on
% U_1=[exp(-1*(PM-k_d)*pi/180*s/freq_r(index)) 0;0 exp(-1*(PM-k_d)*pi/180*s/freq_r(index))];


% % step response:
% L_1=L*U_1;
% M_1=L_1/(I+L_1);
% figure(55)
% step(M)
% hold on
% step(M_1)
% grid on

% for gain margin, find at which frequency phase of L' eig reach -pi

index_1=find(angle_eigs_L(1,:)>-pi-0.2&angle_eigs_L(1,:)<-pi+0.2);
index_2=find(angle_eigs_L(2,:)>-pi-0.2&angle_eigs_L(2,:)<-pi+0.2);
GM_1=(1./Mag_eigs_L(1,index_1));
GM_2=(1./Mag_eigs_L(2,index_2));

if isempty(GM_1)
    GM = min(GM_2);
elseif isempty(GM_2)
    GM = min(GM_1);
elseif min(GM_1)>=min(GM_2)
    GM = min(GM_2);
else
    GM = min(GM_1);
end
Gain_Margin=20*log10(GM)

figure(5)
plot(Leigenvalues(1,1:length(Leigenvalues)),'b')
hold on
plot(Leigenvalues(2,1:length(Leigenvalues)),'r')
grid on
ezplot('x^2+y^2=1')
R=[GM 0;0 GM];
for k=1:length(freq_test)*k_freq_r
   eig_L_1(:,k) = eig(RR_res(:,:,k)*(R));
end
[eig_L_1]=sortloci(eig_L_1);
% figure(41)
plot(eig_L_1(1,1:length(eig_L_1)),'b*')
hold on
plot(eig_L_1(2,1:length(eig_L_1)),'r*')
grid on



% forth case VSI
Zo_vsi_vil = Z_VSI_o_400_fi_4k_fv_150_test_1;                           % then the unstable vsi case
Zo_vsi_vil_s = [Zo_vsi_vil(1,1) 0; 0 Zo_vsi_vil(2,2)];
Yin_vil_pll_cal_frd_s = [Yin_vil_pll_cal_frd(1,1) 0; 0 Yin_vil_pll_cal_frd(2,2)];

% RR is the full GNC return ration
RR=Zo_vsi_vil*Yin_vil_pll_cal_frd;
% rr is the return ratio defined in Rolando's paper
rr=Zo_vsi_vil_s*Yin_vil_pll_cal_frd_s;
RRdd=rr(1,1);
RRdq=rr(1,2);
RRqd=rr(2,1);
RRqq=rr(2,2);
% calculate ACindex defined by Rolando
ACindex=(RRdd-RRqq)*(RRdd-RRqq)/(4*RRqd*RRdq);

% Eigenvalues and Sorting
RR_res=[RR(1,1).ResponseData RR(1,2).ResponseData;RR(2,1).ResponseData RR(2,2).ResponseData];
I = [tf([0 1],[0 1]), 0; 0 tf([0 1],[0 1])];
Resp_I = freqresp(I,freq_test);
M_res=RR_res./(Resp_I+RR_res);
% [RR_res]=sortloci(RR_res);
for k=1:length(freq_test)*k_freq_r
%         Lddeigenvalues(:,k)=eig(Lddresp(:,:,k));
%         Lqqeigenvalues(:,k)=eig(Lqqresp(:,:,k));
    Leigenvalues(:,k)=eig(RR_res(:,:,k));
    Singular_L(:,k) = svd(RR_res(:,:,k));
end
[Leigenvalues]=sortloci(Leigenvalues);
% plot singular value and abs of eigenvalue
figure(11)
semilogx(freq_test(1:length(freq_test)*k_freq_r),Singular_L(1,:),'b*');
hold on
semilogx(freq_test(1:length(freq_test)*k_freq_r),Singular_L(2,:),'b*');
semilogx(freq_test(1:length(freq_test)*k_freq_r),abs(Leigenvalues(1,:)),'k');
semilogx(freq_test(1:length(freq_test)*k_freq_r),abs(Leigenvalues(2,:)),'k');
grid on


% plot the full GNC Characteristic Loci
figHandle=figure(21);
set(figHandle, 'Position', [50, 10, 1000, 800]);
plot(Leigenvalues(1,:),linestyle1,'LineWidth',linewidth)
hold on

% plot(conj(Leigenvalues(1,1:length(freq_test)*k_freq_r)),linestyle1,'LineWidth',linewidth)
plot(Leigenvalues(2,:),linestyle2,'LineWidth',linewidth)
% plot(conj(Leigenvalues(2,1:length(freq_test)*k_freq_r)),linestyle2,'LineWidth',linewidth)

legend({'{\it\lambda}_{1}','{\it\lambda}_{1}_{mirror}','{\it\lambda}_{2}','{\it\lambda}_{2}_{mirror}'},'Fontsize',fontsize,'FontWeight','bold')
plot(-1,0,'r+','LineWidth',linewidth)
% hold off
title('Characteristic Loci of Ldd','Fontsize',fontsize,'FontWeight','bold')
grid
set(gca,'FontSize',fontsize);
ezplot('x^2+y^2=1')

Mag=abs(Leigenvalues(2,:));
Phase=angle(Leigenvalues(2,:))*180/pi;
figure(31)
subplot(2,1,1)
semilogx(freq_test(1:length(freq_test)*k_freq_r),20*log10(Mag),'b*');
grid on
subplot(2,1,2)
semilogx(freq_test(1:length(freq_test)*k_freq_r),Phase,'r*');
hold off
grid on


angle_eigs_L=angle(Leigenvalues);
Mag_eigs_L=abs(Leigenvalues);

% for phase margin, find at which frequency Mag of L' eig reach unity
% circle

index_1=find(Mag_eigs_L(1,:)>0.9&Mag_eigs_L(1,:)<1.1);
index_2=find(Mag_eigs_L(2,:)>0.9&Mag_eigs_L(2,:)<1.1);
PM_1=(-angle_eigs_L(1,index_1)+pi)*180/pi;
PM_2=(-angle_eigs_L(2,index_2)+pi)*180/pi;

if isempty(PM_1)
    [PM Index]= min(PM_2);
elseif isempty(PM_2)
    [PM Index] = min(PM_1);
elseif min(PM_1)>=min(PM_2)
    [PM Index] = min(PM_2);
else
    [PM Index] = min(PM_1);
end
Phase_Margin=PM
% [eig_L]=sortloci(eig_L);
figure(41)
plot(Leigenvalues(1,1:length(Leigenvalues)),'b')
hold on
plot(Leigenvalues(2,1:length(Leigenvalues)),'r')
grid on
ezplot('x^2+y^2=1')

U=[exp(1i*PM*pi/180) 0;0 exp(1i*PM*pi/180)];
for k=1:length(freq_test)*k_freq_r
   eig_L_1(:,k) = eig(RR_res(:,:,k)*(U));
end
[eig_L_1]=sortloci(eig_L_1);
% figure(41)
plot(eig_L_1(1,1:length(eig_L_1)),'b*')
% hold on
plot(eig_L_1(2,1:length(eig_L_1)),'r*')
% grid on
% U_1=[exp(-1*(PM-k_d)*pi/180*s/freq_r(index)) 0;0 exp(-1*(PM-k_d)*pi/180*s/freq_r(index))];


% % step response:
% L_1=L*U_1;
% M_1=L_1/(I+L_1);
% figure(55)
% step(M)
% hold on
% step(M_1)
% grid on

% for gain margin, find at which frequency phase of L' eig reach -pi

index_1=find(angle_eigs_L(1,:)>-pi-0.2&angle_eigs_L(1,:)<-pi+0.2);
index_2=find(angle_eigs_L(2,:)>-pi-0.2&angle_eigs_L(2,:)<-pi+0.2);
GM_1=(1./Mag_eigs_L(1,index_1));
GM_2=(1./Mag_eigs_L(2,index_2));

if isempty(GM_1)
    GM = min(GM_2);
elseif isempty(GM_2)
    GM = min(GM_1);
elseif min(GM_1)>=min(GM_2)
    GM = min(GM_2);
else
    GM = min(GM_1);
end
Gain_Margin=20*log10(GM)

figure(51)
plot(Leigenvalues(1,1:length(Leigenvalues)),'b')
hold on
plot(Leigenvalues(2,1:length(Leigenvalues)),'r')
grid on
ezplot('x^2+y^2=1')
R=[GM 0;0 GM];
for k=1:length(freq_test)*k_freq_r
   eig_L_1(:,k) = eig(RR_res(:,:,k)*(R));
end
[eig_L_1]=sortloci(eig_L_1);
% figure(41)
plot(eig_L_1(1,1:length(eig_L_1)),'b*')
hold on
plot(eig_L_1(2,1:length(eig_L_1)),'r*')
grid on
% 
% 
% 
% first case VSI
Zo_vsi_vil = Z_VSI_o_400_fi_4k_fv_80_test_1;                           % then the unstable vsi case
Zo_vsi_vil_s = [Zo_vsi_vil(1,1) 0; 0 Zo_vsi_vil(2,2)];
Yin_vil_pll_cal_frd_s = [Yin_vil_pll_cal_frd(1,1) 0; 0 Yin_vil_pll_cal_frd(2,2)];

% RR is the full GNC return ration
RR=Zo_vsi_vil*Yin_vil_pll_cal_frd;
% rr is the return ratio defined in Rolando's paper
rr=Zo_vsi_vil_s*Yin_vil_pll_cal_frd_s;
RRdd=rr(1,1);
RRdq=rr(1,2);
RRqd=rr(2,1);
RRqq=rr(2,2);
% calculate ACindex defined by Rolando
ACindex=(RRdd-RRqq)*(RRdd-RRqq)/(4*RRqd*RRdq);

% Eigenvalues and Sorting
RR_res=[RR(1,1).ResponseData RR(1,2).ResponseData;RR(2,1).ResponseData RR(2,2).ResponseData];
I = [tf([0 1],[0 1]), 0; 0 tf([0 1],[0 1])];
Resp_I = freqresp(I,freq_test);
M_res=RR_res./(Resp_I+RR_res);
% [RR_res]=sortloci(RR_res);
for k=1:length(freq_test)*k_freq_r
%         Lddeigenvalues(:,k)=eig(Lddresp(:,:,k));
%         Lqqeigenvalues(:,k)=eig(Lqqresp(:,:,k));
    Leigenvalues(:,k)=eig(RR_res(:,:,k));
    Singular_L(:,k) = svd(RR_res(:,:,k));
end
[Leigenvalues]=sortloci(Leigenvalues);
% plot singular value and abs of eigenvalue
figure(12)
semilogx(freq_test(1:length(freq_test)*k_freq_r),Singular_L(1,:),'b*');
hold on
semilogx(freq_test(1:length(freq_test)*k_freq_r),Singular_L(2,:),'b*');
semilogx(freq_test(1:length(freq_test)*k_freq_r),abs(Leigenvalues(1,:)),'k');
semilogx(freq_test(1:length(freq_test)*k_freq_r),abs(Leigenvalues(2,:)),'k');
grid on


% plot the full GNC Characteristic Loci
figHandle=figure(22);
set(figHandle, 'Position', [50, 10, 1000, 800]);
plot(Leigenvalues(1,:),linestyle1,'LineWidth',linewidth)
hold on

% plot(conj(Leigenvalues(1,1:length(freq_test)*k_freq_r)),linestyle1,'LineWidth',linewidth)
plot(Leigenvalues(2,:),linestyle2,'LineWidth',linewidth)
% plot(conj(Leigenvalues(2,1:length(freq_test)*k_freq_r)),linestyle2,'LineWidth',linewidth)

legend({'{\it\lambda}_{1}','{\it\lambda}_{1}_{mirror}','{\it\lambda}_{2}','{\it\lambda}_{2}_{mirror}'},'Fontsize',fontsize,'FontWeight','bold')
plot(-1,0,'r+','LineWidth',linewidth)
% hold off
title('Characteristic Loci of Ldd','Fontsize',fontsize,'FontWeight','bold')
grid
set(gca,'FontSize',fontsize);
ezplot('x^2+y^2=1')

Mag=abs(Leigenvalues(2,:));
Phase=angle(Leigenvalues(2,:))*180/pi;
figure(32)
subplot(2,1,1)
semilogx(freq_test(1:length(freq_test)*k_freq_r),20*log10(Mag),'b*');
grid on
subplot(2,1,2)
semilogx(freq_test(1:length(freq_test)*k_freq_r),Phase,'r*');
hold off
grid on


angle_eigs_L=angle(Leigenvalues);
Mag_eigs_L=abs(Leigenvalues);

% for phase margin, find at which frequency Mag of L' eig reach unity
% circle

index_1=find(Mag_eigs_L(1,:)>0.99&Mag_eigs_L(1,:)<1.01);
index_2=find(Mag_eigs_L(2,:)>0.99&Mag_eigs_L(2,:)<1.01);
PM_1=(-angle_eigs_L(1,index_1)+pi)*180/pi;
PM_2=(-angle_eigs_L(2,index_2)+pi)*180/pi;

if isempty(PM_1)
    [PM Index]= min(PM_2);
elseif isempty(PM_2)
    [PM Index] = min(PM_1);
elseif min(PM_1)>=min(PM_2)
    [PM Index] = min(PM_2);
else
    [PM Index] = min(PM_1);
end
Phase_Margin=PM
% [eig_L]=sortloci(eig_L);
figure(42)
plot(Leigenvalues(1,1:length(Leigenvalues)),'b')
hold on
plot(Leigenvalues(2,1:length(Leigenvalues)),'r')
grid on
ezplot('x^2+y^2=1')

U=[exp(1i*PM*pi/180) 0;0 exp(1i*PM*pi/180)];
for k=1:length(freq_test)*k_freq_r
   eig_L_1(:,k) = eig(RR_res(:,:,k)*(U));
end
[eig_L_1]=sortloci(eig_L_1);
% figure(41)
plot(eig_L_1(1,1:length(eig_L_1)),'b*')
% hold on
plot(eig_L_1(2,1:length(eig_L_1)),'r*')
% grid on
% U_1=[exp(-1*(PM-k_d)*pi/180*s/freq_r(index)) 0;0 exp(-1*(PM-k_d)*pi/180*s/freq_r(index))];


% % step response:
% L_1=L*U_1;
% M_1=L_1/(I+L_1);
% figure(55)
% step(M)
% hold on
% step(M_1)
% grid on

% for gain margin, find at which frequency phase of L' eig reach -pi

index_1=find(angle_eigs_L(1,:)>-pi-0.2&angle_eigs_L(1,:)<-pi+0.2);
index_2=find(angle_eigs_L(2,:)>-pi-0.2&angle_eigs_L(2,:)<-pi+0.2);
GM_1=(1./Mag_eigs_L(1,index_1));
GM_2=(1./Mag_eigs_L(2,index_2));

if isempty(GM_1)
    GM = min(GM_2);
elseif isempty(GM_2)
    GM = min(GM_1);
elseif min(GM_1)>=min(GM_2)
    GM = min(GM_2);
else
    GM = min(GM_1);
end
Gain_Margin=20*log10(GM)

figure(52)
plot(Leigenvalues(1,1:length(Leigenvalues)),'b')
hold on
plot(Leigenvalues(2,1:length(Leigenvalues)),'r')
grid on
ezplot('x^2+y^2=1')
R=[GM 0;0 GM];
for k=1:length(freq_test)*k_freq_r
   eig_L_1(:,k) = eig(RR_res(:,:,k)*(R));
end
[eig_L_1]=sortloci(eig_L_1);
% figure(41)
plot(eig_L_1(1,1:length(eig_L_1)),'b*')
hold on
plot(eig_L_1(2,1:length(eig_L_1)),'r*')
grid on