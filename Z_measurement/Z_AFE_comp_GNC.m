clc
clear all
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
Bode_O.PhaseWrapping='on';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
linewidth=3; fontsize=18;
linestyle1='b'; linestyle2='g';
linestyle3='r'; linestyle4='k';
linestyle5='m'; linestyle6='c';

linestyle=[linestyle1; linestyle2; linestyle3; linestyle4; linestyle5; linestyle6];
% define GNC plot frequency range
k_freq_r = 0.6;% 0.6;


load ('Z_afe_100p_vsi_fi_4k_fv200_afe_fi_1k_fv_50.mat')
Z_AFE_400_fi_1000_fv_50_pll_1_test = Z.ZDQ/3;%Z.ZDQRAW/3;%

load ('Z_afe_100p_vsi_fi_4k_fv200_afe_fi_1k_fv_100.mat')
Z_AFE_400_fi_1000_fv_100_pll_1_test = Z.ZDQ/3;%Z.ZDQRAW/3;%

load ('Z_afe_100p_vsi_fi_4k_fv200_afe_fi_1k_fv_150.mat')
Z_AFE_400_fi_1000_fv_150_pll_1_test = Z.ZDQ/3;%Z.ZDQRAW/3;%

load ('Z_afe_100p_vsi_fi_4k_fv200_afe_fi_1k_fv_220.mat')
Z_AFE_400_fi_1000_fv_220_pll_1_test = Z.ZDQ/3;%Z.ZDQRAW/3;%

load ('Z_afe_100p_vsi_fi_4k_fv200_afe_fi_1k_fv_260.mat')
Z_AFE_400_fi_1000_fv_260_pll_1_test = Z.ZDQ/3;%Z.ZDQRAW/3;%

load ('Z_afe_100p_vsi_fi_4k_fv200_afe_fi_1k_fv_300.mat')
Z_AFE_400_fi_1000_fv_300_pll_1_test = Z.ZDQ/3;%Z.ZDQRAW/3;%

load ('Z_vsi_100p_afe_fi_1k_fv70_vsi_fi_4k_fv_200.mat')
Z_VSI_o_400_fi_4k_fv_200_test = Z.ZDQ/3;%Z.ZDQRAW/3;%


figure(7)

bode(Z_AFE_400_fi_1000_fv_50_pll_1_test,Bode_O)
hold on
bode(Z_AFE_400_fi_1000_fv_100_pll_1_test,Bode_O)
% bode(Z_AFE_400_fi_1000_fv_150_pll_1_test,Bode_O)
% bode(Z_AFE_400_fi_1000_fv_220_pll_1_test,Bode_O)
bode(Z_AFE_400_fi_1000_fv_260_pll_1_test,Bode_O)
% bode(Z_AFE_400_fi_1000_fv_300_pll_1_test,Bode_O)
bode(Z_VSI_o_400_fi_4k_fv_200_test,Bode_O)
Bode_Darklines(1);

figure(77)

bode(Z_AFE_400_fi_1000_fv_50_pll_1_test(1,1),Bode_O)
hold on
bode(Z_AFE_400_fi_1000_fv_100_pll_1_test(1,1),Bode_O)
% bode(Z_AFE_400_fi_1000_fv_150_pll_1_test(1,1),Bode_O)
bode(Z_AFE_400_fi_1000_fv_220_pll_1_test(1,1),Bode_O)
% bode(Z_AFE_400_fi_1000_fv_260_pll_1_test(1,1),Bode_O)
% bode(Z_AFE_400_fi_1000_fv_300_pll_1_test(1,1),Bode_O)
bode(Z_VSI_o_400_fi_4k_fv_200_test(1,1),Bode_O)
Bode_Darklines(3);

figure(71)

bode(Z_AFE_400_fi_1000_fv_50_pll_1_test(1,1),Bode_O)
hold on
% bode(Z_AFE_400_fi_1000_fv_50_pll_1_test(2,2),Bode_O)
bode(Z_VSI_o_400_fi_4k_fv_200_test(1,1),Bode_O)
% bode(Z_VSI_o_400_fi_4k_fv_200_test(2,2),Bode_O)
Bode_Darklines(3);

figure(711)

bode(Z_AFE_400_fi_1000_fv_50_pll_1_test(2,2),Bode_O)
hold on
bode(Z_VSI_o_400_fi_4k_fv_200_test(2,2),Bode_O)
Bode_Darklines(3);

fighandle=figure(7111);
set(fighandle,'position',[10, 10, 1000, 800])
bode(Z_AFE_400_fi_1000_fv_50_pll_1_test,Bode_O)
hold on
bode(Z_VSI_o_400_fi_4k_fv_200_test,Bode_O)
Bode_Darklines(2);
legend ('Z\_AFE','Z\_VSI')


figure(72)

bode(Z_AFE_400_fi_1000_fv_100_pll_1_test(1,1),Bode_O)
hold on
bode(Z_VSI_o_400_fi_4k_fv_200_test(1,1),Bode_O)
Bode_Darklines(3);

figure(722)

bode(Z_AFE_400_fi_1000_fv_100_pll_1_test(2,2),Bode_O)
hold on
bode(Z_VSI_o_400_fi_4k_fv_200_test(2,2),Bode_O)
Bode_Darklines(3);

fighandle=figure(7222);
set(fighandle,'position',[10, 10, 1000, 800])
bode(Z_AFE_400_fi_1000_fv_100_pll_1_test,Bode_O)
hold on
bode(Z_VSI_o_400_fi_4k_fv_200_test,Bode_O)
Bode_Darklines(2);
legend ('Z\_AFE','Z\_VSI')

figure(73)
bode(Z_AFE_400_fi_1000_fv_220_pll_1_test(1,1),Bode_O)
hold on
bode(Z_VSI_o_400_fi_4k_fv_200_test(1,1),Bode_O)
Bode_Darklines(3);

figure(733)

bode(Z_AFE_400_fi_1000_fv_220_pll_1_test(2,2),Bode_O)
hold on
bode(Z_VSI_o_400_fi_4k_fv_200_test(2,2),Bode_O)
Bode_Darklines(3);

fighandle=figure(7333);
set(fighandle,'position',[10, 10, 1000, 800])
bode(Z_AFE_400_fi_1000_fv_220_pll_1_test,Bode_O)
hold on
bode(Z_VSI_o_400_fi_4k_fv_200_test,Bode_O)
Bode_Darklines(2);
legend ('Z\_AFE','Z\_VSI')
%% GNC plot using measured exprimental results
% to use the exprimental data, first task to do is to reconstruct the test
% data using the dame frequency list. although all the test datas are
% collected using 50 points between 40 Hz and 10 kHz, the frequency value
% of different test may varies for a little bit

freq_test = Z_VSI_o_400_fi_4k_fv_200_test(1,1).Frequency;
Z_VSI_o_400_fi_4k_fv_200_test_1 = Z_VSI_o_400_fi_4k_fv_200_test;
% reconstruct every test data:
Z_dd_resp = Z_AFE_400_fi_1000_fv_50_pll_1_test(1,1).ResponseData;
Z_dd = frd(Z_dd_resp,freq_test);
Z_dq_resp = Z_AFE_400_fi_1000_fv_50_pll_1_test(1,2).ResponseData;
Z_dq = frd(Z_dq_resp,freq_test);
Z_qd_resp = Z_AFE_400_fi_1000_fv_50_pll_1_test(2,1).ResponseData;
Z_qd = frd(Z_qd_resp,freq_test);
Z_qq_resp = Z_AFE_400_fi_1000_fv_50_pll_1_test(2,2).ResponseData;
Z_qq = frd(Z_qq_resp,freq_test);
Z_AFE_400_fi_1000_fv_50_pll_1_test_1 = [Z_dd Z_dq; Z_qd Z_qq];

Z_dd_resp = Z_AFE_400_fi_1000_fv_100_pll_1_test(1,1).ResponseData;
Z_dd = frd(Z_dd_resp,freq_test);
Z_dq_resp = Z_AFE_400_fi_1000_fv_100_pll_1_test(1,2).ResponseData;
Z_dq = frd(Z_dq_resp,freq_test);
Z_qd_resp = Z_AFE_400_fi_1000_fv_100_pll_1_test(2,1).ResponseData;
Z_qd = frd(Z_qd_resp,freq_test);
Z_qq_resp = Z_AFE_400_fi_1000_fv_100_pll_1_test(2,2).ResponseData;
Z_qq = frd(Z_qq_resp,freq_test);
Z_AFE_400_fi_1000_fv_100_pll_1_test_1 = [Z_dd Z_dq; Z_qd Z_qq];


Z_dd_resp = Z_AFE_400_fi_1000_fv_150_pll_1_test(1,1).ResponseData;
Z_dd = frd(Z_dd_resp,freq_test);
Z_dq_resp = Z_AFE_400_fi_1000_fv_150_pll_1_test(1,2).ResponseData;
Z_dq = frd(Z_dq_resp,freq_test);
Z_qd_resp = Z_AFE_400_fi_1000_fv_150_pll_1_test(2,1).ResponseData;
Z_qd = frd(Z_qd_resp,freq_test);
Z_qq_resp = Z_AFE_400_fi_1000_fv_150_pll_1_test(2,2).ResponseData;
Z_qq = frd(Z_qq_resp,freq_test);
Z_AFE_400_fi_1000_fv_150_pll_1_test_1 = [Z_dd Z_dq; Z_qd Z_qq];

Z_dd_resp = Z_AFE_400_fi_1000_fv_220_pll_1_test(1,1).ResponseData;
Z_dd = frd(Z_dd_resp,freq_test);
Z_dq_resp = Z_AFE_400_fi_1000_fv_220_pll_1_test(1,2).ResponseData;
Z_dq = frd(Z_dq_resp,freq_test);
Z_qd_resp = Z_AFE_400_fi_1000_fv_220_pll_1_test(2,1).ResponseData;
Z_qd = frd(Z_qd_resp,freq_test);
Z_qq_resp = Z_AFE_400_fi_1000_fv_220_pll_1_test(2,2).ResponseData;
Z_qq = frd(Z_qq_resp,freq_test);
Z_AFE_400_fi_1000_fv_220_pll_1_test_1 = [Z_dd Z_dq; Z_qd Z_qq];

Z_dd_resp = Z_AFE_400_fi_1000_fv_260_pll_1_test(1,1).ResponseData;
Z_dd = frd(Z_dd_resp,freq_test);
Z_dq_resp = Z_AFE_400_fi_1000_fv_260_pll_1_test(1,2).ResponseData;
Z_dq = frd(Z_dq_resp,freq_test);
Z_qd_resp = Z_AFE_400_fi_1000_fv_260_pll_1_test(2,1).ResponseData;
Z_qd = frd(Z_qd_resp,freq_test);
Z_qq_resp = Z_AFE_400_fi_1000_fv_260_pll_1_test(2,2).ResponseData;
Z_qq = frd(Z_qq_resp,freq_test);
Z_AFE_400_fi_1000_fv_260_pll_1_test_1 = [Z_dd Z_dq; Z_qd Z_qq];

Z_dd_resp = Z_AFE_400_fi_1000_fv_300_pll_1_test(1,1).ResponseData;
Z_dd = frd(Z_dd_resp,freq_test);
Z_dq_resp = Z_AFE_400_fi_1000_fv_300_pll_1_test(1,2).ResponseData;
Z_dq = frd(Z_dq_resp,freq_test);
Z_qd_resp = Z_AFE_400_fi_1000_fv_300_pll_1_test(2,1).ResponseData;
Z_qd = frd(Z_qd_resp,freq_test);
Z_qq_resp = Z_AFE_400_fi_1000_fv_300_pll_1_test(2,2).ResponseData;
Z_qq = frd(Z_qq_resp,freq_test);
Z_AFE_400_fi_1000_fv_300_pll_1_test_1 = [Z_dd Z_dq; Z_qd Z_qq];


% % first case VSI
Zo_vsi_vil = Z_VSI_o_400_fi_4k_fv_200_test_1;                           % then the unstable vsi case
Zo_vsi_vil_s = [Zo_vsi_vil(1,1) 0; 0 Zo_vsi_vil(2,2)];

% second case VSI
Yin_vil_pll_cal_frd = 1/Z_AFE_400_fi_1000_fv_220_pll_1_test_1;           % afe is the same
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
% disable ACindex plot
% figure(3)
% bode(ACindex);
% hold on
% bode(RRdq)
% title('ACindex and Ldd')
% grid on;
% Bode_Darklines(3);

% Eigenvalues and Sorting
RR_res=[RR(1,1).ResponseData RR(1,2).ResponseData;RR(2,1).ResponseData RR(2,2).ResponseData];
rr_res=[rr(1,1).ResponseData rr(1,2).ResponseData;rr(2,1).ResponseData rr(2,2).ResponseData];
for k=1:length(freq_test)*k_freq_r
%         Lddeigenvalues(:,k)=eig(Lddresp(:,:,k));
%         Lqqeigenvalues(:,k)=eig(Lqqresp(:,:,k));
    Leigenvalues(:,k)=eig(RR_res(:,:,k));
    leigenvalues(:,k)=eig(rr_res(:,:,k));
    Singular_L(:,k)=svd(RR_res(:,:,k));
    Singular_l(:,k)=svd(rr_res(:,:,k));
    [U H]=poldec(RR_res(:,:,k));
    eig_U(:,k)=eig(U);
    eig_H(:,k)=eig(H);
end
[Leigenvalues]=sortloci(Leigenvalues);
% a=[flipud(conj(Leigenvalues(1,1:length(freq_test)*k_freq_r))), Leigenvalues(1,1:length(freq_test)*k_freq_r)];
% plot the full GNC Characteristic Loci
figHandle=figure(92);
set(figHandle, 'Position', [50, 10, 1000, 800]);
plot(Leigenvalues(1,1:length(freq_test)*k_freq_r),linestyle1,'LineWidth',linewidth)
hold on
plot(conj(Leigenvalues(1,1:length(freq_test)*k_freq_r)),linestyle1,'LineWidth',linewidth)
plot(Leigenvalues(2,1:length(freq_test)*k_freq_r),linestyle2,'LineWidth',linewidth)
plot(conj(Leigenvalues(2,1:length(freq_test)*k_freq_r)),linestyle2,'LineWidth',linewidth)

legend({'{\it\lambda}_{1}','{\it\lambda}_{1}_{mirror}','{\it\lambda}_{2}','{\it\lambda}_{2}_{mirror}'},'Fontsize',fontsize,'FontWeight','bold')
plot(-1,0,'r+','LineWidth',linewidth)
% hold off
title('Characteristic Loci of Ldd','Fontsize',fontsize,'FontWeight','bold')
grid
set(gca,'FontSize',fontsize);
ezplot('x^2+y^2=1')
%% find the angle of eigenvalues
angle_lamda1_case1=360+angle(Leigenvalues(1,:))'*180/pi;
%% bode plot of eigenvalues
% lamda1_case1 = frd(Leigenvalues(1,1:length(freq_test)*k_freq_r),freq_test(1:length(freq_test)*k_freq_r));
% lamda2_case1 = frd(Leigenvalues(2,1:length(freq_test)*k_freq_r),freq_test(1:length(freq_test)*k_freq_r));
gershgorinplot(RR,Leigenvalues(1,1:length(freq_test)*k_freq_r),Leigenvalues(2,1:length(freq_test)*k_freq_r),k_freq_r,1,921)


% fifth case VSI
Yin_vil_pll_cal_frd = 1/Z_AFE_400_fi_1000_fv_100_pll_1_test_1;           % afe is the same
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
% disable ACindex plot
% figure(3)
% bode(ACindex);
% hold on
% bode(RRdq)
% title('ACindex and Ldd')
% grid on;
% Bode_Darklines(3);

% Eigenvalues and Sorting
RR_res=[RR(1,1).ResponseData RR(1,2).ResponseData;RR(2,1).ResponseData RR(2,2).ResponseData];
rr_res=[rr(1,1).ResponseData rr(1,2).ResponseData;rr(2,1).ResponseData rr(2,2).ResponseData];
for k=1:length(freq_test)*k_freq_r
%         Lddeigenvalues(:,k)=eig(Lddresp(:,:,k));
%         Lqqeigenvalues(:,k)=eig(Lqqresp(:,:,k));
    Leigenvalues(:,k)=eig(RR_res(:,:,k));
    leigenvalues(:,k)=eig(rr_res(:,:,k));
    Singular_L(:,k)=svd(RR_res(:,:,k));
    Singular_l(:,k)=svd(rr_res(:,:,k));
    [U H]=poldec(RR_res(:,:,k));
    eig_U(:,k)=eig(U);
    eig_H(:,k)=eig(H);
end
[Leigenvalues]=sortloci(Leigenvalues);
% plot the full GNC Characteristic Loci
figHandle=figure(93)
set(figHandle, 'Position', [50, 10, 1000, 800]);
plot(Leigenvalues(1,1:length(freq_test)*k_freq_r),linestyle1,'LineWidth',linewidth)
hold on
plot(conj(Leigenvalues(1,1:length(freq_test)*k_freq_r)),linestyle1,'LineWidth',linewidth)
plot(Leigenvalues(2,1:length(freq_test)*k_freq_r),linestyle2,'LineWidth',linewidth)
plot(conj(Leigenvalues(2,1:length(freq_test)*k_freq_r)),linestyle2,'LineWidth',linewidth)

legend({'{\it\lambda}_{1}','{\it\lambda}_{1}_{mirror}','{\it\lambda}_{2}','{\it\lambda}_{2}_{mirror}'},'Fontsize',fontsize,'FontWeight','bold')
plot(-1,0,'r+','LineWidth',linewidth)
% hold off
title('Characteristic Loci of Ldd','Fontsize',fontsize,'FontWeight','bold')
grid
set(gca,'FontSize',fontsize);
ezplot('x^2+y^2=1')
%% find the angle of eigenvalues
angle_lamda1_case2=angle(Leigenvalues(1,:))'*180/pi;
%% bode plot of eigenvalues
% lamda1_case2 = frd(Leigenvalues(1,1:length(freq_test)*k_freq_r),freq_test(1:length(freq_test)*k_freq_r));
% lamda2_case2 = frd(Leigenvalues(2,1:length(freq_test)*k_freq_r),freq_test(1:length(freq_test)*k_freq_r));
gershgorinplot(RR,Leigenvalues(1,1:length(freq_test)*k_freq_r),Leigenvalues(2,1:length(freq_test)*k_freq_r),k_freq_r,1,931)


% sixth case VSI
Yin_vil_pll_cal_frd = 1/Z_AFE_400_fi_1000_fv_50_pll_1_test_1;           % afe is the same
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
% disable ACindex plot
% figure(3)
% bode(ACindex);
% hold on
% bode(RRdq)
% title('ACindex and Ldd')
% grid on;
% Bode_Darklines(3);

% Eigenvalues and Sorting
RR_res=[RR(1,1).ResponseData RR(1,2).ResponseData;RR(2,1).ResponseData RR(2,2).ResponseData];
rr_res=[rr(1,1).ResponseData rr(1,2).ResponseData;rr(2,1).ResponseData rr(2,2).ResponseData];
for k=1:length(freq_test)*k_freq_r
%         Lddeigenvalues(:,k)=eig(Lddresp(:,:,k));
%         Lqqeigenvalues(:,k)=eig(Lqqresp(:,:,k));
    Leigenvalues(:,k)=eig(RR_res(:,:,k));
    leigenvalues(:,k)=eig(rr_res(:,:,k));
    Singular_L(:,k)=svd(RR_res(:,:,k));
    Singular_l(:,k)=svd(rr_res(:,:,k));
    [U H]=poldec(RR_res(:,:,k));
    eig_U(:,k)=eig(U);
    eig_H(:,k)=eig(H);
end
[Leigenvalues]=sortloci(Leigenvalues);
% plot the full GNC Characteristic Loci
figHandle=figure(94)
set(figHandle, 'Position', [50, 10, 1000, 800]);
plot(Leigenvalues(1,1:length(freq_test)*k_freq_r),linestyle1,'LineWidth',linewidth)
hold on
plot(conj(Leigenvalues(1,1:length(freq_test)*k_freq_r)),linestyle1,'LineWidth',linewidth)
plot(Leigenvalues(2,1:length(freq_test)*k_freq_r),linestyle2,'LineWidth',linewidth)
plot(conj(Leigenvalues(2,1:length(freq_test)*k_freq_r)),linestyle2,'LineWidth',linewidth)

legend({'{\it\lambda}_{1}','{\it\lambda}_{1}_{mirror}','{\it\lambda}_{2}','{\it\lambda}_{2}_{mirror}'},'Fontsize',fontsize,'FontWeight','bold')
plot(-1,0,'r+','LineWidth',linewidth)
% hold off
title('Characteristic Loci of Ldd','Fontsize',fontsize,'FontWeight','bold')
grid
set(gca,'FontSize',fontsize);
set(gca,'FontSize',fontsize);
ezplot('x^2+y^2=1')
%% find the angle of eigenvalues
angle_lamda1_case3=angle(Leigenvalues(1,:))'*180/pi;
%% bode plot of eigenvalues
% lamda1_case3 = frd(Leigenvalues(1,1:length(freq_test)*k_freq_r),freq_test(1:length(freq_test)*k_freq_r));
% lamda2_case3 = frd(Leigenvalues(2,1:length(freq_test)*k_freq_r),freq_test(1:length(freq_test)*k_freq_r));
gershgorinplot(RR,Leigenvalues(1,1:length(freq_test)*k_freq_r),Leigenvalues(2,1:length(freq_test)*k_freq_r),k_freq_r,1,941)

% figure(100)
% bode(lamda1_case1,lamda2_case1,Bode_O)
% figure(101)
% bode(lamda1_case2,lamda2_case2,Bode_O)
% figure(102)
% bode(lamda1_case3,lamda2_case3,Bode_O)

%% display the upper frquency for GNC plot
freq_test(length(freq_test)*k_freq_r)/(2*pi)
freq_test_check=freq_test/(2*pi);

% close all



%% plot pole and zero map to check whether to use GNC
% L_c1=Z_VSI_o_400_fi_4k_fv_80_avg/Z_AFE_400_fi_1000_fv_105_pll_1_avg;
% L_c2=Z_VSI_o_400_fi_4k_fv_110_avg/Z_AFE_400_fi_1000_fv_105_pll_1_avg;
% L_c3=Z_VSI_o_400_fi_4k_fv_150_avg/Z_AFE_400_fi_1000_fv_105_pll_1_avg;
% L_c4=Z_VSI_o_400_fi_4k_fv_200_avg/Z_AFE_400_fi_1000_fv_105_pll_1_avg;
% figure(4)
% pzmap(L_c1,L_c2,L_c3,L_c4)

% %% GNC plot using derived small signal model
% n= 1e3;
% w = logspace(-2,4,n);
% w = 2*pi*w;
% Yin_vil_pll_cal_frd = Y_AFE_400_fi_1000_fv_135_pll_1_cal;           % afe is the same
% 
% % first case VSI
% Zo_vsi_vil = Z_VSI_o_400_fi_4k_fv_200_cal;                           % then the unstable vsi case
% Zo_vsi_vil_s = [Zo_vsi_vil(1,1) 0; 0 Zo_vsi_vil(2,2)];
% Yin_vil_pll_cal_frd_s = [Yin_vil_pll_cal_frd(1,1) 0; 0 Yin_vil_pll_cal_frd(2,2)];
% % RR is the full GNC return ration
% RR=Zo_vsi_vil*Yin_vil_pll_cal_frd;
% % rr is the return ratio defined in Rolando's paper
% rr=Zo_vsi_vil_s*Yin_vil_pll_cal_frd_s;
% RRdd=rr(1,1);
% RRdq=rr(1,2);
% RRqd=rr(2,1);
% RRqq=rr(2,2);
% % calculate ACindex defined by Rolando
% ACindex=(RRdd-RRqq)*(RRdd-RRqq)/(4*RRqd*RRdq);
% % disable ACindex plot
% % figure(3)
% % bode(ACindex);
% % hold on
% % bode(RRdq)
% % title('ACindex and Ldd')
% % grid on;
% % Bode_Darklines(3);
% 
% % Eigenvalues and Sorting
% RR_res=[RR(1,1).ResponseData RR(1,2).ResponseData;RR(2,1).ResponseData RR(2,2).ResponseData];
% rr_res=[rr(1,1).ResponseData rr(1,2).ResponseData;rr(2,1).ResponseData rr(2,2).ResponseData];
% for k=1:length(w)*1
% %         Lddeigenvalues(:,k)=eig(Lddresp(:,:,k));
% %         Lqqeigenvalues(:,k)=eig(Lqqresp(:,:,k));
%     Leigenvalues(:,k)=eig(RR_res(:,:,k));
%     leigenvalues(:,k)=eig(rr_res(:,:,k));
% end
% % [Leigenvalues]=sortloci(Leigenvalues);
% % plot the full GNC Characteristic Loci
% figure(8)
% plot(leigenvalues(1,1:length(w)*0.8),linestyle1,'LineWidth',linewidth)
% hold on
% plot(leigenvalues(2,1:length(w)*0.8),linestyle2,'LineWidth',linewidth)
% legend({'{\it\lambda}_{1}','{\it\lambda}_{2}'},'Fontsize',fontsize,'FontWeight','bold')
% plot(-1,0,'r+','LineWidth',linewidth)
% % hold off
% title('Characteristic Loci of Ldd','Fontsize',fontsize,'FontWeight','bold')
% grid
% set(gca,'FontSize',fontsize);
% 
% % second case VSI
% Yin_vil_pll_cal_frd = Y_AFE_400_fi_1000_fv_105_pll_1_cal;           % afe is the same
% Yin_vil_pll_cal_frd_s = [Yin_vil_pll_cal_frd(1,1) 0; 0 Yin_vil_pll_cal_frd(2,2)];
% 
% RR=Zo_vsi_vil*Yin_vil_pll_cal_frd;
% rr=Zo_vsi_vil_s*Yin_vil_pll_cal_frd_s;
% RRdd=rr(1,1);
% RRdq=rr(1,2);
% RRqd=rr(2,1);
% RRqq=rr(2,2);
% 
% ACindex=(RRdd-RRqq)*(RRdd-RRqq)/(4*RRqd*RRdq);
% % figure(3)
% % bode(ACindex);
% % hold on
% % bode(RRdq)
% % title('ACindex and Ldd')
% % grid on;
% % Bode_Darklines(3);
% 
% % Eigenvalues and Sorting
% RR_res=[RR(1,1).ResponseData RR(1,2).ResponseData;RR(2,1).ResponseData RR(2,2).ResponseData];
% rr_res=[rr(1,1).ResponseData rr(1,2).ResponseData;rr(2,1).ResponseData rr(2,2).ResponseData];
% for k=1:length(w)*1
% %         Lddeigenvalues(:,k)=eig(Lddresp(:,:,k));
% %         Lqqeigenvalues(:,k)=eig(Lqqresp(:,:,k));
%     Leigenvalues(:,k)=eig(RR_res(:,:,k));
%     leigenvalues(:,k)=eig(rr_res(:,:,k));
% end
% % [Leigenvalues]=sortloci(Leigenvalues);
% % Characteristic Loci
% plot(leigenvalues(1,1:length(w)*0.8),linestyle2,'LineWidth',linewidth)
% % hold on
% plot(leigenvalues(2,1:length(w)*0.8),linestyle2,'LineWidth',linewidth)
% legend({'{\it\lambda}_{1}','{\it\lambda}_{2}'},'Fontsize',fontsize,'FontWeight','bold')
% plot(-1,0,'r+','LineWidth',linewidth)
% % hold off
% title('Characteristic Loci of Ldd','Fontsize',fontsize,'FontWeight','bold')
% grid
% set(gca,'FontSize',fontsize);
% 
% % third case VSI
% Yin_vil_pll_cal_frd = Y_AFE_400_fi_1000_fv_70_pll_1_cal;           % afe is the same
% Yin_vil_pll_cal_frd_s = [Yin_vil_pll_cal_frd(1,1) 0; 0 Yin_vil_pll_cal_frd(2,2)];
% 
% RR=Zo_vsi_vil*Yin_vil_pll_cal_frd;
% rr=Zo_vsi_vil_s*Yin_vil_pll_cal_frd_s;
% RRdd=rr(1,1);
% RRdq=rr(1,2);
% RRqd=rr(2,1);
% RRqq=rr(2,2);
% 
% ACindex=(RRdd-RRqq)*(RRdd-RRqq)/(4*RRqd*RRdq);
% % figure(3)
% % bode(ACindex);
% % hold on
% % bode(RRdq)
% % title('ACindex and Ldd')
% % grid on;
% % Bode_Darklines(3);
% 
% % Eigenvalues and Sorting
% RR_res=[RR(1,1).ResponseData RR(1,2).ResponseData;RR(2,1).ResponseData RR(2,2).ResponseData];
% rr_res=[rr(1,1).ResponseData rr(1,2).ResponseData;rr(2,1).ResponseData rr(2,2).ResponseData];
% for k=1:length(w)*1
% %         Lddeigenvalues(:,k)=eig(Lddresp(:,:,k));
% %         Lqqeigenvalues(:,k)=eig(Lqqresp(:,:,k));
%     Leigenvalues(:,k)=eig(RR_res(:,:,k));
%     leigenvalues(:,k)=eig(rr_res(:,:,k));
% end
% % [Leigenvalues]=sortloci(Leigenvalues);
% % Characteristic Loci
% plot(leigenvalues(1,1:length(w)*0.8),linestyle3,'LineWidth',linewidth)
% % hold on
% plot(leigenvalues(2,1:length(w)*0.8),linestyle2,'LineWidth',linewidth)
% legend({'{\it\lambda}_{1}','{\it\lambda}_{2}'},'Fontsize',fontsize,'FontWeight','bold')
% plot(-1,0,'r+','LineWidth',linewidth)
% % hold off
% title('Characteristic Loci of Ldd','Fontsize',fontsize,'FontWeight','bold')
% grid
% set(gca,'FontSize',fontsize);
% 
% % fourth case VSI
% Yin_vil_pll_cal_frd = Y_AFE_400_fi_1000_fv_35_pll_1_cal;           % afe is the same
% Yin_vil_pll_cal_frd_s = [Yin_vil_pll_cal_frd(1,1) 0; 0 Yin_vil_pll_cal_frd(2,2)];
% 
% RR=Zo_vsi_vil*Yin_vil_pll_cal_frd;
% rr=Zo_vsi_vil_s*Yin_vil_pll_cal_frd_s;
% RRdd=rr(1,1);
% RRdq=rr(1,2);
% RRqd=rr(2,1);
% RRqq=rr(2,2);
% 
% ACindex=(RRdd-RRqq)*(RRdd-RRqq)/(4*RRqd*RRdq);
% % figure(3)
% % bode(ACindex);
% % hold on
% % bode(RRdq)
% % title('ACindex and Ldd')
% % grid on;
% % Bode_Darklines(3);
% 
% % Eigenvalues and Sorting
% RR_res=[RR(1,1).ResponseData RR(1,2).ResponseData;RR(2,1).ResponseData RR(2,2).ResponseData];
% rr_res=[rr(1,1).ResponseData rr(1,2).ResponseData;rr(2,1).ResponseData rr(2,2).ResponseData];
% for k=1:length(w)*1
% %         Lddeigenvalues(:,k)=eig(Lddresp(:,:,k));
% %         Lqqeigenvalues(:,k)=eig(Lqqresp(:,:,k));
%     Leigenvalues(:,k)=eig(RR_res(:,:,k));
%     leigenvalues(:,k)=eig(rr_res(:,:,k));
% end
% % [Leigenvalues]=sortloci(Leigenvalues);
% % Characteristic Loci
% 
% plot(leigenvalues(1,1:length(w)*0.8),linestyle4,'LineWidth',linewidth)
% % hold on
% plot(leigenvalues(2,1:length(w)*0.8),linestyle2,'LineWidth',linewidth)
% legend({'{\it\lambda}_{1}','{\it\lambda}_{2}'},'Fontsize',fontsize,'FontWeight','bold')
% plot(-1,0,'r+','LineWidth',linewidth)
% hold off
% title('Characteristic Loci of Ldd','Fontsize',fontsize,'FontWeight','bold')
% grid on
% set(gca,'FontSize',fontsize);