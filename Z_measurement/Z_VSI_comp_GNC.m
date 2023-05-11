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
k_freq_r = 1;%1;%0.6;

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
% load('Z_VSI_400_fi_4000_fv_200_avg.mat')
% Z_VSI_o_400_fi_4k_fv_200_cal = Zo_vil_cal_frd;
% Z_VSI_o_400_fi_4k_fv_200_avg = Zo_vil_avg_sim;


load ('Z_afe_100p_vsi_fi_4k_fv200_afe_fi_1k_fv_70.mat')
Z_AFE_400_fi_1000_fv_70_pll_1_test = Z.ZDQRAW/3;


figure(1)
bode(Z_VSI_o_400_fi_4k_fv_250_test,Bode_O);
hold on
% bode(Z_VSI_o_400_fi_4k_fv_200_test,Bode_O);
% bode(Z_VSI_o_400_fi_4k_fv_150_test,Bode_O);
bode(Z_VSI_o_400_fi_4k_fv_110_test,Bode_O);
% bode(Z_VSI_o_400_fi_4k_fv_100_test,Bode_O);
bode(Z_VSI_o_400_fi_4k_fv_80_test,Bode_O);
bode(Z_AFE_400_fi_1000_fv_70_pll_1_test,Bode_O)
Bode_Darklines(3);



figure(3)
bode(Z_VSI_o_400_fi_4k_fv_250_test(1,1),Bode_O);
hold on
% bode(Z_VSI_o_400_fi_4k_fv_200_test(1,1),Bode_O);
% bode(Z_VSI_o_400_fi_4k_fv_150_test(1,1),Bode_O);
bode(Z_VSI_o_400_fi_4k_fv_110_test(1,1),Bode_O);
% bode(Z_VSI_o_400_fi_4k_fv_100_test(1,1),Bode_O);
bode(Z_VSI_o_400_fi_4k_fv_80_test(1,1),Bode_O);
bode(Z_AFE_400_fi_1000_fv_70_pll_1_test(1,1),Bode_O)
Bode_Darklines(1);


figure(31)
bode(Z_VSI_o_400_fi_4k_fv_250_test(1,1),Bode_O);
hold on
bode(Z_AFE_400_fi_1000_fv_70_pll_1_test(1,1),Bode_O)
Bode_Darklines(1);

figure(311)
bode(Z_VSI_o_400_fi_4k_fv_250_test(2,2),Bode_O);
hold on
bode(Z_AFE_400_fi_1000_fv_70_pll_1_test(2,2),Bode_O)
Bode_Darklines(1);

fighandle=figure(3111);
set(fighandle,'position',[10, 10, 1000, 800])
bode(Z_VSI_o_400_fi_4k_fv_250_test,Bode_O);
hold on
bode(Z_AFE_400_fi_1000_fv_70_pll_1_test,Bode_O)
Bode_Darklines(2);
legend ('Z\_VSI','Z\_AFE')

figure(32)
bode(Z_VSI_o_400_fi_4k_fv_150_test(1,1),Bode_O);
hold on
bode(Z_AFE_400_fi_1000_fv_70_pll_1_test(1,1),Bode_O)
Bode_Darklines(1);

figure(322)
bode(Z_VSI_o_400_fi_4k_fv_150_test(2,2),Bode_O);
hold on
bode(Z_AFE_400_fi_1000_fv_70_pll_1_test(2,2),Bode_O)
Bode_Darklines(1);

fighandle=figure(3222);
set(fighandle,'position',[10, 10, 1000, 800])
bode(Z_VSI_o_400_fi_4k_fv_150_test,Bode_O);
hold on
bode(Z_AFE_400_fi_1000_fv_70_pll_1_test,Bode_O)
Bode_Darklines(2);
legend ('Z\_VSI','Z\_AFE')

figure(33)
bode(Z_VSI_o_400_fi_4k_fv_80_test(1,1),Bode_O);
hold on
bode(Z_AFE_400_fi_1000_fv_70_pll_1_test(1,1),Bode_O)
Bode_Darklines(1);

figure(333)
bode(Z_VSI_o_400_fi_4k_fv_80_test(2,2),Bode_O);
hold on
bode(Z_AFE_400_fi_1000_fv_70_pll_1_test(2,2),Bode_O)
Bode_Darklines(1);

fighandle=figure(3333);
set(fighandle,'position',[10, 10, 1000, 800])
bode(Z_VSI_o_400_fi_4k_fv_80_test,Bode_O);
hold on
bode(Z_AFE_400_fi_1000_fv_70_pll_1_test,Bode_O)
Bode_Darklines(2);
legend ('Z\_VSI','Z\_AFE')
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
% Yin_vil_pll_cal_frd = Y_AFE_400_fi_1000_fv_105_pll_1_cal;           % afe is the same
% 
% % first case VSI
% Zo_vsi_vil = Z_VSI_o_400_fi_4k_fv_80_cal;                           % then the unstable vsi case
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
% figure(5)
% plot(Leigenvalues(1,1:length(w)*0.8),linestyle1,'LineWidth',linewidth)
% hold on
% plot(Leigenvalues(2,1:length(w)*0.8),linestyle2,'LineWidth',linewidth)
% legend({'{\it\lambda}_{1}','{\it\lambda}_{2}'},'Fontsize',fontsize,'FontWeight','bold')
% plot(-1,0,'r+','LineWidth',linewidth)
% % hold off
% title('Characteristic Loci of Ldd','Fontsize',fontsize,'FontWeight','bold')
% grid
% set(gca,'FontSize',fontsize);
% 
% % second case VSI
% Zo_vsi_vil = Z_VSI_o_400_fi_4k_fv_110_cal;
% Zo_vsi_vil_s = [Zo_vsi_vil(1,1) 0; 0 Zo_vsi_vil(2,2)];
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
% plot(Leigenvalues(1,1:length(w)*0.8),linestyle2,'LineWidth',linewidth)
% % hold on
% plot(Leigenvalues(2,1:length(w)*0.8),linestyle2,'LineWidth',linewidth)
% legend({'{\it\lambda}_{1}','{\it\lambda}_{2}'},'Fontsize',fontsize,'FontWeight','bold')
% plot(-1,0,'r+','LineWidth',linewidth)
% % hold off
% title('Characteristic Loci of Ldd','Fontsize',fontsize,'FontWeight','bold')
% grid
% set(gca,'FontSize',fontsize);
% 
% % third case VSI
% Zo_vsi_vil = Z_VSI_o_400_fi_4k_fv_150_cal;
% Zo_vsi_vil_s = [Zo_vsi_vil(1,1) 0; 0 Zo_vsi_vil(2,2)];
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
% plot(Leigenvalues(1,1:length(w)*0.8),linestyle3,'LineWidth',linewidth)
% % hold on
% plot(Leigenvalues(2,1:length(w)*0.8),linestyle2,'LineWidth',linewidth)
% legend({'{\it\lambda}_{1}','{\it\lambda}_{2}'},'Fontsize',fontsize,'FontWeight','bold')
% plot(-1,0,'r+','LineWidth',linewidth)
% % hold off
% title('Characteristic Loci of Ldd','Fontsize',fontsize,'FontWeight','bold')
% grid
% set(gca,'FontSize',fontsize);
% 
% % fourth case VSI
% Zo_vsi_vil = Z_VSI_o_400_fi_4k_fv_200_cal;
% Zo_vsi_vil_s = [Zo_vsi_vil(1,1) 0; 0 Zo_vsi_vil(2,2)];
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
% plot(Leigenvalues(1,1:length(w)*0.8),linestyle4,'LineWidth',linewidth)
% % hold on
% plot(Leigenvalues(2,1:length(w)*0.8),linestyle2,'LineWidth',linewidth)
% legend({'{\it\lambda}_{1}','{\it\lambda}_{2}'},'Fontsize',fontsize,'FontWeight','bold')
% plot(-1,0,'r+','LineWidth',linewidth)
% hold off
% title('Characteristic Loci of Ldd','Fontsize',fontsize,'FontWeight','bold')
% grid on
% set(gca,'FontSize',fontsize);


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
    principle_r1(:,k)=[max(eig_H(:,k))*exp(1i*max(angle(eig_U(:,k)))); max(eig_H(:,k))*exp(1i*min(angle(eig_U(:,k))))];
    principle_r2(:,k)=[min(eig_H(:,k))*exp(1i*min(angle(eig_U(:,k)))); min(eig_H(:,k))*exp(1i*max(angle(eig_U(:,k))))];
end
[Leigenvalues]=sortloci(Leigenvalues);
Mag_Leig=abs(Leigenvalues);
Arg_Leig=angle(Leigenvalues);
Arg_U=angle(eig_U);
[principle_r1]=sortloci(principle_r1);
[principle_r2]=sortloci(principle_r2);
% plot the full GNC Characteristic Loci
figHandle=figure(51)
set(figHandle, 'Position', [50, 10, 1000, 800]);
plot(Leigenvalues(1,1:length(freq_test)*k_freq_r),linestyle1,'LineWidth',linewidth)
hold on
plot(principle_r1(1,:),linestyle3)
plot(principle_r1(2,:),linestyle3)
plot(principle_r2(1,:),linestyle4)
plot(principle_r2(2,:),linestyle4)
% plot(conj(Leigenvalues(1,1:length(freq_test)*k_freq_r)),linestyle1,'LineWidth',linewidth)
plot(Leigenvalues(2,1:length(freq_test)*k_freq_r),linestyle2,'LineWidth',linewidth)
% plot(conj(Leigenvalues(2,1:length(freq_test)*k_freq_r)),linestyle2,'LineWidth',linewidth)

legend({'{\it\lambda}_{1}','{\it\lambda}_{1}_{mirror}','{\it\lambda}_{2}','{\it\lambda}_{2}_{mirror}'},'Fontsize',fontsize,'FontWeight','bold')
plot(-1,0,'r+','LineWidth',linewidth)
% hold off
title('Characteristic Loci of Ldd','Fontsize',fontsize,'FontWeight','bold')
grid
set(gca,'FontSize',fontsize);
ezplot('x^2+y^2=1')

figHandle=figure(511)
set(figHandle, 'Position', [50, 10, 1000, 800]);
plot(Mag_Leig(1,:))
hold on
plot(Mag_Leig(2,:))
plot(eig_H(1,:),'r')
plot(eig_H(2,:),'r')
plot(Singular_L(1,:),'g')
plot(Singular_L(2,:),'g')
figHandle=figure(5111)
set(figHandle, 'Position', [50, 10, 1000, 800]);
plot(Arg_Leig(1,:))
hold on
plot(Arg_Leig(2,:))
plot(Arg_U(1,:),'r')
plot(Arg_U(2,:),'r')

%% find the angle of eigenvalues
angle_lamda1_case1=360+angle(Leigenvalues(1,:))'*180/pi;
%% bode plot of eigenvalues
% lamda1_case1 = frd(Leigenvalues(1,:),freq_test*1:length(freq_test)*k_freq_r);
% lamda2_case1 = frd(Leigenvalues(2,:),freq_test*1:length(freq_test)*k_freq_r);
gershgorinplot(RR,k_freq_r,1,501)

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
    principle_r1(:,k)=[max(eig_H(:,k))*exp(1i*max(angle(eig_U(:,k)))); max(eig_H(:,k))*exp(1i*min(angle(eig_U(:,k))))];
    principle_r2(:,k)=[min(eig_H(:,k))*exp(1i*min(angle(eig_U(:,k)))); min(eig_H(:,k))*exp(1i*max(angle(eig_U(:,k))))];
end
[Leigenvalues]=sortloci(Leigenvalues);
Mag_Leig=abs(Leigenvalues);
Arg_Leig=angle(Leigenvalues);
Arg_U=angle(eig_U);
[principle_r1]=sortloci(principle_r1);
[principle_r2]=sortloci(principle_r2);
% plot the full GNC Characteristic Loci
figHandle=figure(52);
set(figHandle, 'Position', [50, 10, 1000, 800]);
plot(Leigenvalues(1,1:length(freq_test)*k_freq_r),linestyle1,'LineWidth',linewidth)
hold on
plot(principle_r1(1,:),linestyle3)
plot(principle_r1(2,:),linestyle3)
plot(principle_r2(1,:),linestyle4)
plot(principle_r2(2,:),linestyle4)
% plot(conj(Leigenvalues(1,1:length(freq_test)*k_freq_r)),linestyle1,'LineWidth',linewidth)
plot(Leigenvalues(2,1:length(freq_test)*k_freq_r),linestyle2,'LineWidth',linewidth)
% plot(conj(Leigenvalues(2,1:length(freq_test)*k_freq_r)),linestyle2,'LineWidth',linewidth)

legend({'{\it\lambda}_{1}','{\it\lambda}_{1}_{mirror}','{\it\lambda}_{2}','{\it\lambda}_{2}_{mirror}'},'Fontsize',fontsize,'FontWeight','bold')
plot(-1,0,'r+','LineWidth',linewidth)
% hold off
title('Characteristic Loci of Ldd','Fontsize',fontsize,'FontWeight','bold')
grid
set(gca,'FontSize',fontsize);
set(gca,'FontSize',fontsize);
ezplot('x^2+y^2=1')

figHandle=figure(522)
set(figHandle, 'Position', [50, 10, 1000, 800]);
plot(Mag_Leig(1,:))
hold on
plot(Mag_Leig(2,:))
plot(eig_H(1,:),'r')
plot(eig_H(2,:),'r')
plot(Singular_L(1,:),'g')
plot(Singular_L(2,:),'g')
figHandle=figure(5222)
set(figHandle, 'Position', [50, 10, 1000, 800]);
plot(Arg_Leig(1,:))
hold on
plot(Arg_Leig(2,:))
plot(Arg_U(1,:),'r')
plot(Arg_U(2,:),'r')


%% find the angle of eigenvalues
angle_lamda1_case2=angle(Leigenvalues(1,:))'*180/pi;
%% bode plot of eigenvalues
% lamda1_case2 = frd(Leigenvalues(1,:),freq_test);
% lamda2_case2 = frd(Leigenvalues(2,:),freq_test);
gershgorinplot(RR,k_freq_r,1,521)

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
% disable ACindex plot
% figure(3)
% bode(ACindex);
% hold on
% bode(RRdq)
% title('ACindex and Ldd')
% grid on;
% Bode_Darklines(3);
clear Leigenvalues leigenvalues Singular_L Singular_l U H eig_U eig_H principle_r1 principle_r2
% Eigenvalues and Sorting
RR_res=[RR(1,1).ResponseData RR(1,2).ResponseData;RR(2,1).ResponseData RR(2,2).ResponseData];
rr_res=[rr(1,1).ResponseData rr(1,2).ResponseData;rr(2,1).ResponseData rr(2,2).ResponseData];
% [RR_res]=sortloci(RR_res);
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
    principle_r1(:,k)=[max(eig_H(:,k))*exp(1i*max(angle(eig_U(:,k)))); max(eig_H(:,k))*exp(1i*min(angle(eig_U(:,k))))];
    principle_r2(:,k)=[min(eig_H(:,k))*exp(1i*min(angle(eig_U(:,k)))); min(eig_H(:,k))*exp(1i*max(angle(eig_U(:,k))))];
end
[Leigenvalues]=sortloci(Leigenvalues);
Mag_Leig=abs(Leigenvalues);
Arg_Leig=angle(Leigenvalues);
Arg_U=angle(eig_U);
[principle_r1]=sortloci(principle_r1);
[principle_r2]=sortloci(principle_r2);
% plot the full GNC Characteristic Loci
figHandle=figure(53);
set(figHandle, 'Position', [50, 10, 1000, 800]);
plot(Leigenvalues(1,:),linestyle1,'LineWidth',linewidth)
hold on
plot(principle_r1(1,:),linestyle3)
plot(principle_r1(2,:),linestyle3)
plot(principle_r2(1,:),linestyle4)
plot(principle_r2(2,:),linestyle4)

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

figHandle=figure(533)
set(figHandle, 'Position', [50, 10, 1000, 800]);
plot(Mag_Leig(1,:))
hold on
plot(Mag_Leig(2,:))
plot(eig_H(1,:),'r')
plot(eig_H(2,:),'r')
plot(Singular_L(1,:),'g')
plot(Singular_L(2,:),'g')
figHandle=figure(5333)
set(figHandle, 'Position', [50, 10, 1000, 800]);
plot(Arg_Leig(1,:))
hold on
plot(Arg_Leig(2,:))
plot(Arg_U(1,:),'r')
plot(Arg_U(2,:),'r')
%% find the angle of eigenvalues
angle_lamda1_case3=angle(Leigenvalues(1,:))'*180/pi;
%% bode plot of eigenvalues
% lamda1_case3 = frd(Leigenvalues(1,:),freq_test);
% lamda2_case3 = frd(Leigenvalues(2,:),freq_test);
gershgorinplot(RR,k_freq_r,1,531)
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