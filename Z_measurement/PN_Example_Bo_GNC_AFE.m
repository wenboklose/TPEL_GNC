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
k_freq_r = 1;% 0.6;


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


k_d = +5;
k_g = 1.1;

% % first case VSI
Zo_vsi_vil = Z_VSI_o_400_fi_4k_fv_200_test_1;                           % then the unstable vsi case

% second case VSI
Yin_vil_pll_cal_frd = 1/Z_AFE_400_fi_1000_fv_220_pll_1_test_1;           % afe is the same

AZ = 1/sqrt(2)*[1 1i;1 -1i];
Z_S_pn = AZ*Zo_vsi_vil*inv(AZ);
Y_L_pn = AZ*Yin_vil_pll_cal_frd*inv(AZ);

% RR is the full GNC return ration
RR=Zo_vsi_vil*Yin_vil_pll_cal_frd;

RR_pn=Z_S_pn*Y_L_pn;

% Eigenvalues and Sorting
RR_res=[RR(1,1).ResponseData RR(1,2).ResponseData;RR(2,1).ResponseData RR(2,2).ResponseData];
RR_pn_res=[RR_pn(1,1).ResponseData RR_pn(1,2).ResponseData; RR_pn(2,1).ResponseData RR_pn(2,2).ResponseData];

for k=1:length(freq_test)*k_freq_r
%         Lddeigenvalues(:,k)=eig(Lddresp(:,:,k));
%         Lqqeigenvalues(:,k)=eig(Lqqresp(:,:,k));
    Leigenvalues(:,k)=eig(RR_res(:,:,k));
    Leigenvalues_pn(:,k)=eig(RR_pn_res(:,:,k));
    
end
[Leigenvalues]=sortloci(Leigenvalues);
[Leigenvalues_pn]=sortloci(Leigenvalues_pn);
% plot the full GNC Characteristic Loci
figHandle=figure(2);
set(figHandle, 'Position', [50, 10, 1000, 800]);
plot(Leigenvalues(1,:),linestyle1,'LineWidth',linewidth)
hold on
% plot(conj(Leigenvalues(1,1:length(freq_test)*k_freq_r)),linestyle1,'LineWidth',linewidth)
plot(Leigenvalues(2,:),linestyle2,'LineWidth',linewidth)
plot(Leigenvalues_pn(1,:),linestyle3,'LineWidth',linewidth)
hold on
plot(Leigenvalues_pn(2,:),linestyle4,'LineWidth',linewidth)
% plot(conj(Leigenvalues(2,1:length(freq_test)*k_freq_r)),linestyle2,'LineWidth',linewidth)

% legend({'{\it\lambda}_{1}','{\it\lambda}_{1}_{mirror}','{\it\lambda}_{2}','{\it\lambda}_{2}_{mirror}'},'Fontsize',fontsize,'FontWeight','bold')
plot(-1,0,'r+','LineWidth',linewidth)
% hold off
title('Characteristic Loci of Ldd','Fontsize',fontsize,'FontWeight','bold')
grid
set(gca,'FontSize',fontsize);
ezplot('x^2+y^2=1')


% fifth case VSI
Yin_vil_pll_cal_frd = 1/Z_AFE_400_fi_1000_fv_100_pll_1_test_1;           % afe is the same
AZ = 1/sqrt(2)*[1 1i;1 -1i];
Y_L_pn = AZ*Yin_vil_pll_cal_frd*inv(AZ);

% RR is the full GNC return ration
RR=Zo_vsi_vil*Yin_vil_pll_cal_frd;
RR_pn=Z_S_pn*Y_L_pn;
% Eigenvalues and Sorting
RR_res=[RR(1,1).ResponseData RR(1,2).ResponseData;RR(2,1).ResponseData RR(2,2).ResponseData];
RR_pn_res=[RR_pn(1,1).ResponseData RR_pn(1,2).ResponseData; RR_pn(2,1).ResponseData RR_pn(2,2).ResponseData];

for k=1:length(freq_test)*k_freq_r
%         Lddeigenvalues(:,k)=eig(Lddresp(:,:,k));
%         Lqqeigenvalues(:,k)=eig(Lqqresp(:,:,k));
    Leigenvalues(:,k)=eig(RR_res(:,:,k));
    Leigenvalues_pn(:,k)=eig(RR_pn_res(:,:,k));

end
[Leigenvalues]=sortloci(Leigenvalues);
[Leigenvalues_pn]=sortloci(Leigenvalues_pn);

% plot the full GNC Characteristic Loci
figHandle=figure(21);
set(figHandle, 'Position', [50, 10, 1000, 800]);
plot(Leigenvalues(1,:),linestyle1,'LineWidth',linewidth)
hold on

% plot(conj(Leigenvalues(1,1:length(freq_test)*k_freq_r)),linestyle1,'LineWidth',linewidth)
plot(Leigenvalues(2,:),linestyle2,'LineWidth',linewidth)
% plot(conj(Leigenvalues(2,1:length(freq_test)*k_freq_r)),linestyle2,'LineWidth',linewidth)
plot(Leigenvalues_pn(1,:),linestyle3,'LineWidth',linewidth)
hold on
plot(Leigenvalues_pn(2,:),linestyle4,'LineWidth',linewidth)

legend({'{\it\lambda}_{1}','{\it\lambda}_{1}_{mirror}','{\it\lambda}_{2}','{\it\lambda}_{2}_{mirror}'},'Fontsize',fontsize,'FontWeight','bold')
plot(-1,0,'r+','LineWidth',linewidth)
% hold off
title('Characteristic Loci of Ldd','Fontsize',fontsize,'FontWeight','bold')
grid
set(gca,'FontSize',fontsize);
ezplot('x^2+y^2=1')


% sixth case VSI
Yin_vil_pll_cal_frd = 1/Z_AFE_400_fi_1000_fv_50_pll_1_test_1;           % afe is the same
AZ = 1/sqrt(2)*[1 1i;1 -1i];
Y_L_pn = AZ*Yin_vil_pll_cal_frd*inv(AZ);

% RR is the full GNC return ration
RR=Zo_vsi_vil*Yin_vil_pll_cal_frd;
RR_pn=Z_S_pn*Y_L_pn;

% Eigenvalues and Sorting
RR_res=[RR(1,1).ResponseData RR(1,2).ResponseData;RR(2,1).ResponseData RR(2,2).ResponseData];
RR_pn_res=[RR_pn(1,1).ResponseData RR_pn(1,2).ResponseData; RR_pn(2,1).ResponseData RR_pn(2,2).ResponseData];
for k=1:length(freq_test)*k_freq_r
%         Lddeigenvalues(:,k)=eig(Lddresp(:,:,k));
%         Lqqeigenvalues(:,k)=eig(Lqqresp(:,:,k));
    Leigenvalues(:,k)=eig(RR_res(:,:,k));
    Leigenvalues_pn(:,k)=eig(RR_pn_res(:,:,k));
end
[Leigenvalues]=sortloci(Leigenvalues);
[Leigenvalues_pn]=sortloci(Leigenvalues_pn);

% plot the full GNC Characteristic Loci
figHandle=figure(22);
set(figHandle, 'Position', [50, 10, 1000, 800]);
plot(Leigenvalues(1,:),linestyle1,'LineWidth',linewidth)
hold on

% plot(conj(Leigenvalues(1,1:length(freq_test)*k_freq_r)),linestyle1,'LineWidth',linewidth)
plot(Leigenvalues(2,:),linestyle2,'LineWidth',linewidth)
% plot(conj(Leigenvalues(2,1:length(freq_test)*k_freq_r)),linestyle2,'LineWidth',linewidth)
plot(Leigenvalues_pn(1,:),linestyle3,'LineWidth',linewidth)
hold on
plot(Leigenvalues_pn(2,:),linestyle4,'LineWidth',linewidth)

legend({'{\it\lambda}_{1}','{\it\lambda}_{1}_{mirror}','{\it\lambda}_{2}','{\it\lambda}_{2}_{mirror}'},'Fontsize',fontsize,'FontWeight','bold')
plot(-1,0,'r+','LineWidth',linewidth)
% hold off
title('Characteristic Loci of Ldd','Fontsize',fontsize,'FontWeight','bold')
grid
set(gca,'FontSize',fontsize);
ezplot('x^2+y^2=1')
