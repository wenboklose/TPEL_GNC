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
linewidth=2; fontsize=18;
linestyle1='-rs'; linestyle2='-ks';
linestyle3='r'; linestyle4='k';
linestyle5='m'; linestyle6='c';

% linestyle=[linestyle1; linestyle2; linestyle3; linestyle4; linestyle5; linestyle6];

% define GNC plot frequency range
k_freq_r = 1;%0.86/4;

load ('Zin_AFE_G1_exp.mat')
Zin_AFE_G1_exp = Z.ZDQ/3;%Z.ZDQRAW/3;%
load ('Zout_VSI_G1_c1_exp.mat')
Zout_VSI_G1_c1_exp = Z.ZDQ/3;%Z.ZDQRAW/3;%
load ('Zout_VSI_G1_c2_exp.mat')
Zout_VSI_G1_c2_exp = Z.ZDQ/3;%Z.ZDQRAW/3;%
load ('Zout_VSI_G1_c3_exp.mat')
Zout_VSI_G1_c3_exp = Z.ZDQ/3;%Z.ZDQRAW/3;%

%% load average model simulation results:

load ('Zin_AFE_G1_avg.mat')
Zin_AFE_G1_avg = Zin_vl_pll_avg_sim;
load ('Zout_VSI_G1_c1_avg_v1.mat')
Zout_VSI_G1_c1_avg = Zo_vil_avg_sim;%Z.ZDQRAW/3;%
load ('Zout_VSI_G1_c2_avg_v1.mat')
Zout_VSI_G1_c2_avg = Zo_vil_avg_sim;%Z.ZDQRAW/3;%
load ('Zout_VSI_G1_c3_avg_v1.mat')
Zout_VSI_G1_c3_avg = Zo_vil_avg_sim;%Z.ZDQRAW/3;%

%% GNC plot using measured exprimental results
freq_exp = Zin_AFE_G1_exp(1,1).Frequency;

Z_dd_resp = Zin_AFE_G1_exp(1,1).ResponseData;
Z_dd = frd(Z_dd_resp,freq_exp);
Z_dq_resp = Zin_AFE_G1_exp(1,2).ResponseData;
Z_dq = frd(Z_dq_resp,freq_exp);
Z_qd_resp = Zin_AFE_G1_exp(2,1).ResponseData;
Z_qd = frd(Z_qd_resp,freq_exp);
Z_qq_resp = Zin_AFE_G1_exp(2,2).ResponseData;
Z_qq = frd(Z_qq_resp,freq_exp);
Zin_AFE_G1_exp = [Z_dd Z_dq; Z_qd Z_qq];
Yin_AFE_G1_exp = 1/Zin_AFE_G1_exp;           % afe is the same

% reconstruct every test data:
Z_dd_resp = Zout_VSI_G1_c1_exp(1,1).ResponseData;
Z_dd = frd(Z_dd_resp,freq_exp);
Z_dq_resp = Zout_VSI_G1_c1_exp(1,2).ResponseData;
Z_dq = frd(Z_dq_resp,freq_exp);
Z_qd_resp = Zout_VSI_G1_c1_exp(2,1).ResponseData;
Z_qd = frd(Z_qd_resp,freq_exp);
Z_qq_resp = Zout_VSI_G1_c1_exp(2,2).ResponseData;
Z_qq = frd(Z_qq_resp,freq_exp);
Zout_VSI_G1_c1_exp = [Z_dd Z_dq; Z_qd Z_qq];

Z_dd_resp = Zout_VSI_G1_c2_exp(1,1).ResponseData;
Z_dd = frd(Z_dd_resp,freq_exp);
Z_dq_resp = Zout_VSI_G1_c2_exp(1,2).ResponseData;
Z_dq = frd(Z_dq_resp,freq_exp);
Z_qd_resp = Zout_VSI_G1_c2_exp(2,1).ResponseData;
Z_qd = frd(Z_qd_resp,freq_exp);
Z_qq_resp = Zout_VSI_G1_c2_exp(2,2).ResponseData;
Z_qq = frd(Z_qq_resp,freq_exp);
Zout_VSI_G1_c2_exp = [Z_dd Z_dq; Z_qd Z_qq];

Z_dd_resp = Zout_VSI_G1_c3_exp(1,1).ResponseData;
Z_dd = frd(Z_dd_resp,freq_exp);
Z_dq_resp = Zout_VSI_G1_c3_exp(1,2).ResponseData;
Z_dq = frd(Z_dq_resp,freq_exp);
Z_qd_resp = Zout_VSI_G1_c3_exp(2,1).ResponseData;
Z_qd = frd(Z_qd_resp,freq_exp);
Z_qq_resp = Zout_VSI_G1_c3_exp(2,2).ResponseData;
Z_qq = frd(Z_qq_resp,freq_exp);
Zout_VSI_G1_c3_exp = [Z_dd Z_dq; Z_qd Z_qq];

% Zin_AFE for group 1
Yin_vil_pll_exp = Yin_AFE_G1_exp;           % afe is the same
% case 1
Zo_vsi_vil = Zout_VSI_G1_c1_exp;                           % then the unstable vsi case
% RR is the full GNC return ration
RR=Zo_vsi_vil*Yin_vil_pll_exp;

RR_res=[RR(1,1).ResponseData RR(1,2).ResponseData;RR(2,1).ResponseData RR(2,2).ResponseData];
% find loci
for k=1:length(freq_exp)*k_freq_r
    Leigenvalues(:,k)=eig(RR_res(:,:,k));
end
[Leigenvalues]=sortloci(Leigenvalues);
% plot the full GNC Characteristic Loci
figHandle=figure(1);
set(figHandle, 'Position', [50, 10, 1000, 800]);
plot(Leigenvalues(1,:),linestyle1,'LineWidth',linewidth)
hold on
plot(Leigenvalues(2,:),linestyle2,'LineWidth',linewidth)

n=400;
f=logspace(1.60206,4,n);
% f=logspace(1,4,n);
freq_avg=f*2*pi;
Zo_vsi_vil = Zout_VSI_G1_c1_avg;
Yin_vil_pll_avg = 1/Zin_AFE_G1_avg;
PM = -5;
PM1 = 10;
GM = 0.980;
U=[exp(-1i*PM*pi/180) 0;0 exp(-1i*PM*pi/180)];
R=[GM 0;0 GM];
% RR is the full GNC return ration
clear RR
% RR is the full GNC return ration
RR=Zo_vsi_vil*Yin_vil_pll_avg;
% % closed-loop transfer function for the system
% I = [tf([0 1],[0 1]) tf([0 0],[0 1]);tf([0 0],[0 1]) tf([0 1],[0 1])];
% tf_cl_c1=I/(I+RR);
% tf_cl=tf_cl_c1;
% % simulation of the interface response
% sim('Vdq_response');
% Vdq_c1_sim = Vdq_sim;
clear RR_res
RR_res1=freqresp(RR,freq_avg);

for k=1:length(freq_avg)*k_freq_r
    Leigenvalues1(:,k)=eig(RR_res1(:,:,k));
end
[Leigenvalues1]=sortloci(Leigenvalues1);
Leigenvalues1(1,:)=Leigenvalues1(1,:)*1*exp(-1i*PM*pi/180);
Leigenvalues1(2,:)=Leigenvalues1(2,:)*GM/1*exp(-1i*PM1*pi/180);
set(figHandle, 'Position', [50, 10, 1000, 800]);
plot(Leigenvalues1(1,:),linestyle3,'LineWidth',linewidth)
plot(Leigenvalues1(2,:),linestyle4,'LineWidth',linewidth)
legend({'{\it\lambda}_{1\_case1\_exp}','{\it\lambda}_{2\_case1\_exp}','{\it\lambda}_{1\_case1\_avg}','{\it\lambda}_{2\_case1\_avg}'},'Fontsize',fontsize,'FontWeight','bold')
title('Characteristic Loci of Ldd','Fontsize',fontsize,'FontWeight','bold')
plot(-1,0,'g+','LineWidth',3)
grid
set(gca,'FontSize',fontsize);
ezplot('x^2+y^2=1')

% % case 2
% Zo_vsi_vil = Zout_VSI_G1_c2_exp;                           % then the unstable vsi case
% 
% % RR is the full GNC return ration
% RR=Zo_vsi_vil*Yin_vil_pll_exp;
% 
% clear RR_res
% RR_res=[RR(1,1).ResponseData RR(1,2).ResponseData;RR(2,1).ResponseData RR(2,2).ResponseData];
% 
% for k=1:length(freq_exp)*k_freq_r
%     Leigenvalues(:,k)=eig(RR_res(:,:,k));
% end
% [Leigenvalues]=sortloci(Leigenvalues);
% 
% % plot the full GNC Characteristic Loci
% figHandle=figure(2);
% set(figHandle, 'Position', [50, 10, 1000, 800]);
% plot(Leigenvalues(1,:),linestyle1,'LineWidth',linewidth)
% hold on
% plot(Leigenvalues(2,:),linestyle2,'LineWidth',linewidth)
% 
% Zo_vsi_vil = Zout_VSI_G1_c2_avg;                           % then the unstable vsi case
% PM = 1;
% PM1 = 10;
% GM = 0.98;
% U=[exp(-1i*PM*pi/180) 0;0 exp(-1i*PM1*pi/180)];
% R=[GM 0;0 GM];
% % RR is the full GNC return ration
% RR=Zo_vsi_vil*Yin_vil_pll_avg;
% % % closed-loop transfer function for the system
% % I = [tf([0 1],[0 1]) tf([0 0],[0 1]);tf([0 0],[0 1]) tf([0 1],[0 1])];
% % tf_cl_c2=I/(I+RR);
% % tf_cl=tf_cl_c2;
% % % simulation of the interface response
% % sim('Vdq_response');
% % Vdq_c2_sim = Vdq_sim;
% clear RR_res
% RR_res1=freqresp(RR,freq_avg);
% % [RR(1,1).ResponseData RR(1,2).ResponseData;RR(2,1).ResponseData RR(2,2).ResponseData];
% 
% for k=1:length(freq_avg)*k_freq_r
%     Leigenvalues2(:,k)=eig(RR_res1(:,:,k));
% end
% [Leigenvalues2]=sortloci(Leigenvalues2);
% Leigenvalues2(1,:)=Leigenvalues2(1,:)*GM*exp(-1i*PM*pi/180);
% Leigenvalues2(2,:)=Leigenvalues2(2,:)*GM*exp(-1i*PM1*pi/180);
% plot(Leigenvalues2(1,:),linestyle3,'LineWidth',linewidth)
% plot(Leigenvalues2(2,:),linestyle4,'LineWidth',linewidth)
% legend({'{\it\lambda}_{1\_case2\_exp}','{\it\lambda}_{2\_case2\_exp}','{\it\lambda}_{1\_case2\_avg}','{\it\lambda}_{2\_case2\_avg}'},'Fontsize',fontsize,'FontWeight','bold')
% title('Characteristic Loci of Ldd','Fontsize',fontsize,'FontWeight','bold')
% plot(-1,0,'g+','LineWidth',3)
% grid
% set(gca,'FontSize',fontsize);
% ezplot('x^2+y^2=1')
% 
% % case 3
% Zo_vsi_vil = Zout_VSI_G1_c3_exp;                           % then the unstable vsi case
% 
% % RR is the full GNC return ration
% RR=Zo_vsi_vil*Yin_vil_pll_exp;
% 
% clear RR_res
% RR_res=[RR(1,1).ResponseData RR(1,2).ResponseData;RR(2,1).ResponseData RR(2,2).ResponseData];
% 
% for k=1:length(freq_exp)*k_freq_r
%     Leigenvalues(:,k)=eig(RR_res(:,:,k));
% end
% [Leigenvalues]=sortloci(Leigenvalues);
% 
% % plot the full GNC Characteristic Loci
% figHandle=figure(3);
% set(figHandle, 'Position', [50, 10, 1000, 800]);
% plot(Leigenvalues(1,:),linestyle1,'LineWidth',linewidth)
% hold on
% plot(Leigenvalues(2,:),linestyle2,'LineWidth',linewidth)
% 
% Zo_vsi_vil = Zout_VSI_G1_c3_avg;                           % then the unstable vsi case
% PM = 3;
% PM1 = 10;
% GM = 0.90;
% U=[exp(-1i*PM*pi/180) 0;0 exp(-1i*PM1*pi/180)];
% R=[GM 0;0 GM];
% % RR is the full GNC return ration
% RR=Zo_vsi_vil*Yin_vil_pll_avg;
% % % closed-loop transfer function for the system
% % I = [tf([0 1],[0 1]) tf([0 0],[0 1]);tf([0 0],[0 1]) tf([0 1],[0 1])];
% % tf_cl_c3=I/(I+RR);
% % % simulation of the interface response
% % tf_cl=tf_cl_c3;
% % sim('Vdq_response');
% % Vdq_c3_sim = Vdq_sim;
% clear RR_res
% RR_res1=freqresp(RR,freq_avg);
% for k=1:length(freq_avg)*k_freq_r
%     Leigenvalues3(:,k)=eig(RR_res1(:,:,k));
% end
% [Leigenvalues3]=sortloci(Leigenvalues3);
% Leigenvalues3(1,:)=Leigenvalues3(1,:)*GM*exp(-1i*PM*pi/180);
% Leigenvalues3(2,:)=Leigenvalues3(2,:)*GM*exp(-1i*PM1*pi/180);
% plot(Leigenvalues3(1,:),linestyle3,'LineWidth',linewidth)
% plot(Leigenvalues3(2,:),linestyle4,'LineWidth',linewidth)
% legend({'{\it\lambda}_{1\_case3\_exp}','{\it\lambda}_{2\_case3\_exp}','{\it\lambda}_{1\_case3\_avg}','{\it\lambda}_{2\_case3\_avg}'},'Fontsize',fontsize,'FontWeight','bold')
% title('Characteristic Loci of Ldd','Fontsize',fontsize,'FontWeight','bold')
% plot(-1,0,'g+','LineWidth',3)
% grid
% set(gca,'FontSize',fontsize);
% ezplot('x^2+y^2=1')

% 
% %% using simulation data:
% %% GNC plot using measured STASU simulation results
% 
% % freq_avg = freq_exp;
% % 
% % n=400;
% % f=logspace(1,4,n);
% % % f=40:.05:1e3;
% % freq_avg=f*2*pi;
%     
% Zo_vsi_vil = Zout_VSI_G1_c1_avg;
% Yin_vil_pll_exp = 1/Zin_AFE_G1_avg;
% 
% % RR is the full GNC return ration
% RR=Zo_vsi_vil*Yin_vil_pll_exp;
% 
% RR_res=freqresp(RR,freq_avg);
% % [RR(1,1).ResponseData RR(1,2).ResponseData;RR(2,1).ResponseData RR(2,2).ResponseData];
% 
% for k=1:length(freq_avg)*k_freq_r
% %         Lddeigenvalues(:,k)=eig(Lddresp(:,:,k));
% %         Lqqeigenvalues(:,k)=eig(Lqqresp(:,:,k));
%     Leigenvalues(:,k)=eig(RR_res(:,:,k));
% 
% end
% [Leigenvalues]=sortloci(Leigenvalues);
% 
% 
% % plot the full GNC Characteristic Loci
% figHandle=figure(11);
% set(figHandle, 'Position', [50, 10, 1000, 800]);
% plot(Leigenvalues(1,:),linestyle1,'LineWidth',linewidth)
% hold on
% 
% plot(Leigenvalues(2,:),linestyle2,'LineWidth',linewidth)
% 
% legend({'{\it\lambda}_{1}','{\it\lambda}_{2}'},'Fontsize',fontsize,'FontWeight','bold')
% plot(-1,0,'r+','LineWidth',linewidth)
% % hold off
% title('Characteristic Loci of Ldd','Fontsize',fontsize,'FontWeight','bold')
% grid
% set(gca,'FontSize',fontsize);
% ezplot('x^2+y^2=1')
% 
% % case 2
% Zo_vsi_vil = Zout_VSI_G1_c2_avg;                           % then the unstable vsi case
% 
% % RR is the full GNC return ration
% RR=Zo_vsi_vil*Yin_vil_pll_exp;
% 
% RR_res=freqresp(RR,freq_avg);
% % [RR(1,1).ResponseData RR(1,2).ResponseData;RR(2,1).ResponseData RR(2,2).ResponseData];
% 
% for k=1:length(freq_avg)*k_freq_r
% %         Lddeigenvalues(:,k)=eig(Lddresp(:,:,k));
% %         Lqqeigenvalues(:,k)=eig(Lqqresp(:,:,k));
%     Leigenvalues(:,k)=eig(RR_res(:,:,k));
% 
% end
% [Leigenvalues]=sortloci(Leigenvalues);
% 
% 
% % plot the full GNC Characteristic Loci
% figHandle=figure(21);
% set(figHandle, 'Position', [50, 10, 1000, 800]);
% plot(Leigenvalues(1,:),linestyle1,'LineWidth',linewidth)
% hold on
% 
% plot(Leigenvalues(2,:),linestyle2,'LineWidth',linewidth)
% 
% legend({'{\it\lambda}_{1}','{\it\lambda}_{2}'},'Fontsize',fontsize,'FontWeight','bold')
% plot(-1,0,'r+','LineWidth',linewidth)
% % hold off
% title('Characteristic Loci of Ldd','Fontsize',fontsize,'FontWeight','bold')
% grid
% set(gca,'FontSize',fontsize);
% ezplot('x^2+y^2=1')
% 
% 
% % case 3
% Zo_vsi_vil = Zout_VSI_G1_c3_avg;                           % then the unstable vsi case
% 
% % RR is the full GNC return ration
% RR=Zo_vsi_vil*Yin_vil_pll_exp;
% 
% 
% RR_res=freqresp(RR,freq_avg);
% % [RR(1,1).ResponseData RR(1,2).ResponseData;RR(2,1).ResponseData RR(2,2).ResponseData];
% 
% for k=1:length(freq_avg)*k_freq_r
% %         Lddeigenvalues(:,k)=eig(Lddresp(:,:,k));
% %         Lqqeigenvalues(:,k)=eig(Lqqresp(:,:,k));
%     Leigenvalues(:,k)=eig(RR_res(:,:,k));
% 
% end
% [Leigenvalues]=sortloci(Leigenvalues);
% 
% 
% % plot the full GNC Characteristic Loci
% figHandle=figure(31);
% set(figHandle, 'Position', [50, 10, 1000, 800]);
% plot(Leigenvalues(1,:),linestyle1,'LineWidth',linewidth)
% hold on
% 
% plot(Leigenvalues(2,:),linestyle2,'LineWidth',linewidth)
% 
% legend({'{\it\lambda}_{1}','{\it\lambda}_{2}'},'Fontsize',fontsize,'FontWeight','bold')
% plot(-1,0,'r+','LineWidth',linewidth)
% % hold off
% title('Characteristic Loci of Ldd','Fontsize',fontsize,'FontWeight','bold')
% grid
% set(gca,'FontSize',fontsize);
% ezplot('x^2+y^2=1')
% 
% bode plots for all the impedances

figHandle=figure(12);
set(figHandle, 'Position', [50, 10, 1000, 800]);
bode(Zout_VSI_G1_c1_exp,Bode_O)
hold on
bode(Zin_AFE_G1_exp,Bode_O)
Bode_Darklines(2)

figHandle=figure(22);
set(figHandle, 'Position', [50, 10, 1000, 800]);
bode(Zout_VSI_G1_c2_exp,Bode_O)
hold on
bode(Zin_AFE_G1_exp,Bode_O)
Bode_Darklines(2)

figHandle=figure(32);
set(figHandle, 'Position', [50, 10, 1000, 800]);
bode(Zout_VSI_G1_c3_exp,Bode_O)
hold on
bode(Zin_AFE_G1_exp,Bode_O)
Bode_Darklines(2)

figHandle=figure(42);
set(figHandle, 'Position', [50, 10, 1000, 800]);
bode(Zin_AFE_G1_exp,Zin_AFE_G1_avg,Bode_O)
Bode_Darklines(2)

figHandle=figure(52);
set(figHandle, 'Position', [50, 10, 1000, 800]);
bode(Zin_AFE_G1_avg,Zout_VSI_G1_c1_avg,Zout_VSI_G1_c2_avg,Zout_VSI_G1_c3_avg,Bode_O)
Bode_Darklines(2)

figHandle=figure(62);
set(figHandle, 'Position', [50, 10, 1000, 800]);
plot(Leigenvalues1(1,:),linestyle3,'LineWidth',linewidth)
hold on
plot(Leigenvalues1(2,:),linestyle4,'LineWidth',linewidth)
plot(Leigenvalues2(1,:),linestyle3,'LineWidth',linewidth)
plot(Leigenvalues2(2,:),linestyle4,'LineWidth',linewidth)
plot(Leigenvalues3(1,:),linestyle3,'LineWidth',linewidth)
plot(Leigenvalues3(2,:),linestyle4,'LineWidth',linewidth)

legend({'{\it\lambda}_{1\_case1\_avg}','{\it\lambda}_{2\_case1\_avg}','{\it\lambda}_{1\_case2\_avg}','{\it\lambda}_{2\_case2\_avg}','{\it\lambda}_{1\_case3\_avg}','{\it\lambda}_{2\_case3\_avg}'},'Fontsize',fontsize,'FontWeight','bold')
title('Characteristic Loci of Ldd','Fontsize',fontsize,'FontWeight','bold')
plot(-1,0,'r+','LineWidth',linewidth)
grid
set(gca,'FontSize',fontsize);
ezplot('x^2+y^2=1')

% plot the interface voltage response for three cases
figHandle=figure(72);
set(figHandle, 'Position', [50, 10, 1000, 400]);
plot(Vdq_c1_sim.time,Vdq_c1_sim.signals(1,1).values,'r')
hold on
plot(Vdq_c1_sim.time,Vdq_c1_sim.signals(1,2).values,'k')
xlim([0 0.3])
ylim([-10 120])
xlabel('Time (s)','FontSize',20)
ylabel('Voltage (V)','FontSize',20)
set(findobj(get(gca,'Children'),'LineWidth',0.5),'LineWidth',2);
set(gca,'FontSize',20)
set(gca,'Fontname','Times New Roman')
grid on

figHandle=figure(82);
set(figHandle, 'Position', [50, 10, 1000, 400]);
plot(Vdq_c2_sim.time,Vdq_c2_sim.signals(1,1).values,'r')
hold on
plot(Vdq_c2_sim.time,Vdq_c2_sim.signals(1,2).values,'k')
xlim([0 0.3])
ylim([-10 120])
xlabel('Time (s)','FontSize',20)
ylabel('Voltage (V)','FontSize',20)
set(findobj(get(gca,'Children'),'LineWidth',0.5),'LineWidth',2);
set(gca,'FontSize',20)
set(gca,'Fontname','Times New Roman')
grid on

figHandle=figure(92);
set(figHandle, 'Position', [50, 10, 1000, 400]);
plot(Vdq_c3_sim.time,Vdq_c3_sim.signals(1,1).values,'r')
hold on
plot(Vdq_c3_sim.time,Vdq_c3_sim.signals(1,2).values,'k')
xlim([0 0.3])
ylim([-50 150])
xlabel('Time (s)','FontSize',20)
ylabel('Voltage (V)','FontSize',20)
set(findobj(get(gca,'Children'),'LineWidth',0.5),'LineWidth',2);
set(gca,'FontSize',20)
set(gca,'Fontname','Times New Roman')
grid on