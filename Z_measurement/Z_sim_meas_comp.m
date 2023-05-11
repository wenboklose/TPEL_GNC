%% compare simulation results and expreiments results

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
linewidth=3; fontsize=18;
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

load ('Z_afe_100p_vsi_fi_4k_fv200_afe_fi_1k_fv_70.mat')
Z_AFE_400_fi_1k_fv_70_pll_1_test = Z.ZDQRAW/3;


cd ..
cd ..
cd VSI/07_01_12_STASU_deliverable/Models/sw_model_SC
load ('Z_VSI_fi_4k_fv_80_R_12_50p.mat')
Z_VSI_o_400_fi_4k_fv_80_sim = Z.ZDQSRAW/3;
load ('Z_VSI_fi_4k_fv_150_R_12_50p.mat')
Z_VSI_o_400_fi_4k_fv_150_sim = Z.ZDQSRAW/3;
load ('Z_VSI_fi_4k_fv_200_R_12_50p.mat')
Z_VSI_o_400_fi_4k_fv_200_sim = Z.ZDQSRAW/3;
load ('Z_VSI_fi_4k_fv_250_R_12_50p.mat')
Z_VSI_o_400_fi_4k_fv_250_sim = Z.ZDQSRAW/3;
cd ..
cd ..
cd ..
cd ..
cd VSI_AFE/Z_measurement

cd ..
cd ..
cd VSI/avg_model
load ('Z_VSI_fi_4k_fv_80_R_12_avg.mat')
Z_VSI_o_400_fi_4k_fv_80_avg = Zo_vil_cal_frd;
load ('Z_VSI_fi_4k_fv_150_R_12_avg.mat')
Z_VSI_o_400_fi_4k_fv_150_avg = Zo_vil_cal_frd;
load ('Z_VSI_fi_4k_fv_200_R_12_avg.mat')
Z_VSI_o_400_fi_4k_fv_200_avg = Zo_vil_cal_frd;
load ('Z_VSI_fi_4k_fv_250_R_12_avg.mat')
Z_VSI_o_400_fi_4k_fv_250_avg = Zo_vil_cal_frd;
cd ..
cd ..
cd VSI_AFE/Z_measurement

cd ..
cd ..
cd AFE/07_01_12_STASU_deliverable/Models/sw_model_SV
load ('Z_AFE_fi_4k_fv_70_50p.mat')
Z_AFE_400_fi_1k_fv_70_pll_1_sim = Z.ZDQLRAW/3;
cd ..
cd ..
cd ..
cd avg
load('Z_AFE_fi_1k_fv_70_pll_1_avg.mat');
Z_AFE_400_fi_1k_fv_70_pll_1_avg = Zin_vl_pll_avg_sim;
cd ..
cd ..
cd VSI_AFE/Z_measurement

fighandle=figure(1)
set(fighandle,'position',[10, 10, 1000, 800])
bode(Z_VSI_o_400_fi_4k_fv_80_test,Bode_O);
hold on
bode(Z_VSI_o_400_fi_4k_fv_80_sim,Bode_O);
bode(Z_VSI_o_400_fi_4k_fv_80_avg,Bode_O)
Bode_Darklines(3);

fighandle=figure(2)
set(fighandle,'position',[10, 10, 1000, 800])
bode(Z_VSI_o_400_fi_4k_fv_150_test,Bode_O);
hold on
bode(Z_VSI_o_400_fi_4k_fv_150_sim,Bode_O);
bode(Z_VSI_o_400_fi_4k_fv_150_avg,Bode_O)
Bode_Darklines(3);

fighandle=figure(3)
set(fighandle,'position',[10, 10, 1000, 800])
bode(Z_VSI_o_400_fi_4k_fv_250_test,Bode_O);
hold on
bode(Z_VSI_o_400_fi_4k_fv_250_sim,Bode_O);
bode(Z_VSI_o_400_fi_4k_fv_250_avg,Bode_O)
Bode_Darklines(3);

fighandle=figure(4)
set(fighandle,'position',[10, 10, 1000, 800])
bode(Z_AFE_400_fi_1k_fv_70_pll_1_test,Bode_O);
hold on
bode(Z_AFE_400_fi_1k_fv_70_pll_1_sim,Bode_O);
bode(Z_AFE_400_fi_1k_fv_70_pll_1_avg,Bode_O)
Bode_Darklines(3);