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
Bode_O.XLim={[10 1e4]};
Bode_O.XLimMode={'manual'};
Bode_O.FreqUnits='Hz';
Bode_O.PhaseWrapping='on';

load Z_AFE_ab_PI.mat
Z_AFE_ab_PI_cal=3*Zin_vil_pll_cal_frd;
Y_AFE_ab_PI_cal=1/Zin_vil_pll_cal_frd/3;

load Z_AFE_ab_PR.mat
Z_AFE_ab_PR_cal=3*Zin_vil_pll_cal_frd;
Y_AFE_ab_PR_cal=1/Zin_vil_pll_cal_frd/3;

load Z_AFE_ab_PI_v1.mat
Z_AFE_ab_PI_cal_v1=3*Zin_vil_pll_cal_frd;
Y_AFE_ab_PI_cal_v1=1/Zin_vil_pll_cal_frd/3;

load Z_AFE_ab_PR_v1.mat
Z_AFE_ab_PR_cal_v1=3*Zin_vil_pll_cal_frd;
Y_AFE_ab_PR_cal_v1=1/Zin_vil_pll_cal_frd/3;
figure(1)
bode(Z_AFE_ab_PI_cal,Z_AFE_ab_PI_cal_v1,Bode_O)
figure(2)
bode(Z_AFE_ab_PR_cal,Z_AFE_ab_PR_cal_v1,Bode_O)