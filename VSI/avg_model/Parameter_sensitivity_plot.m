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
Bode_O.XLim={[1 1e4]};
Bode_O.XLimMode={'manual'};
Bode_O.FreqUnits='Hz';
Bode_O.PhaseWrapping='off';

load Zo_vsi_L_1_C_32.mat
Zo_VSI_L_1_C_32=Zo_vil_cal_frd;
Zo_VSI_L_1_C_32_dd=Zo_vil_cal_frd(1,1);
load Zo_vsi_L_1_C_0_95_32.mat
Zo_VSI_L_1_C_30=Zo_vil_cal_frd;
Zo_VSI_L_1_C_30_dd=Zo_vil_cal_frd(1,1);
load Zo_vsi_L_1_C_1_05_32.mat
Zo_VSI_L_1_C_34=Zo_vil_cal_frd;
Zo_VSI_L_1_C_34_dd=Zo_vil_cal_frd(1,1);

figure(1)
bode(Zo_VSI_L_1_C_30_dd,Zo_VSI_L_1_C_32_dd,Zo_VSI_L_1_C_34_dd,Bode_O)

load Zo_vsi_L_0_9_C_32.mat
Zo_VSI_L_0_9_C_32=Zo_vil_cal_frd;
Zo_VSI_L_0_9_C_32_dd=Zo_vil_cal_frd(1,1);
load Zo_vsi_L_1_1_C_32.mat
Zo_VSI_L_1_1_C_32=Zo_vil_cal_frd;
Zo_VSI_L_1_1_C_32_dd=Zo_vil_cal_frd(1,1);

figure(2)
bode(Zo_VSI_L_0_9_C_32_dd,Zo_VSI_L_1_C_32_dd,Zo_VSI_L_1_1_C_32_dd,Bode_O)