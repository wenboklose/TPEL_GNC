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

load Z_cal_control_1.mat
Zin_vil_pll_100_v_200_i_2k=Zin_vil_pll_cal_frd;
Yin_vil_pll_100_v_200_i_2k=1/Zin_vil_pll_cal_frd;
load Z_cal_control_2.mat
Zin_vil_pll_100_v_200_i_1_6k=Zin_vil_pll_cal_frd;
Yin_vil_pll_100_v_200_i_1_6k=1/Zin_vil_pll_cal_frd;
load Z_cal_control_3.mat
Zin_vil_pll_100_v_200_i_1k=Zin_vil_pll_cal_frd;
Yin_vil_pll_100_v_200_i_1k=1/Zin_vil_pll_cal_frd;
load Z_cal_control_4.mat
Zin_vil_pll_100_v_100_i_2k=Zin_vil_pll_cal_frd;
Yin_vil_pll_100_v_100_i_2k=1/Zin_vil_pll_cal_frd;
load Z_cal_control_5.mat
Zin_vil_pll_100_v_50_i_2k=Zin_vil_pll_cal_frd;
Yin_vil_pll_100_v_50_i_2k=1/Zin_vil_pll_cal_frd;
load Z_cal_control_6.mat
Zin_vil_pll_300_v_200_i_2k=Zin_vil_pll_cal_frd;
Yin_vil_pll_300_v_200_i_2k=1/Zin_vil_pll_cal_frd;
load Z_cal_control_7.mat
Zin_vil_pll_50_v_200_i_2k=Zin_vil_pll_cal_frd;
Yin_vil_pll_50_v_200_i_2k=1/Zin_vil_pll_cal_frd;

figure(1)
bode(Yin_vil_pll_100_v_200_i_1k,Yin_vil_pll_100_v_200_i_1_6k,Yin_vil_pll_100_v_200_i_2k,Bode_O)
Bode_Darklines(3)
figure(2)
bode(Zin_vil_pll_100_v_200_i_1k,Zin_vil_pll_100_v_200_i_1_6k,Zin_vil_pll_100_v_200_i_2k,Bode_O)
Bode_Darklines(3)
figure(3)
bode(Zin_vil_pll_100_v_200_i_2k,Zin_vil_pll_100_v_100_i_2k,Zin_vil_pll_100_v_50_i_2k,Bode_O)
Bode_Darklines(3)
figure(4)
bode(Zin_vil_pll_50_v_200_i_2k,Zin_vil_pll_100_v_200_i_2k,Zin_vil_pll_300_v_200_i_2k,Bode_O)
Bode_Darklines(3)
figure(5)
bode(Yin_vil_pll_50_v_200_i_2k,Yin_vil_pll_100_v_200_i_2k,Yin_vil_pll_300_v_200_i_2k,Bode_O)
Bode_Darklines(3)