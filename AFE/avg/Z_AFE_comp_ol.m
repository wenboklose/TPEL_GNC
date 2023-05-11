%% simulatio results comparison with calculation
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

load('ZAFE_in_pll_1.mat')
load('ZAFE_in_pll_10.mat')
load('ZAFE_in_pll_100.mat')
load('ZAFE_in_pll_400.mat')
load('ZAFE_in_pll_800.mat')

figure(1)
bode(Zin_pll_1_avg_sim,'b',Bode_O)
hold on
bode(Zin_pll_1_cal,'-or',Bode_O)
bode(Zin_pll_10_avg_sim,'b',Bode_O)
bode(Zin_pll_10_cal,'-or',Bode_O)
bode(Zin_pll_100_avg_sim,'b',Bode_O)
bode(Zin_pll_100_cal,'-or',Bode_O)
bode(Zin_pll_400_avg_sim,'b',Bode_O)
bode(Zin_pll_400_cal,'-or',Bode_O)
bode(Zin_pll_800_avg_sim,'b',Bode_O)
bode(Zin_pll_800_cal,'-or',Bode_O)
Bode_Darklines(3)

%% load switching model simulation
load('Z_AFE_pll_100_sw.mat')
fit=fit*2*pi;
ZAFE_in_dd=frd(Zldd,fit);
ZAFE_in_dq=frd(Zldq,fit);
ZAFE_in_qd=frd(Zlqd,fit);
ZAFE_in_qq=frd(Zlqq,fit);
Zin_pll_100_sw_sim=[ZAFE_in_dd ZAFE_in_dq; ZAFE_in_qd ZAFE_in_qq];
load('Z_AFE_pll_400_sw.mat')
fit=fit*2*pi;
ZAFE_in_dd=frd(Zldd,fit);
ZAFE_in_dq=frd(Zldq,fit);
ZAFE_in_qd=frd(Zlqd,fit);
ZAFE_in_qq=frd(Zlqq,fit);
Zin_pll_400_sw_sim=[ZAFE_in_dd ZAFE_in_dq; ZAFE_in_qd ZAFE_in_qq];
figure(6)
bode(Zin_pll_100_sw_sim,'b',Bode_O)
hold on
bode(Zin_pll_100_cal,'-or',Bode_O)
Bode_Darklines(3)
figure(7)
bode(Zin_pll_400_sw_sim,'b',Bode_O)
hold on
bode(Zin_pll_400_cal,'-or',Bode_O)
Bode_Darklines(3)
%% load switching model simulation series injection method
load('Z_in_AFE_pll_100_sj.mat')
load('Z_in_AFE_pll_400_sj.mat')
load('Z_in_AFE_pll_800_sj.mat')
Zin_pll_100_cal_sw_sim_sj=Z_in_AFE_pll_100_sj;
Zin_pll_400_cal_sw_sim_sj=Z_in_AFE_pll_400_sj;
Zin_pll_800_cal_sw_sim_sj=Z_in_AFE_pll_800_sj;

figure(8)
bode(Zin_pll_100_cal_sw_sim_sj,'b',Bode_O)
hold on
bode(Zin_pll_100_cal,'-or',Bode_O)
Bode_Darklines(3)

figure(9)
bode(Zin_pll_400_cal_sw_sim_sj,'b',Bode_O)
hold on
bode(Zin_pll_400_cal,'-or',Bode_O)
Bode_Darklines(3)

figure(10)
bode(Zin_pll_800_cal_sw_sim_sj,'b',Bode_O)
hold on
bode(Zin_pll_800_cal,'-or',Bode_O)
Bode_Darklines(3)