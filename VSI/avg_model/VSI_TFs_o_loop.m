%% AFE parameters and Control to output transfer functions:
clc
clear all
close all

P=bodeoptions;
P.XLabel.FontSize=16;
P.YLabel.FontSize=16;
P.TickLabel.FontSize=16;
P.Title.FontSize=16;
P.title.String=' ';
P.Grid='on';
P.XLim={[20 1e4]};
P.XLimMode={'manual'};
P.FreqUnits='Hz';
P.PhaseWrapping='off';

%% for phase to netural derivation
s=tf([1,0],[0,1]);
omega=400*2*pi;
L=1.1e-3;
Vdc=270;
Vod=99.6;
Voq=0;
Rac=1500;
C=41.3e-6;
Cdc=150e-6;
Rl=0.15;
fsw=20e3;
Tdelay = 1.5/fsw;
% tfdelay = (1-0.5*Tdelay*s+(s*Tdelay)^2/12)/(1+0.5*Tdelay*s+(s*Tdelay)^2/12);
tfdelay = (1-Tdelay*s)/(1+Tdelay*s);
tfidd_hl = (L*C^2*Rac^2*s^3+2*s^2*C*Rac*L+s*L+L*s*C^2*Rac^2*omega^2+C^2*Rac^2*Rl*s^2+2*s*C*Rac*Rl+s*C*Rac^2+C^2*Rac^2*omega^2*Rl+Rac+Rl)*Vdc/(L^2*C^2*Rac^2*s^4+2*L^2*C*Rac*s^3+s^2*L^2+2*L^2*s^2*C^2*Rac^2*omega^2+2*L^2*C*Rac*omega^2*s+omega^2*L^2+L^2*C^2*Rac^2*omega^4+2*L*C^2*Rac^2*Rl*s^3+2*L*s^2*C*Rac^2+4*L*C*Rac*Rl*s^2+2*L*s*C^2*Rac^2*omega^2*Rl+2*L*s*Rac+2*s*L*Rl-2*L*omega^2*C*Rac^2+C^2*Rac^2*Rl^2*s^2+2*s*C*Rac^2*Rl+2*C*Rac*Rl^2*s+Rl^2+Rac^2+2*Rac*Rl+C^2*Rac^2*omega^2*Rl^2);
tfidd_nl_delay =  tfidd_hl*tfdelay;
tfidq_nl = (L*C^2*Rac^2*s^2+2*L*s*C*Rac+L+L*C^2*Rac^2*omega^2-C*Rac^2)*omega*Vdc/(L^2*C^2*Rac^2*s^4+2*L^2*C*Rac*s^3+s^2*L^2+2*L^2*s^2*C^2*Rac^2*omega^2+2*L^2*C*Rac*omega^2*s+omega^2*L^2+L^2*C^2*Rac^2*omega^4+2*L*C^2*Rac^2*Rl*s^3+2*L*s^2*C*Rac^2+4*L*C*Rac*Rl*s^2+2*L*s*C^2*Rac^2*omega^2*Rl+2*L*s*Rac+2*s*L*Rl-2*L*omega^2*C*Rac^2+C^2*Rac^2*Rl^2*s^2+2*s*C*Rac^2*Rl+2*C*Rac*Rl^2*s+Rl^2+Rac^2+2*Rac*Rl+C^2*Rac^2*omega^2*Rl^2);
tfidq_nl_delay = tfidq_nl*tfdelay;
tfiqd_nl = -(L*C^2*Rac^2*s^2+2*L*s*C*Rac+L+L*C^2*Rac^2*omega^2-C*Rac^2)*omega*Vdc/(L^2*C^2*Rac^2*s^4+2*L^2*C*Rac*s^3+s^2*L^2+2*L^2*s^2*C^2*Rac^2*omega^2+2*L^2*C*Rac*omega^2*s+omega^2*L^2+L^2*C^2*Rac^2*omega^4+2*L*C^2*Rac^2*Rl*s^3+2*L*s^2*C*Rac^2+4*L*C*Rac*Rl*s^2+2*L*s*C^2*Rac^2*omega^2*Rl+2*L*s*Rac+2*s*L*Rl-2*L*omega^2*C*Rac^2+C^2*Rac^2*Rl^2*s^2+2*s*C*Rac^2*Rl+2*C*Rac*Rl^2*s+Rl^2+Rac^2+2*Rac*Rl+C^2*Rac^2*omega^2*Rl^2);
tfiqd_nl_delay = tfiqd_nl*tfdelay;
tfiqq_nl = (L*C^2*Rac^2*s^3+2*s^2*C*Rac*L+s*L+L*s*C^2*Rac^2*omega^2+C^2*Rac^2*Rl*s^2+2*s*C*Rac*Rl+s*C*Rac^2+C^2*Rac^2*omega^2*Rl+Rac+Rl)*Vdc/(L^2*C^2*Rac^2*s^4+2*L^2*C*Rac*s^3+s^2*L^2+2*L^2*s^2*C^2*Rac^2*omega^2+2*L^2*C*Rac*omega^2*s+omega^2*L^2+L^2*C^2*Rac^2*omega^4+2*L*C^2*Rac^2*Rl*s^3+2*L*s^2*C*Rac^2+4*L*C*Rac*Rl*s^2+2*L*s*C^2*Rac^2*omega^2*Rl+2*L*s*Rac+2*s*L*Rl-2*L*omega^2*C*Rac^2+C^2*Rac^2*Rl^2*s^2+2*s*C*Rac^2*Rl+2*C*Rac*Rl^2*s+Rl^2+Rac^2+2*Rac*Rl+C^2*Rac^2*omega^2*Rl^2);
tfiqq_nl_delay = tfiqq_nl*tfdelay;

tfvdd_ln = (s^2*C*Rac*L+s*L-L*C*Rac*omega^2+s*C*Rac*Rl+Rl+Rac)*Rac*Vdc/(L^2*C^2*Rac^2*s^4+2*L^2*C*Rac*s^3+s^2*L^2+2*L^2*s^2*C^2*Rac^2*omega^2+2*L^2*C*Rac*omega^2*s+omega^2*L^2+L^2*C^2*Rac^2*omega^4+2*L*C^2*Rac^2*Rl*s^3+2*L*s^2*C*Rac^2+4*L*C*Rac*Rl*s^2+2*L*s*C^2*Rac^2*omega^2*Rl+2*L*s*Rac+2*s*L*Rl-2*L*omega^2*C*Rac^2+C^2*Rac^2*Rl^2*s^2+2*s*C*Rac^2*Rl+2*C*Rac*Rl^2*s+Rl^2+Rac^2+2*Rac*Rl+C^2*Rac^2*omega^2*Rl^2);
tfvdd_ln_delay = tfvdd_ln*tfdelay;
tfvdq_ln = (2*L*s*C*Rac+L+C*Rac*Rl)*omega*Rac*Vdc/(L^2*C^2*Rac^2*s^4+2*L^2*C*Rac*s^3+s^2*L^2+2*L^2*s^2*C^2*Rac^2*omega^2+2*L^2*C*Rac*omega^2*s+omega^2*L^2+L^2*C^2*Rac^2*omega^4+2*L*C^2*Rac^2*Rl*s^3+2*L*s^2*C*Rac^2+4*L*C*Rac*Rl*s^2+2*L*s*C^2*Rac^2*omega^2*Rl+2*L*s*Rac+2*s*L*Rl-2*L*omega^2*C*Rac^2+C^2*Rac^2*Rl^2*s^2+2*s*C*Rac^2*Rl+2*C*Rac*Rl^2*s+Rl^2+Rac^2+2*Rac*Rl+C^2*Rac^2*omega^2*Rl^2);
tfvdq_ln_delay = tfvdq_ln*tfdelay;
tfvqd_ln = -(2*L*s*C*Rac+L+C*Rac*Rl)*omega*Rac*Vdc/(L^2*C^2*Rac^2*s^4+2*L^2*C*Rac*s^3+s^2*L^2+2*L^2*s^2*C^2*Rac^2*omega^2+2*L^2*C*Rac*omega^2*s+omega^2*L^2+L^2*C^2*Rac^2*omega^4+2*L*C^2*Rac^2*Rl*s^3+2*L*s^2*C*Rac^2+4*L*C*Rac*Rl*s^2+2*L*s*C^2*Rac^2*omega^2*Rl+2*L*s*Rac+2*s*L*Rl-2*L*omega^2*C*Rac^2+C^2*Rac^2*Rl^2*s^2+2*s*C*Rac^2*Rl+2*C*Rac*Rl^2*s+Rl^2+Rac^2+2*Rac*Rl+C^2*Rac^2*omega^2*Rl^2);
tfvqd_ln_delay = tfvdq_ln*tfdelay;
tfvqq_ln = (s^2*C*Rac*L+s*L-L*C*Rac*omega^2+s*C*Rac*Rl+Rl+Rac)*Rac*Vdc/(L^2*C^2*Rac^2*s^4+2*L^2*C*Rac*s^3+s^2*L^2+2*L^2*s^2*C^2*Rac^2*omega^2+2*L^2*C*Rac*omega^2*s+omega^2*L^2+L^2*C^2*Rac^2*omega^4+2*L*C^2*Rac^2*Rl*s^3+2*L*s^2*C*Rac^2+4*L*C*Rac*Rl*s^2+2*L*s*C^2*Rac^2*omega^2*Rl+2*L*s*Rac+2*s*L*Rl-2*L*omega^2*C*Rac^2+C^2*Rac^2*Rl^2*s^2+2*s*C*Rac^2*Rl+2*C*Rac*Rl^2*s+Rl^2+Rac^2+2*Rac*Rl+C^2*Rac^2*omega^2*Rl^2);
tfvqq_ln_delay = tfvqq_ln*tfdelay;

% Rac = 15e10;
% Zoutdd_ln = (L^2*C*Rac*s^3+s^2*L^2+L^2*C*Rac*omega^2*s+omega^2*L^2+2*L*C*Rac*Rl*s^2+2*s*L*Rl+L*s*Rac+C*Rac*Rl^2*s+Rac*Rl+Rl^2)*Rac/(L^2*C^2*Rac^2*s^4+2*L^2*C*Rac*s^3+s^2*L^2+2*L^2*s^2*C^2*Rac^2*omega^2+2*L^2*C*Rac*omega^2*s+omega^2*L^2+L^2*C^2*Rac^2*omega^4+2*L*C^2*Rac^2*Rl*s^3+2*L*s^2*C*Rac^2+4*L*C*Rac*Rl*s^2+2*L*s*C^2*Rac^2*omega^2*Rl+2*L*s*Rac+2*s*L*Rl-2*L*omega^2*C*Rac^2+C^2*Rac^2*Rl^2*s^2+2*s*C*Rac^2*Rl+2*C*Rac*Rl^2*s+Rl^2+Rac^2+2*Rac*Rl+C^2*Rac^2*omega^2*Rl^2);
% Zoutdd_ln_delay = Zoutdd_ln*tfdelay;
% Zoutdq_ln = omega*(C*s^2*L^2+2*C*s*L*Rl+C*Rl^2+C*omega^2*L^2-L)*Rac^2/(L^2*C^2*Rac^2*s^4+2*L^2*C*Rac*s^3+s^2*L^2+2*L^2*s^2*C^2*Rac^2*omega^2+2*L^2*C*Rac*omega^2*s+omega^2*L^2+L^2*C^2*Rac^2*omega^4+2*L*C^2*Rac^2*Rl*s^3+2*L*s^2*C*Rac^2+4*L*C*Rac*Rl*s^2+2*L*s*C^2*Rac^2*omega^2*Rl+2*L*s*Rac+2*s*L*Rl-2*L*omega^2*C*Rac^2+C^2*Rac^2*Rl^2*s^2+2*s*C*Rac^2*Rl+2*C*Rac*Rl^2*s+Rl^2+Rac^2+2*Rac*Rl+C^2*Rac^2*omega^2*Rl^2);
% Zoutdq_ln_delay =  Zoutdq_ln*tfdelay;
% Zoutqd_ln = -omega*(C*s^2*L^2+2*C*s*L*Rl+C*Rl^2+C*omega^2*L^2-L)*Rac^2/(L^2*C^2*Rac^2*s^4+2*L^2*C*Rac*s^3+s^2*L^2+2*L^2*s^2*C^2*Rac^2*omega^2+2*L^2*C*Rac*omega^2*s+omega^2*L^2+L^2*C^2*Rac^2*omega^4+2*L*C^2*Rac^2*Rl*s^3+2*L*s^2*C*Rac^2+4*L*C*Rac*Rl*s^2+2*L*s*C^2*Rac^2*omega^2*Rl+2*L*s*Rac+2*s*L*Rl-2*L*omega^2*C*Rac^2+C^2*Rac^2*Rl^2*s^2+2*s*C*Rac^2*Rl+2*C*Rac*Rl^2*s+Rl^2+Rac^2+2*Rac*Rl+C^2*Rac^2*omega^2*Rl^2);
% Zoutqd_ln_delay = Zoutdq_ln*tfdelay;
% Zoutqq_ln = (L^2*C*Rac*s^3+s^2*L^2+L^2*C*Rac*omega^2*s+omega^2*L^2+2*L*C*Rac*Rl*s^2+2*s*L*Rl+L*s*Rac+C*Rac*Rl^2*s+Rac*Rl+Rl^2)*Rac/(L^2*C^2*Rac^2*s^4+2*L^2*C*Rac*s^3+s^2*L^2+2*L^2*s^2*C^2*Rac^2*omega^2+2*L^2*C*Rac*omega^2*s+omega^2*L^2+L^2*C^2*Rac^2*omega^4+2*L*C^2*Rac^2*Rl*s^3+2*L*s^2*C*Rac^2+4*L*C*Rac*Rl*s^2+2*L*s*C^2*Rac^2*omega^2*Rl+2*L*s*Rac+2*s*L*Rl-2*L*omega^2*C*Rac^2+C^2*Rac^2*Rl^2*s^2+2*s*C*Rac^2*Rl+2*C*Rac*Rl^2*s+Rl^2+Rac^2+2*Rac*Rl+C^2*Rac^2*omega^2*Rl^2);
% Zoutqq_ln_delay = Zoutqq_ln*tfdelay;
% Rac = 15;
% Zin_ln = Rac^2*Vdc^2*(L^2*C^2*Rac^2*s^4+2*L^2*C*Rac*s^3+s^2*L^2+2*L^2*s^2*C^2*Rac^2*omega^2+2*L^2*C*Rac*omega^2*s+omega^2*L^2+L^2*C^2*Rac^2*omega^4+2*L*C^2*Rac^2*Rl*s^3+2*L*s^2*C*Rac^2+4*L*C*Rac*Rl*s^2+2*L*s*C^2*Rac^2*omega^2*Rl+2*L*s*Rac+2*s*L*Rl-2*L*omega^2*C*Rac^2+C^2*Rac^2*Rl^2*s^2+2*s*C*Rac^2*Rl+2*C*Rac*Rl^2*s+Rl^2+Rac^2+2*Rac*Rl+C^2*Rac^2*omega^2*Rl^2)/(-omega^4*L^2*C^3*Voq^2*Rac^4*s+omega^4*L^3*C^4*Voq^2*Rac^4*s^3+2*omega^2*L^3*Vod^2*s^2*C*Rac-2*omega^2*L^2*C^3*Voq^2*Rac^4*s^3+3*Rl^2*Voq^2*Rac+omega^2*L^3*Vod^2*C^2*Rac^2*s^3+2*Rl*Voq^2*L*C^2*Rac^3*s^3+2*Rl^2*Voq^2*s^2*C*Rac*L+4*Rl*Voq^2*s^2*C*Rac^2*L+omega^2*L^3*Voq^2*C^2*Rac^2*s^3+2*omega^2*L^3*Voq^2*s^2*C*Rac+2*omega^4*L^3*Voq^2*s*C^2*Rac^2-omega^2*L^2*Voq^2*s*C*Rac^2+2*omega^4*L^2*Voq^2*C^2*Rac^2*Rl+2*Rl*Vod^2*L*C^2*Rac^3*s^3+Rl^2*Vod^2*L*C^2*Rac^2*s^3+4*Rl*Vod^2*s^2*C*Rac^2*L+2*Rl^2*Vod^2*s^2*C*Rac*L+Rl^3*omega^2*C^4*Voq^2*Rac^4*s^2+2*Rl^3*omega^2*C^3*Voq^2*Rac^3*s+Rl^2*omega^2*C^3*Voq^2*Rac^4*s-2*omega^2*L^2*C^3*Vod^2*Rac^4*s^3+omega^4*L^3*C^4*Vod^2*Rac^4*s^3-4*omega^2*L^2*C^2*Vod^2*Rac^3*s^2+2*omega^4*L^3*C^3*Vod^2*Rac^3*s^2-omega^2*L^2*C*Vod^2*Rac^2*s+2*omega^4*L^3*C^2*Vod^2*Rac^2*s-omega^4*L^2*C^3*Vod^2*Rac^4*s+omega^6*L^3*C^4*Vod^2*Rac^4*s-2*omega^4*L*C^3*Vod^2*Rac^4*Rl+omega^6*L^2*C^4*Vod^2*Rac^4*Rl-2*omega^2*L*C*Vod^2*Rac^2*Rl+2*omega^4*L^2*C^2*Vod^2*Rac^2*Rl+2*omega^4*L^3*C^3*Voq^2*Rac^3*s^2-4*omega^2*Rac^3*Voq^2*L^2*C^2*s^2-2*omega^2*L*Voq^2*Rac^2*Rl*C-Rac^4*Vod^2*s*C^2*omega^2*L-2*Rac^4*Voq^2*L*C^3*omega^4*Rl+omega^6*L^2*C^4*Voq^2*Rac^4*Rl+Rl^2*Voq^2*L*C^2*Rac^2*s^3+omega^2*L^2*Voq^2*s^2*C^2*Rac^2*Rl+2*omega^4*L^2*Voq^2*s*C^3*Rac^3*Rl-2*omega^2*L*Voq^2*s*C^2*Rac^3*Rl+Rl*Vod^2*s^2*C^2*Rac^2*L^2*omega^2+2*Rl^2*omega^2*C^3*Voq^2*Rac^3*s^2*L+2*Rl^2*omega^2*C^2*Voq^2*Rac^2*s*L+2*omega^2*L^2*C*Vod^2*Rac*s*Rl+2*omega^4*L^2*C^3*Vod^2*Rac^3*s*Rl+omega^4*L^2*C^4*Voq^2*Rac^4*s^2*Rl-2*omega^2*Rac^4*Voq^2*L*C^3*s^2*Rl+2*omega^2*L^2*Voq^2*s*C*Rac*Rl-2*Rl*Vod^2*L*s*C^2*Rac^3*omega^2+2*Rl^2*Vod^2*L*s*C^2*Rac^2*omega^2+Rl^2*omega^2*C^4*Voq^2*Rac^4*L*s^3+Rl^2*omega^4*C^4*Voq^2*Rac^4*L*s-2*omega^2*L*C^3*Vod^2*Rac^4*Rl*s^2+2*omega^2*L*C^3*Vod^2*Rac^3*Rl^2*s^2+omega^4*L^2*C^4*Vod^2*Rac^4*Rl*s^2+omega^2*L^3*Voq^2*s+omega^2*L^2*Voq^2*Rac+omega^2*L^2*Voq^2*Rl+Rac^2*Vod^2*s*L+Rac^4*Vod^2*s*C+Rl^2*Vod^2*s*L+omega^2*L^2*Vod^2*Rac+omega^2*L^2*Vod^2*Rl+Rl^2*Voq^2*s*L+Rac^2*Voq^2*s*L+Rac^4*Voq^2*s*C+omega^2*L^3*Vod^2*s+3*Rl*Voq^2*Rac^2+3*Rac^2*Vod^2*Rl+3*Rac*Vod^2*Rl^2+Rac^3*Vod^2+Rl^3*Vod^2+Rac^3*Voq^2+2*Rac*Vod^2*s*L*Rl+4*Rac^3*Vod^2*s*C*Rl-2*omega^2*Rac^3*Voq^2*C*L+omega^2*Rac^4*Voq^2*C^2*Rl+Rac^4*Vod^2*L*C^2*s^3+2*Rac^3*Vod^2*s^2*C*L+Rac^4*Vod^2*C^2*Rl*s^2+2*Rac^3*Vod^2*C^2*Rl^2*s^2+5*Rac^2*Vod^2*s*C*Rl^2+Rac^4*Vod^2*C^2*omega^2*Rl+3*Rac^3*Vod^2*C^2*omega^2*Rl^2+Rl^3*Vod^2*C^2*Rac^2*s^2+2*Rl^3*Vod^2*s*C*Rac+2*Rl^3*Vod^2*C^2*Rac^2*omega^2+Rl^3*omega^4*C^4*Voq^2*Rac^4+2*Rl^3*omega^2*C^2*Voq^2*Rac^2+Rac^3*Voq^2*L^2*C^2*omega^4+omega^4*L^2*Vod^2*C^2*Rac^3-2*Rac^3*Vod^2*omega^2*L*C+2*Rl*Voq^2*s*L*Rac+4*Rac^3*Voq^2*s*C*Rl+3*Rl^2*omega^2*C^2*Voq^2*Rac^3+Rl^3*Voq^2*C^2*Rac^2*s^2+2*Rl^2*Voq^2*C^2*Rac^3*s^2+2*Rl^3*Voq^2*s*C*Rac+5*Rl^2*Voq^2*s*C*Rac^2+Rac^4*Voq^2*L*C^2*s^3+2*Rac^3*Voq^2*s^2*C*L+Rac^4*Voq^2*C^2*Rl*s^2+Rl^3*omega^4*C^4*Vod^2*Rac^4+Rl^3*Voq^2+omega^6*L^3*C^4*Voq^2*Rac^4*s+Rl^3*omega^2*C^4*Vod^2*Rac^4*s^2+2*Rl^3*omega^2*C^3*Vod^2*Rac^3*s+Rl^2*omega^2*C^3*Vod^2*Rac^4*s-Rac^4*Voq^2*s*C^2*omega^2*L+Rl^2*omega^2*C^4*Vod^2*Rac^4*L*s^3+Rl^2*omega^4*C^4*Vod^2*Rac^4*L*s);
% Zin_ln_delay = Zin_ln*tfdelay;
% ZinwCdc_ln = 1/(1/Zin_ln+Cdc*s);
% ZinwCdc_ln_delay = 1/(1/Zin_ln_delay+Cdc*s);

%% for line to line derivation
% Vod=99.6*sqrt(3);
% tfidd_ll = (1/3)*(L*C^2*Rac^2*s^3+2*s^2*C*Rac*L+s*L+L*s*C^2*Rac^2*omega^2+C^2*Rac^2*Rl*s^2+s*C*Rac^2+2*s*C*Rac*Rl+Rac+Rl+C^2*Rac^2*omega^2*Rl)*Vdc/(L^2*C^2*Rac^2*s^4+2*L^2*C*Rac*s^3+s^2*L^2+2*L^2*s^2*C^2*Rac^2*omega^2+2*L^2*C*Rac*omega^2*s+L^2*C^2*Rac^2*omega^4+omega^2*L^2+2*L*C^2*Rac^2*Rl*s^3+2*L*s^2*C*Rac^2+4*L*C*Rac*Rl*s^2+2*L*s*C^2*Rac^2*omega^2*Rl+2*s*L*Rl+2*L*s*Rac-2*L*omega^2*C*Rac^2+C^2*Rac^2*Rl^2*s^2+2*C*Rac*Rl^2*s+2*s*C*Rac^2*Rl+2*Rac*Rl+Rl^2+C^2*Rac^2*omega^2*Rl^2+Rac^2);
% 
% tfidq_ll = (1/3)*(L*C^2*Rac^2*s^2+2*L*s*C*Rac+L+L*C^2*Rac^2*omega^2-C*Rac^2)*omega*Vdc/(L^2*C^2*Rac^2*s^4+2*L^2*C*Rac*s^3+s^2*L^2+2*L^2*s^2*C^2*Rac^2*omega^2+2*L^2*C*Rac*omega^2*s+L^2*C^2*Rac^2*omega^4+omega^2*L^2+2*L*C^2*Rac^2*Rl*s^3+2*L*s^2*C*Rac^2+4*L*C*Rac*Rl*s^2+2*L*s*C^2*Rac^2*omega^2*Rl+2*s*L*Rl+2*L*s*Rac-2*L*omega^2*C*Rac^2+C^2*Rac^2*Rl^2*s^2+2*C*Rac*Rl^2*s+2*s*C*Rac^2*Rl+2*Rac*Rl+Rl^2+C^2*Rac^2*omega^2*Rl^2+Rac^2);
% 
% tfiqd_ll = -(1/3)*(L*C^2*Rac^2*s^2+2*L*s*C*Rac+L+L*C^2*Rac^2*omega^2-C*Rac^2)*omega*Vdc/(L^2*C^2*Rac^2*s^4+2*L^2*C*Rac*s^3+s^2*L^2+2*L^2*s^2*C^2*Rac^2*omega^2+2*L^2*C*Rac*omega^2*s+L^2*C^2*Rac^2*omega^4+omega^2*L^2+2*L*C^2*Rac^2*Rl*s^3+2*L*s^2*C*Rac^2+4*L*C*Rac*Rl*s^2+2*L*s*C^2*Rac^2*omega^2*Rl+2*s*L*Rl+2*L*s*Rac-2*L*omega^2*C*Rac^2+C^2*Rac^2*Rl^2*s^2+2*C*Rac*Rl^2*s+2*s*C*Rac^2*Rl+2*Rac*Rl+Rl^2+C^2*Rac^2*omega^2*Rl^2+Rac^2);
% 
% tfiqq_ll = (1/3)*(L*C^2*Rac^2*s^3+2*s^2*C*Rac*L+s*L+L*s*C^2*Rac^2*omega^2+C^2*Rac^2*Rl*s^2+s*C*Rac^2+2*s*C*Rac*Rl+Rac+Rl+C^2*Rac^2*omega^2*Rl)*Vdc/(L^2*C^2*Rac^2*s^4+2*L^2*C*Rac*s^3+s^2*L^2+2*L^2*s^2*C^2*Rac^2*omega^2+2*L^2*C*Rac*omega^2*s+L^2*C^2*Rac^2*omega^4+omega^2*L^2+2*L*C^2*Rac^2*Rl*s^3+2*L*s^2*C*Rac^2+4*L*C*Rac*Rl*s^2+2*L*s*C^2*Rac^2*omega^2*Rl+2*s*L*Rl+2*L*s*Rac-2*L*omega^2*C*Rac^2+C^2*Rac^2*Rl^2*s^2+2*C*Rac*Rl^2*s+2*s*C*Rac^2*Rl+2*Rac*Rl+Rl^2+C^2*Rac^2*omega^2*Rl^2+Rac^2);
% 
% tfvdd_ll = (s^2*C*Rac*L+s*L-L*C*Rac*omega^2+s*C*Rac*Rl+Rl+Rac)*Rac*Vdc/(L^2*C^2*Rac^2*s^4+2*L^2*C*Rac*s^3+s^2*L^2+2*L^2*s^2*C^2*Rac^2*omega^2+2*L^2*C*Rac*omega^2*s+L^2*C^2*Rac^2*omega^4+omega^2*L^2+2*L*C^2*Rac^2*Rl*s^3+2*L*s^2*C*Rac^2+4*L*C*Rac*Rl*s^2+2*L*s*C^2*Rac^2*omega^2*Rl+2*s*L*Rl+2*L*s*Rac-2*L*omega^2*C*Rac^2+C^2*Rac^2*Rl^2*s^2+2*C*Rac*Rl^2*s+2*s*C*Rac^2*Rl+2*Rac*Rl+Rl^2+C^2*Rac^2*omega^2*Rl^2+Rac^2);
% 
% tfvdq_ll = (2*L*s*C*Rac+L+C*Rac*Rl)*omega*Rac*Vdc/(L^2*C^2*Rac^2*s^4+2*L^2*C*Rac*s^3+s^2*L^2+2*L^2*s^2*C^2*Rac^2*omega^2+2*L^2*C*Rac*omega^2*s+L^2*C^2*Rac^2*omega^4+omega^2*L^2+2*L*C^2*Rac^2*Rl*s^3+2*L*s^2*C*Rac^2+4*L*C*Rac*Rl*s^2+2*L*s*C^2*Rac^2*omega^2*Rl+2*s*L*Rl+2*L*s*Rac-2*L*omega^2*C*Rac^2+C^2*Rac^2*Rl^2*s^2+2*C*Rac*Rl^2*s+2*s*C*Rac^2*Rl+2*Rac*Rl+Rl^2+C^2*Rac^2*omega^2*Rl^2+Rac^2);
% 
% tfvqd_ll = -(2*L*s*C*Rac+L+C*Rac*Rl)*omega*Rac*Vdc/(L^2*C^2*Rac^2*s^4+2*L^2*C*Rac*s^3+s^2*L^2+2*L^2*s^2*C^2*Rac^2*omega^2+2*L^2*C*Rac*omega^2*s+L^2*C^2*Rac^2*omega^4+omega^2*L^2+2*L*C^2*Rac^2*Rl*s^3+2*L*s^2*C*Rac^2+4*L*C*Rac*Rl*s^2+2*L*s*C^2*Rac^2*omega^2*Rl+2*s*L*Rl+2*L*s*Rac-2*L*omega^2*C*Rac^2+C^2*Rac^2*Rl^2*s^2+2*C*Rac*Rl^2*s+2*s*C*Rac^2*Rl+2*Rac*Rl+Rl^2+C^2*Rac^2*omega^2*Rl^2+Rac^2);
% 
% tfvqq_ll = (s^2*C*Rac*L+s*L-L*C*Rac*omega^2+s*C*Rac*Rl+Rl+Rac)*Rac*Vdc/(L^2*C^2*Rac^2*s^4+2*L^2*C*Rac*s^3+s^2*L^2+2*L^2*s^2*C^2*Rac^2*omega^2+2*L^2*C*Rac*omega^2*s+L^2*C^2*Rac^2*omega^4+omega^2*L^2+2*L*C^2*Rac^2*Rl*s^3+2*L*s^2*C*Rac^2+4*L*C*Rac*Rl*s^2+2*L*s*C^2*Rac^2*omega^2*Rl+2*s*L*Rl+2*L*s*Rac-2*L*omega^2*C*Rac^2+C^2*Rac^2*Rl^2*s^2+2*C*Rac*Rl^2*s+2*s*C*Rac^2*Rl+2*Rac*Rl+Rl^2+C^2*Rac^2*omega^2*Rl^2+Rac^2);
% 
% Zoutdd_ll = (3*(L^2*C*Rac*s^3+s^2*L^2+L^2*C*Rac*omega^2*s+omega^2*L^2+2*L*C*Rac*Rl*s^2+2*s*L*Rl+L*s*Rac+C*Rac*Rl^2*s+Rac*Rl+Rl^2))*Rac/(L^2*C^2*Rac^2*s^4+2*L^2*C*Rac*s^3+s^2*L^2+2*L^2*s^2*C^2*Rac^2*omega^2+2*L^2*C*Rac*omega^2*s+L^2*C^2*Rac^2*omega^4+omega^2*L^2+2*L*C^2*Rac^2*Rl*s^3+2*L*s^2*C*Rac^2+4*L*C*Rac*Rl*s^2+2*L*s*C^2*Rac^2*omega^2*Rl+2*s*L*Rl+2*L*s*Rac-2*L*omega^2*C*Rac^2+C^2*Rac^2*Rl^2*s^2+2*C*Rac*Rl^2*s+2*s*C*Rac^2*Rl+2*Rac*Rl+Rl^2+C^2*Rac^2*omega^2*Rl^2+Rac^2);
% 
% Zoutdq_ll = 3*omega*(C*s^2*L^2+2*C*s*L*Rl+C*Rl^2+C*omega^2*L^2-L)*Rac^2/(L^2*C^2*Rac^2*s^4+2*L^2*C*Rac*s^3+s^2*L^2+2*L^2*s^2*C^2*Rac^2*omega^2+2*L^2*C*Rac*omega^2*s+L^2*C^2*Rac^2*omega^4+omega^2*L^2+2*L*C^2*Rac^2*Rl*s^3+2*L*s^2*C*Rac^2+4*L*C*Rac*Rl*s^2+2*L*s*C^2*Rac^2*omega^2*Rl+2*s*L*Rl+2*L*s*Rac-2*L*omega^2*C*Rac^2+C^2*Rac^2*Rl^2*s^2+2*C*Rac*Rl^2*s+2*s*C*Rac^2*Rl+2*Rac*Rl+Rl^2+C^2*Rac^2*omega^2*Rl^2+Rac^2);
% 
% Zoutqd_ll = -3*omega*(C*s^2*L^2+2*C*s*L*Rl+C*Rl^2+C*omega^2*L^2-L)*Rac^2/(L^2*C^2*Rac^2*s^4+2*L^2*C*Rac*s^3+s^2*L^2+2*L^2*s^2*C^2*Rac^2*omega^2+2*L^2*C*Rac*omega^2*s+L^2*C^2*Rac^2*omega^4+omega^2*L^2+2*L*C^2*Rac^2*Rl*s^3+2*L*s^2*C*Rac^2+4*L*C*Rac*Rl*s^2+2*L*s*C^2*Rac^2*omega^2*Rl+2*s*L*Rl+2*L*s*Rac-2*L*omega^2*C*Rac^2+C^2*Rac^2*Rl^2*s^2+2*C*Rac*Rl^2*s+2*s*C*Rac^2*Rl+2*Rac*Rl+Rl^2+C^2*Rac^2*omega^2*Rl^2+Rac^2);
% 
% Zoutqq_ll = (3*(L^2*C*Rac*s^3+s^2*L^2+L^2*C*Rac*omega^2*s+omega^2*L^2+2*L*C*Rac*Rl*s^2+2*s*L*Rl+L*s*Rac+C*Rac*Rl^2*s+Rac*Rl+Rl^2))*Rac/(L^2*C^2*Rac^2*s^4+2*L^2*C*Rac*s^3+s^2*L^2+2*L^2*s^2*C^2*Rac^2*omega^2+2*L^2*C*Rac*omega^2*s+L^2*C^2*Rac^2*omega^4+omega^2*L^2+2*L*C^2*Rac^2*Rl*s^3+2*L*s^2*C*Rac^2+4*L*C*Rac*Rl*s^2+2*L*s*C^2*Rac^2*omega^2*Rl+2*s*L*Rl+2*L*s*Rac-2*L*omega^2*C*Rac^2+C^2*Rac^2*Rl^2*s^2+2*C*Rac*Rl^2*s+2*s*C*Rac^2*Rl+2*Rac*Rl+Rl^2+C^2*Rac^2*omega^2*Rl^2+Rac^2);
% 
% Zin_ll = Rac^2*Vdc^2*(L^2*C^2*Rac^2*s^4+2*L^2*C*Rac*s^3+s^2*L^2+2*L^2*s^2*C^2*Rac^2*omega^2+2*L^2*C*Rac*omega^2*s+omega^2*L^2+L^2*C^2*Rac^2*omega^4+2*L*C^2*Rac^2*Rl*s^3+2*L*s^2*C*Rac^2+4*L*C*Rac*Rl*s^2+2*L*s*C^2*Rac^2*omega^2*Rl+2*L*s*Rac+2*s*L*Rl-2*L*omega^2*C*Rac^2+C^2*Rac^2*Rl^2*s^2+2*s*C*Rac^2*Rl+2*C*Rac*Rl^2*s+Rl^2+Rac^2+2*Rac*Rl+C^2*Rac^2*omega^2*Rl^2)/(-omega^4*L^2*C^3*Voq^2*Rac^4*s+omega^4*L^3*C^4*Voq^2*Rac^4*s^3+2*omega^2*L^3*Vod^2*s^2*C*Rac-2*omega^2*L^2*C^3*Voq^2*Rac^4*s^3+3*Rl^2*Voq^2*Rac+omega^2*L^3*Vod^2*C^2*Rac^2*s^3+2*Rl*Voq^2*L*C^2*Rac^3*s^3+2*Rl^2*Voq^2*s^2*C*Rac*L+4*Rl*Voq^2*s^2*C*Rac^2*L+omega^2*L^3*Voq^2*C^2*Rac^2*s^3+2*omega^2*L^3*Voq^2*s^2*C*Rac+2*omega^4*L^3*Voq^2*s*C^2*Rac^2-omega^2*L^2*Voq^2*s*C*Rac^2+2*omega^4*L^2*Voq^2*C^2*Rac^2*Rl+2*Rl*Vod^2*L*C^2*Rac^3*s^3+Rl^2*Vod^2*L*C^2*Rac^2*s^3+4*Rl*Vod^2*s^2*C*Rac^2*L+2*Rl^2*Vod^2*s^2*C*Rac*L+Rl^3*omega^2*C^4*Voq^2*Rac^4*s^2+2*Rl^3*omega^2*C^3*Voq^2*Rac^3*s+Rl^2*omega^2*C^3*Voq^2*Rac^4*s-2*omega^2*L^2*C^3*Vod^2*Rac^4*s^3+omega^4*L^3*C^4*Vod^2*Rac^4*s^3-4*omega^2*L^2*C^2*Vod^2*Rac^3*s^2+2*omega^4*L^3*C^3*Vod^2*Rac^3*s^2-omega^2*L^2*C*Vod^2*Rac^2*s+2*omega^4*L^3*C^2*Vod^2*Rac^2*s-omega^4*L^2*C^3*Vod^2*Rac^4*s+omega^6*L^3*C^4*Vod^2*Rac^4*s-2*omega^4*L*C^3*Vod^2*Rac^4*Rl+omega^6*L^2*C^4*Vod^2*Rac^4*Rl-2*omega^2*L*C*Vod^2*Rac^2*Rl+2*omega^4*L^2*C^2*Vod^2*Rac^2*Rl+2*omega^4*L^3*C^3*Voq^2*Rac^3*s^2-4*omega^2*Rac^3*Voq^2*L^2*C^2*s^2-2*omega^2*L*Voq^2*Rac^2*Rl*C-Rac^4*Vod^2*s*C^2*omega^2*L-2*Rac^4*Voq^2*L*C^3*omega^4*Rl+omega^6*L^2*C^4*Voq^2*Rac^4*Rl+Rl^2*Voq^2*L*C^2*Rac^2*s^3+omega^2*L^2*Voq^2*s^2*C^2*Rac^2*Rl+2*omega^4*L^2*Voq^2*s*C^3*Rac^3*Rl-2*omega^2*L*Voq^2*s*C^2*Rac^3*Rl+Rl*Vod^2*s^2*C^2*Rac^2*L^2*omega^2+2*Rl^2*omega^2*C^3*Voq^2*Rac^3*s^2*L+2*Rl^2*omega^2*C^2*Voq^2*Rac^2*s*L+2*omega^2*L^2*C*Vod^2*Rac*s*Rl+2*omega^4*L^2*C^3*Vod^2*Rac^3*s*Rl+omega^4*L^2*C^4*Voq^2*Rac^4*s^2*Rl-2*omega^2*Rac^4*Voq^2*L*C^3*s^2*Rl+2*omega^2*L^2*Voq^2*s*C*Rac*Rl-2*Rl*Vod^2*L*s*C^2*Rac^3*omega^2+2*Rl^2*Vod^2*L*s*C^2*Rac^2*omega^2+Rl^2*omega^2*C^4*Voq^2*Rac^4*L*s^3+Rl^2*omega^4*C^4*Voq^2*Rac^4*L*s-2*omega^2*L*C^3*Vod^2*Rac^4*Rl*s^2+2*omega^2*L*C^3*Vod^2*Rac^3*Rl^2*s^2+omega^4*L^2*C^4*Vod^2*Rac^4*Rl*s^2+omega^2*L^3*Voq^2*s+omega^2*L^2*Voq^2*Rac+omega^2*L^2*Voq^2*Rl+Rac^2*Vod^2*s*L+Rac^4*Vod^2*s*C+Rl^2*Vod^2*s*L+omega^2*L^2*Vod^2*Rac+omega^2*L^2*Vod^2*Rl+Rl^2*Voq^2*s*L+Rac^2*Voq^2*s*L+Rac^4*Voq^2*s*C+omega^2*L^3*Vod^2*s+3*Rl*Voq^2*Rac^2+3*Rac^2*Vod^2*Rl+3*Rac*Vod^2*Rl^2+Rac^3*Vod^2+Rl^3*Vod^2+Rac^3*Voq^2+2*Rac*Vod^2*s*L*Rl+4*Rac^3*Vod^2*s*C*Rl-2*omega^2*Rac^3*Voq^2*C*L+omega^2*Rac^4*Voq^2*C^2*Rl+Rac^4*Vod^2*L*C^2*s^3+2*Rac^3*Vod^2*s^2*C*L+Rac^4*Vod^2*C^2*Rl*s^2+2*Rac^3*Vod^2*C^2*Rl^2*s^2+5*Rac^2*Vod^2*s*C*Rl^2+Rac^4*Vod^2*C^2*omega^2*Rl+3*Rac^3*Vod^2*C^2*omega^2*Rl^2+Rl^3*Vod^2*C^2*Rac^2*s^2+2*Rl^3*Vod^2*s*C*Rac+2*Rl^3*Vod^2*C^2*Rac^2*omega^2+Rl^3*omega^4*C^4*Voq^2*Rac^4+2*Rl^3*omega^2*C^2*Voq^2*Rac^2+Rac^3*Voq^2*L^2*C^2*omega^4+omega^4*L^2*Vod^2*C^2*Rac^3-2*Rac^3*Vod^2*omega^2*L*C+2*Rl*Voq^2*s*L*Rac+4*Rac^3*Voq^2*s*C*Rl+3*Rl^2*omega^2*C^2*Voq^2*Rac^3+Rl^3*Voq^2*C^2*Rac^2*s^2+2*Rl^2*Voq^2*C^2*Rac^3*s^2+2*Rl^3*Voq^2*s*C*Rac+5*Rl^2*Voq^2*s*C*Rac^2+Rac^4*Voq^2*L*C^2*s^3+2*Rac^3*Voq^2*s^2*C*L+Rac^4*Voq^2*C^2*Rl*s^2+Rl^3*omega^4*C^4*Vod^2*Rac^4+Rl^3*Voq^2+omega^6*L^3*C^4*Voq^2*Rac^4*s+Rl^3*omega^2*C^4*Vod^2*Rac^4*s^2+2*Rl^3*omega^2*C^3*Vod^2*Rac^3*s+Rl^2*omega^2*C^3*Vod^2*Rac^4*s-Rac^4*Voq^2*s*C^2*omega^2*L+Rl^2*omega^2*C^4*Vod^2*Rac^4*L*s^3+Rl^2*omega^4*C^4*Vod^2*Rac^4*L*s);
% 
% ZinwCdc_ll = 1/(1/Zin_ll+Cdc*s);
% 
% load TFidqdq_VSI.mat
% % load TFidqdq_VSI_delay.mat
% tfidd_matlab = tfidqdq_VSI(1,1);
% tfidq_matlab = tfidqdq_VSI(1,2);
% tfiqd_matlab = tfidqdq_VSI(2,1);
% tfiqq_matlab = tfidqdq_VSI(2,2);
% 
% load TFvdqdq_VSI.mat
% % load TFvdqdq_VSI_delay.mat
% tfvddd_matlab = tfvddq_VSI(1,1);
% tfvddq_matlab = tfvddq_VSI(1,2);
% tfvqdd_matlab = tfvddq_VSI(2,1);
% tfvqdq_matlab = tfvddq_VSI(2,2);
% 
% load ZoutVSI_o_loop.mat
% % load ZoutVSI_o_loop_delay.mat
% Zoutdd_matlab = ZoutVSI(1,1);
% Zoutdq_matlab = ZoutVSI(1,2);
% Zoutqd_matlab = ZoutVSI(2,1);
% Zoutqq_matlab = ZoutVSI(2,2);
% 
% load ZinVSI_o_loop.mat
% % load ZinVSI_o_loop_delay.mat
% Zin_matlab = ZinVSI;
% ZinwCdc_matlab = 1/(1/Zin_matlab+Cdc*s);

% figure(1);
% bode(tfidd_ln,P);
% hold on
% % bode(tfidd_ln_delay,P);
% bode(tfidd_matlab,P);
% Bode_Darklines(1);
% hold off
% 
% figure(2);
% bode(tfidq_ln,P);
% hold on
% % bode(tfidq_ln_delay,P);
% bode(tfidq_matlab,P);
% Bode_Darklines(1);
% hold off
% 
% figure(3);
% bode(tfiqd_ln,P);
% hold on
% % bode(tfiqd_ln_delay,P);
% bode(tfiqd_matlab,P);
% Bode_Darklines(1);
% hold off
% 
% figure(4);
% bode(tfiqq_ln,P);
% hold on
% % bode(tfiqq_ln_delay,P);
% bode(tfiqq_matlab,P);
% Bode_Darklines(1);
% hold off
% 
% figure(5);
% bode(tfvdd_ln,P);
% hold on
% % bode(tfvdd_ln_delay,P);
% bode(tfvddd_matlab,P);
% Bode_Darklines(1);
% hold off
% 
% figure(6);
% bode(tfvdq_ln,P);
% hold on
% % bode(tfvdq_ln_delay,P);
% bode(tfvddq_matlab,P);
% Bode_Darklines(1);
% hold off
% 
% figure(7);
% bode(tfvqd_ln,P);
% hold on
% % bode(tfvqd_ln_delay,P);
% bode(tfvqdd_matlab,P);
% Bode_Darklines(1);
% hold off
% 
% figure(8);
% bode(tfvqq_ln,P);
% hold on
% % bode(tfvqq_ln_delay,P);
% bode(tfvqdq_matlab,P);
% Bode_Darklines(1);
% hold off
% 
% figure(9);
% bode(Zoutdd_ln,P);
% hold on
% % bode(Zoutdd_ln_delay,P);
% bode(Zoutdd_matlab,P);
% Bode_Darklines(1);
% hold off
% 
% figure(10);
% bode(Zoutdq_ln,P);
% hold on
% % bode(Zoutdq_ln_delay,P);
% bode(Zoutdq_matlab,P);
% Bode_Darklines(1);
% hold off
% 
% figure(11);
% bode(Zoutqd_ln,P);
% hold on
% % bode(Zoutqd_ln_delay,P);
% bode(Zoutqd_matlab,P);
% Bode_Darklines(1);
% hold off
% 
% figure(12);
% bode(Zoutqq_ln,P);
% hold on
% % bode(Zoutqq_ln_delay,P);
% bode(Zoutqq_matlab,P);
% Bode_Darklines(1);
% hold off
% 
% figure(13);
% bode(Zin_ln,P);
% hold on
% % bode(Zin_ln_delay,P);
% bode(Zin_matlab,P);
% Bode_Darklines(1);
% hold off
% 
% figure(14);
% bode(ZinwCdc_ln,P);
% hold on
% % bode(ZinwCdc_ln_delay,P);
% bode(ZinwCdc_matlab,P);
% Bode_Darklines(1);
% hold off
% 
% fprintf('Old current controller design:\n')
% fi=100;
% Kpi=2*pi*fi*L/Vdc 
% Kii=2*pi*fi*(Rl)/Vdc

%% design the controller when load is inf, and then check the phase
% margin and bandwidth at the norminal load
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Design of conventional PI controllers
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Current loop: transfer function to be compensated:1/(sLp+RLp)

kpi=2*pi*800*L/Vdc*0.5
kii=2*pi*800*(Rl)/Vdc*10
tf_i_c=tf([kpi  kii],[1  0]);

Rac=1.5e3;
tfidd_ll = (L*C^2*Rac^2*s^3+2*s^2*C*Rac*L+s*L+L*s*C^2*Rac^2*omega^2+C^2*Rac^2*Rl*s^2+2*s*C*Rac*Rl+s*C*Rac^2+C^2*Rac^2*omega^2*Rl+Rac+Rl)*Vdc/(L^2*C^2*Rac^2*s^4+2*L^2*C*Rac*s^3+s^2*L^2+2*L^2*s^2*C^2*Rac^2*omega^2+2*L^2*C*Rac*omega^2*s+omega^2*L^2+L^2*C^2*Rac^2*omega^4+2*L*C^2*Rac^2*Rl*s^3+2*L*s^2*C*Rac^2+4*L*C*Rac*Rl*s^2+2*L*s*C^2*Rac^2*omega^2*Rl+2*L*s*Rac+2*s*L*Rl-2*L*omega^2*C*Rac^2+C^2*Rac^2*Rl^2*s^2+2*s*C*Rac^2*Rl+2*C*Rac*Rl^2*s+Rl^2+Rac^2+2*Rac*Rl+C^2*Rac^2*omega^2*Rl^2);
tfidd_ll_del=tfidd_ll*tfdelay;		  				% 1/(sLp+RLp)
tfvdd_ll = (s^2*C*Rac*L+s*L-L*C*Rac*omega^2+s*C*Rac*Rl+Rl+Rac)*Rac*Vdc/(L^2*C^2*Rac^2*s^4+2*L^2*C*Rac*s^3+s^2*L^2+2*L^2*s^2*C^2*Rac^2*omega^2+2*L^2*C*Rac*omega^2*s+omega^2*L^2+L^2*C^2*Rac^2*omega^4+2*L*C^2*Rac^2*Rl*s^3+2*L*s^2*C*Rac^2+4*L*C*Rac*Rl*s^2+2*L*s*C^2*Rac^2*omega^2*Rl+2*L*s*Rac+2*s*L*Rl-2*L*omega^2*C*Rac^2+C^2*Rac^2*Rl^2*s^2+2*s*C*Rac^2*Rl+2*C*Rac*Rl^2*s+Rl^2+Rac^2+2*Rac*Rl+C^2*Rac^2*omega^2*Rl^2);
tfvdd_ll_del=tfvdd_ll*tfdelay;
Rac=15;
tfidd_nl = (L*C^2*Rac^2*s^3+2*s^2*C*Rac*L+s*L+L*s*C^2*Rac^2*omega^2+C^2*Rac^2*Rl*s^2+2*s*C*Rac*Rl+s*C*Rac^2+C^2*Rac^2*omega^2*Rl+Rac+Rl)*Vdc/(L^2*C^2*Rac^2*s^4+2*L^2*C*Rac*s^3+s^2*L^2+2*L^2*s^2*C^2*Rac^2*omega^2+2*L^2*C*Rac*omega^2*s+omega^2*L^2+L^2*C^2*Rac^2*omega^4+2*L*C^2*Rac^2*Rl*s^3+2*L*s^2*C*Rac^2+4*L*C*Rac*Rl*s^2+2*L*s*C^2*Rac^2*omega^2*Rl+2*L*s*Rac+2*s*L*Rl-2*L*omega^2*C*Rac^2+C^2*Rac^2*Rl^2*s^2+2*s*C*Rac^2*Rl+2*C*Rac*Rl^2*s+Rl^2+Rac^2+2*Rac*Rl+C^2*Rac^2*omega^2*Rl^2);
tfidd_nl_del=tfidd_nl*tfdelay;
tfvdd_nl = (s^2*C*Rac*L+s*L-L*C*Rac*omega^2+s*C*Rac*Rl+Rl+Rac)*Rac*Vdc/(L^2*C^2*Rac^2*s^4+2*L^2*C*Rac*s^3+s^2*L^2+2*L^2*s^2*C^2*Rac^2*omega^2+2*L^2*C*Rac*omega^2*s+omega^2*L^2+L^2*C^2*Rac^2*omega^4+2*L*C^2*Rac^2*Rl*s^3+2*L*s^2*C*Rac^2+4*L*C*Rac*Rl*s^2+2*L*s*C^2*Rac^2*omega^2*Rl+2*L*s*Rac+2*s*L*Rl-2*L*omega^2*C*Rac^2+C^2*Rac^2*Rl^2*s^2+2*s*C*Rac^2*Rl+2*C*Rac*Rl^2*s+Rl^2+Rac^2+2*Rac*Rl+C^2*Rac^2*omega^2*Rl^2);
tfvdd_nl_del=tfvdd_nl*tfdelay;

Rac=5;
tfidd_hl = (L*C^2*Rac^2*s^3+2*s^2*C*Rac*L+s*L+L*s*C^2*Rac^2*omega^2+C^2*Rac^2*Rl*s^2+2*s*C*Rac*Rl+s*C*Rac^2+C^2*Rac^2*omega^2*Rl+Rac+Rl)*Vdc/(L^2*C^2*Rac^2*s^4+2*L^2*C*Rac*s^3+s^2*L^2+2*L^2*s^2*C^2*Rac^2*omega^2+2*L^2*C*Rac*omega^2*s+omega^2*L^2+L^2*C^2*Rac^2*omega^4+2*L*C^2*Rac^2*Rl*s^3+2*L*s^2*C*Rac^2+4*L*C*Rac*Rl*s^2+2*L*s*C^2*Rac^2*omega^2*Rl+2*L*s*Rac+2*s*L*Rl-2*L*omega^2*C*Rac^2+C^2*Rac^2*Rl^2*s^2+2*s*C*Rac^2*Rl+2*C*Rac*Rl^2*s+Rl^2+Rac^2+2*Rac*Rl+C^2*Rac^2*omega^2*Rl^2);
tfidd_hl_del=tfidd_hl*tfdelay;
tfvdd_hl = (s^2*C*Rac*L+s*L-L*C*Rac*omega^2+s*C*Rac*Rl+Rl+Rac)*Rac*Vdc/(L^2*C^2*Rac^2*s^4+2*L^2*C*Rac*s^3+s^2*L^2+2*L^2*s^2*C^2*Rac^2*omega^2+2*L^2*C*Rac*omega^2*s+omega^2*L^2+L^2*C^2*Rac^2*omega^4+2*L*C^2*Rac^2*Rl*s^3+2*L*s^2*C*Rac^2+4*L*C*Rac*Rl*s^2+2*L*s*C^2*Rac^2*omega^2*Rl+2*L*s*Rac+2*s*L*Rl-2*L*omega^2*C*Rac^2+C^2*Rac^2*Rl^2*s^2+2*s*C*Rac^2*Rl+2*C*Rac*Rl^2*s+Rl^2+Rac^2+2*Rac*Rl+C^2*Rac^2*omega^2*Rl^2);
tfvdd_hl_del=tfvdd_hl*tfdelay;

figure(1)
bode(tfidd_ll_del,tfidd_nl_del,tfidd_hl_del,P);
Bode_Darklines(3);
grid on


% %% light load loop gain
% T_i_ll=tfidd_ll_del*tf_i_c;					% open loop gain
% WW_ll=feedback(T_i_ll,1,-1);						% transfer function bw i_ref and i
% figure(2)
% step(WW_ll)										% step response
% figure(3)
% bode(T_i_ll,P)
% grid on
% Bode_Darklines(3);
% title('Current Control Loop Gain')
%% normal load loop gain
T_i_nl=tfidd_nl_del*tf_i_c;
tf_i_cl=feedback(T_i_nl,1,-1);						% transfer function bw i_ref and i
figure(4)
step(tf_i_cl)										% step response
figure(5)
bode(T_i_nl,P)
grid on
Bode_Darklines(3);
title('Current Control Loop Gain')
figure(6)
bode(tf_i_cl,P)
grid on
Bode_Darklines(3);
title('Current Control Loop Gain')
% %% high load loop gain
% T_i_hl=tfidd_hl_del*tf_i_c;
% WW_hl=feedback(T_i_hl,1,-1);						% transfer function bw i_ref and i
% figure(7)
% step(WW_hl)										% step response
% figure(8)
% bode(T_i_hl,P)
% grid on
% Bode_Darklines(3)
% title('Current Control Loop Gain')

