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
Bode_O.XLim={[1 1e8]};
Bode_O.XLimMode={'manual'};
Bode_O.FreqUnits='Hz';
Bode_O.PhaseWrapping='off';

f=400;
omega=2*pi*f;

nsamples=100;
K_norm=1;C1_norm=75e-12;C2_norm=39e-12;R1_norm=15e3;R2_norm=15e3;
Nom=[K_norm C1_norm C2_norm R1_norm R2_norm];
Nk=nsamples;
Nc=size(Nom,2);
randn('state',sum(100*clock));
%
Tn=zeros(Nc,Nk); % Reserve space; decreases run time
%
% Create tolerance array T;
% 
Tr=0.1;Te=0.05;Tf=0.001;
T=[-Tr -Te -Te -Tf -Tf;Tr Te Te Tf Tf];
% T below for fig 20 and fig 21
%T=[ -Tr -Tr -Tr -Tr -Te -Te;0.03 Tr 0.03 Tr Te Te];
% Convert to tolerance multipliers
for k=1:Nk
   for w=1:Nc
      Tn(w,k)=Nom(w)*(((T(2,w)-T(1,w))/6)*(randn+3)+T(1,w)+1);
   end % end w loop
end % end k loop
K_rec=Tn(1,:);
C1_rec=Tn(2,:);
C2_rec=Tn(3,:);
R1_rec=Tn(4,:);
R2_rec=Tn(5,:);
for i=1:1:nsamples
    
%         C1=75e-12;C2=39e-12;R1=15e3;R2=15e3;
        C1=C1_rec(i);C2=C2_rec(i);R1=R1_rec(i);R2=R2_rec(i);
        K=K_rec(i);
        tf_filter = K*tf([0 1],[C1*C2*R1*R2 C2*(R1+R2) 1]);
        tf_filter_dq_s = tf(JF_DQFromABC(tf_filter,omega));
        tf_filter_dq{i} = [tf_filter_dq_s(1,1) 1*tf_filter_dq_s(1,2); -1*tf_filter_dq_s(1,2) tf_filter_dq_s(2,2)];
    
end

% calculation and simulation comparison
figure(1)
hold on
for i=1:1:nsamples*1
    bode(tf_filter_dq{i}(1,1),Bode_O)
end