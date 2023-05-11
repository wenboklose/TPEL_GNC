%% verification of Bar-on example:
clc
clear all
close all

k_d = -5;
k_g = 1.1;

s=tf('s');

% % system 1
% den = s^4+1.75*s^3+7.5*s^2+4*s+8;
% num11 = 0.0625*s+0.25;
% num12 = s^2+s+4;
% num21 = 0.25*s^2+0.1875*s+0.75;
% num22 = s+4;

% system 2
den = (s+1)*(s+2);
num11 = -47*s+2;
num12 = 56*s;
num21 = -42*s;
num22 = 50*s+2;
k1 = 0;
k2 = -0.0;
Compensator = [1+k1 0;0 1+k2];
% num11 = s^2+s+4;
% num12 = 0.0625*s+0.25;
% num21 = 0.25*s^2+0.1875*s+0.75;
% num22 = s+4;
L = Compensator*1/den*[num11 num12; num21 num22];
I = [tf([0 1],[0 1]), 0; 0 tf([0 1],[0 1])];
M=L/(I+L);
figure(1)
nyquist(L)
grid on
hold on
% nyquist(M)
freq_r = logspace(-5, 5,5000)*2*pi;
Resp_L = freqresp(L,freq_r);
Resp_M = freqresp(M,freq_r);
Resp_I = freqresp(I,freq_r);
for k=1:length(freq_r)
   eig_L(:,k) = eig(Resp_L(:,:,k));
   eig_L_I(:,k) = eig(Resp_L(:,:,k));
   [V,D]=eig(Resp_L(:,:,k));
    C_V(k)=dot(V(:,1),V(:,2));
    [V,D]=eig(Resp_L(:,:,k)'*Resp_L(:,:,k));
    C_V_s(k)=dot(V(:,1),V(:,2));
   Singular_L(:,k) = svd(Resp_L(:,:,k));
   [U H]=poldec(Resp_L(:,:,k));
   eig_U(:,k)=eig(U);
   eig_H(:,k)=eig(H);
   Singular_M(:,k) = svd(Resp_M(:,:,k));
   [U H]=poldec(Resp_M(:,:,k));
   eig_U_M(:,k)=eig(U);
   eig_H_M(:,k)=eig(H);
end
figure(2)
semilogx(freq_r,Singular_L(1,:),'b*');
hold on
semilogx(freq_r,Singular_L(2,:),'b*');
semilogx(freq_r,eig_H(1,:),'r-');
semilogx(freq_r,eig_H(2,:),'r-');
semilogx(freq_r,abs(eig_L(1,:)),'k');
semilogx(freq_r,abs(eig_L(2,:)),'k');
grid on
figure(3)
semilogx(freq_r,Singular_M(1,:),'b*');
hold on
semilogx(freq_r,Singular_M(2,:),'b*');
semilogx(freq_r,eig_H_M(1,:),'r-');
semilogx(freq_r,eig_H_M(2,:),'r-');
grid on
figure(4)
plot(eig_L(1,1:length(eig_L)),'b')
hold on
plot(eig_L(2,1:length(eig_L)),'r')
grid on
ezplot('x^2+y^2=1')
Mag=abs(eig_L(2,:));
Phase=angle(eig_L(2,:))*180/pi;
figure(41)
subplot(2,1,1)
semilogx(freq_r,20*log10(Mag),'b*');
grid on
subplot(2,1,2)
semilogx(freq_r,Phase,'r*');
hold off
grid on
Mag=abs(eig_L(1,:));
Phase=angle(eig_L(1,:))*180/pi;
figure(411)
subplot(2,1,1)
semilogx(freq_r,20*log10(Mag),'b*');
grid on
subplot(2,1,2)
semilogx(freq_r,Phase,'r*');
hold off
grid on


[M_max_1 index_1]=max(Singular_M(1,:));
[M_max_2 index_2]=max(Singular_M(2,:));
if M_max_1>=M_max_2
    GM = 1/M_max_1;
    freq_crit = freq_r(index_1); 
    index = index_1;
else
    GM = 1/M_max_2;
    freq_crit = freq_r(index_2);
    index = index_2;
end
[U,S,V]=svd(Resp_M(:,:,index));
um1=U(:,1);vm1=V(:,1);
PM=acosd(real(-um1'*vm1));
% GM = 1/M_max_2
%% now, form a perturbation which will destabilize the system
D = -1*GM*vm1*um1';
% D=[0.0822+0.1281*j -0.0125-0.0960*j;0.0141-0.1191*j 0.1122-0.1733*j];
[U R]=poldec(D);
eig(R);
[S V]=eig(U);
angle(eig(U))*180/pi;
max(abs(angle(eig(U))))*180/pi;

figure(42)
plot(eig_L_I(1,1:length(eig_L)),'b')
hold on
plot(eig_L_I(2,1:length(eig_L)),'r')
grid on
plot(conj(eig_L_I(1,1:length(eig_L))),'b+')
plot(conj(eig_L_I(2,1:length(eig_L))),'r+')
for k=1:length(freq_r)
   eig_L(:,k) = eig(Resp_L(:,:,k)*(Resp_I(:,:,k)+D));
end
[eig_L]=sortloci(eig_L);

plot(eig_L(1,1:length(eig_L)),'r')
hold on
plot(eig_L(2,1:length(eig_L)),'b')
grid on
plot(conj(eig_L(1,1:length(eig_L))),'r*')
plot(conj(eig_L(2,1:length(eig_L))),'b*')


Mag=abs(eig_L(2,:));
Phase=angle(eig_L(2,:))*180/pi;
figure(43)
subplot(2,1,1)
semilogx(freq_r,20*log10(Mag),'b*');
grid on
subplot(2,1,2)
semilogx(freq_r,Phase,'r*');
hold off
grid on

figure(44)
step(M)
grid on

% 
% D_1=[0.0822+0.1281*j -0.0125-0.0960*j;0.0141-0.1191*j 0.1122-0.1733*j]
% [U R]=poldec(D_1)
% eig(R)
% [S V]=eig(U)
% angle(eig(U))*180/pi
% max(abs(angle(eig(U))))*180/pi
% 
% for k=1:length(freq_r)
%    eig_L(:,k) = eig(Resp_L(:,:,k)*(Resp_I(:,:,k)+D_1));
% end
% % figure(42)
% plot(eig_L(1,1:length(eig_L)),'b*')
% % hold on
% plot(eig_L(2,1:length(eig_L)),'r')
% % grid on




%% this part is for multiply perturbation: trying to define phase and gain margin using eigenvalue of matrix L
% first, get the magnitude and angle of L's eigenvalue:
for k=1:length(freq_r)
   eig_L(:,k) = eig(Resp_L(:,:,k));
end

angle_eigs_L=angle(eig_L);
Mag_eigs_L=abs(eig_L);

% for phase margin, find at which frequency Mag of L' eig reach unity
% circle

index_1=find(Mag_eigs_L(1,:)>0.999&Mag_eigs_L(1,:)<1.002);
index_2=find(Mag_eigs_L(2,:)>0.99&Mag_eigs_L(2,:)<1.01);
PM_1=(angle_eigs_L(1,index_1)+pi)*180/pi;
PM_2=(angle_eigs_L(2,index_2)+pi)*180/pi;

if isempty(PM_1)
    [PM Index]= min(PM_2);
    Index = index_2;
elseif isempty(PM_2)
    [PM Index] = min(PM_1);
    Index = index_1;
elseif min(PM_1)>=min(PM_2)
    [PM Index] = min(PM_2);
    Index = index_2;
else
    [PM Index] = min(PM_1);
    Index = index_1;
end
freq_r(Index)/(2*pi)
Phase_Margin=PM
eig_L(:,Index)
for k=1:length(freq_r)
   eig_L(:,k) = eig(Resp_L(:,:,k));
end
% [eig_L]=sortloci(eig_L);
figure(5)
plot(eig_L(1,1:length(eig_L)),'b')
hold on
plot(eig_L(2,1:length(eig_L)),'r')
grid on
ezplot('x^2+y^2=1')

U=[exp(-1i*PM*pi/180) 0;0 exp(-1i*PM*pi/180)];
for k=1:length(freq_r)
   eig_L_1(:,k) = eig(Resp_L(:,:,k)*(U));
end
[eig_L_1]=sortloci(eig_L_1);
% figure(41)
plot(eig_L_1(1,1:length(eig_L_1)),'b*')
% hold on
plot(eig_L_1(2,1:length(eig_L_1)),'r*')
% grid on
U_1=[exp(-1*(PM-k_d)*pi/180*s/freq_r(index)) 0;0 exp(-1*(PM-k_d)*pi/180*s/freq_r(index))];


% step response:
L_1=L*U_1;
M_1=L_1/(I+L_1);
figure(55)
step(M)
hold on
step(M_1)
grid on


U_2=[exp(-1*(PM+k_d)*pi/180*s/freq_r(index)) 0;0 exp(-1*(PM+k_d)*pi/180*s/freq_r(index))];


% step response:
L_2=L*U_2;
M_2=L_2/(I+L_2);
figure(551)
step(M)
hold on
step(M_2)
grid on


U_1=[exp(-1i*(PM-k_d)*pi/180) 0;0 exp(-1i*(PM-k_d)*pi/180)];
figure(51)
for k=1:length(freq_r)
   eig_L_2(:,k) = eig(Resp_L(:,:,k)*U_1);
end
[eig_L_2]=sortloci(eig_L_2);
U_2=[exp(-1i*(PM+k_d)*pi/180) 0;0 exp(-1i*(PM+k_d)*pi/180)];
figure(51)
for k=1:length(freq_r)
   eig_L_3(:,k) = eig(Resp_L(:,:,k)*U_2);
end
[eig_L_3]=sortloci(eig_L_3);
% figure(41)
plot(eig_L(1,1:length(eig_L)),'b')
hold on
plot(eig_L(2,1:length(eig_L)),'r')

plot(eig_L_1(1,1:length(eig_L_1)),'b*')
hold on
plot(eig_L_1(2,1:length(eig_L_1)),'r*')
grid on

plot(eig_L_2(1,1:length(eig_L_2)),'bo')
hold on
plot(eig_L_2(2,1:length(eig_L_2)),'ro')
grid on

plot(eig_L_3(1,1:length(eig_L_2)),'bo')
hold on
plot(eig_L_3(2,1:length(eig_L_2)),'ro')
grid on
ezplot('x^2+y^2=1')
% for gain margin, find at which frequency phase of L' eig reach -pi

index_1=find(angle_eigs_L(1,:)>-pi-0.02&angle_eigs_L(1,:)<-pi+0.02);
index_2=find(angle_eigs_L(2,:)>-pi-0.02&angle_eigs_L(2,:)<-pi+0.02);
GM_1=(1./Mag_eigs_L(1,index_1));
GM_2=(1./Mag_eigs_L(2,index_2));

if isempty(GM_1)
    [GM Index] = min(GM_2);
    Index = index_2;
elseif isempty(GM_2)
    [GM Index] = min(GM_1);
    Index = index_1;
elseif min(GM_1)>=min(GM_2)
    [GM Index] = min(GM_2);
    Index = index_2;
else
    [GM Index] = min(GM_1);
    Index = index_1;
end
freq_r(Index)/(2*pi)
Gain_Margin=20*log10(GM)
eig_L(:,Index)
for k=1:length(freq_r)
   eig_L(:,k) = eig(Resp_L(:,:,k));
end
% [eig_L]=sortloci(eig_L);
figure(6)
plot(eig_L(1,1:length(eig_L)),'b')
hold on
plot(eig_L(2,1:length(eig_L)),'r')
grid on
ezplot('x^2+y^2=1')
R=[GM 0;0 GM];
for k=1:length(freq_r)
   eig_L_1(:,k) = eig(Resp_L(:,:,k)*(R));
end
[eig_L_1]=sortloci(eig_L_1);
% figure(41)
plot(eig_L_1(1,1:length(eig_L_1)),'b*')
hold on
plot(eig_L_1(2,1:length(eig_L_1)),'r*')
grid on

R_1=[k_g*GM 0;0 k_g*GM];
% step response:
L_2=L*R_1;
M_2=L_2/(I+L_2);
figure(66)
step(M)
hold on
step(M_2)
grid on

R_2=[1/k_g*GM 0;0 1/k_g*GM];
% step response:
L_3=L*R_2;
M_3=L_3/(I+L_3);
figure(661)
step(M)
hold on
step(M_3)
grid on


figure(61)
for k=1:length(freq_r)
   eig_L_2(:,k) = eig(Resp_L(:,:,k)*R_1);
end
[eig_L_2]=sortloci(eig_L_2);

R_2=[1/k_g*GM 0;0 1/k_g*GM];
for k=1:length(freq_r)
   eig_L_3(:,k) = eig(Resp_L(:,:,k)*R_2);
end
[eig_L_2]=sortloci(eig_L_2);

% figure(41)
plot(eig_L(1,1:length(eig_L)),'b')
hold on
plot(eig_L(2,1:length(eig_L)),'r')

plot(eig_L_1(1,1:length(eig_L_1)),'b*')
hold on
plot(eig_L_1(2,1:length(eig_L_1)),'r*')
grid on

plot(eig_L_2(1,1:length(eig_L_2)),'bo')
hold on
plot(eig_L_2(2,1:length(eig_L_2)),'ro')
grid on

plot(eig_L_3(1,1:length(eig_L_2)),'bo')
hold on
plot(eig_L_3(2,1:length(eig_L_2)),'ro')
grid on
ezplot('x^2+y^2=1')