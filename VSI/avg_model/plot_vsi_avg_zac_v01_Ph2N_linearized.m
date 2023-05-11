H=tf(Model_sys);

Vavg1=[H(1,1) H(1,2);
    H(2,1) H(2,2);];

Ilavg1=[H(3,1) H(3,2);
    H(4,1) H(4,2);];

Isavg1=[H(5,1) H(5,2);
    H(6,1) H(6,2);];
Zsavg1=Vavg1/Isavg1; Zlavg1=Vavg1/Ilavg1;

% H=tf(Model2_sys);
% 
% Vavg2=[H(1,1) H(1,2);
%     H(2,1) H(2,2);];
% 
% Ilavg2=[H(3,1) H(3,2);
%     H(4,1) H(4,2);];
% 
% Isavg2=[H(5,1) H(5,2);
%     H(6,1) H(6,2);];
% Zsavg2=Vavg2/Isavg2; Zlavg2=Vavg2/Ilavg2;

% wmin=2*pi*1; wmax=2*pi*10e3;
% figure(1);
% bode(Zsavg1,{wmin, wmax},Zsavg2,{wmin, wmax});
% 
% figure(2);
% bode(Zlavg1,{wmin, wmax},Zlavg2,{wmin, wmax});


fresp=40:1:10000; 
% fresp=fplot;
wresp=2*pi*fresp;
Vfr=freqresp(Vavg1,wresp); Ilfr=freqresp(Ilavg1,wresp); Isfr=freqresp(Isavg1,wresp); 
Nfr=length(fresp); Zsfr=zeros(2,2,Nfr); Zlfr=zeros(2,2,Nfr);
for iterator=1:length(Vfr),
    Vfrit=[Vfr(1,1,iterator) Vfr(1,2,iterator);
        Vfr(2,1,iterator) Vfr(2,2,iterator);];
    Ilfrit=[Ilfr(1,1,iterator) Ilfr(1,2,iterator);
        Ilfr(2,1,iterator) Ilfr(2,2,iterator);];
    Isfrit=[Isfr(1,1,iterator) Isfr(1,2,iterator);
        Isfr(2,1,iterator) Isfr(2,2,iterator);];
    Zsfr(:,:,iterator)=Vfrit/Isfrit;
    Zlfr(:,:,iterator)=Vfrit/Ilfrit;
end

linewidth=4; fontsize=24;

linestyle1='--b'; linestyle2='-g';
linestyle3='--r'; linestyle4='-m';


figure(5);
clf;
subplot(2,1,1);
semilogx(fresp,20*log10(squeeze(abs(Zsfr(1,1,:)))),linestyle1,'LineWidth',linewidth);
hold on;
semilogx(fresp,20*log10(squeeze(abs(Zsfr(2,2,:)))),linestyle2,'LineWidth',linewidth);
hold off;
grid on;
xlabel('frequency [Hz]','FontSize',fontsize);
ylabel('magnitude [dB]','FontSize',fontsize);
set(gca,'FontSize',fontsize);
subplot(2,1,2);
semilogx(fresp,unwrap(phase(squeeze(Zsfr(1,1,:))))*180/pi,linestyle1,'LineWidth',linewidth);
hold on;
semilogx(fresp,unwrap(phase(squeeze(Zsfr(2,2,:))))*180/pi,linestyle2,'LineWidth',linewidth);
hold off;
grid on;
xlabel('frequency [Hz]','FontSize',fontsize);
ylabel('phase [deg]','FontSize',fontsize);
set(gca,'FontSize',fontsize);


figure(6);
clf;
subplot(2,1,1);
semilogx(fresp,20*log10(squeeze(abs(Zsfr(1,2,:)))),linestyle1,'LineWidth',linewidth);
hold on;
semilogx(fresp,20*log10(squeeze(abs(Zsfr(2,1,:)))),linestyle2,'LineWidth',linewidth);
hold off;
grid on;
xlabel('frequency [Hz]','FontSize',fontsize);
ylabel('magnitude [dB]','FontSize',fontsize);
set(gca,'FontSize',fontsize);
subplot(2,1,2);
semilogx(fresp,unwrap(phase(squeeze(Zsfr(1,2,:))))*180/pi,linestyle1,'LineWidth',linewidth);
hold on;
semilogx(fresp,unwrap(phase(squeeze(Zsfr(2,1,:))))*180/pi,linestyle2,'LineWidth',linewidth);
hold off;
grid on;
xlabel('frequency [Hz]','FontSize',fontsize);
ylabel('phase [deg]','FontSize',fontsize);
set(gca,'FontSize',fontsize);


figure(7);
clf;
subplot(2,1,1);
semilogx(fresp,20*log10(squeeze(abs(Zlfr(1,1,:)))),linestyle3,'LineWidth',linewidth);
hold on;
semilogx(fresp,20*log10(squeeze(abs(Zlfr(2,2,:)))),linestyle4,'LineWidth',linewidth);
hold off;
grid on;
xlabel('frequency [Hz]','FontSize',fontsize);
ylabel('magnitude [dB]','FontSize',fontsize);
set(gca,'FontSize',fontsize);
subplot(2,1,2);
semilogx(fresp,unwrap(phase(squeeze(Zlfr(1,1,:))))*180/pi,linestyle3,'LineWidth',linewidth);
hold on;
semilogx(fresp,unwrap(phase(squeeze(Zlfr(2,2,:))))*180/pi,linestyle4,'LineWidth',linewidth);
hold off;
grid on;
xlabel('frequency [Hz]','FontSize',fontsize);
ylabel('phase [deg]','FontSize',fontsize);
set(gca,'FontSize',fontsize);


figure(8);
clf;
subplot(2,1,1);
semilogx(fresp,20*log10(squeeze(abs(Zlfr(1,2,:)))),linestyle3,'LineWidth',linewidth);
hold on;
semilogx(fresp,20*log10(squeeze(abs(Zlfr(2,1,:)))),linestyle4,'LineWidth',linewidth);
hold off;
grid on;
xlabel('frequency [Hz]','FontSize',fontsize);
ylabel('magnitude [dB]','FontSize',fontsize);
set(gca,'FontSize',fontsize);
subplot(2,1,2);
semilogx(fresp,unwrap(phase(squeeze(Zlfr(1,2,:))))*180/pi,linestyle3,'LineWidth',linewidth);
hold on;
semilogx(fresp,unwrap(phase(squeeze(Zlfr(2,1,:))))*180/pi,linestyle4,'LineWidth',linewidth);
hold off;
grid on;
xlabel('frequency [Hz]','FontSize',fontsize);
ylabel('phase [deg]','FontSize',fontsize);
set(gca,'FontSize',fontsize);


ZVSI=Zsfr;
ZoutVSI=Zsavg1;
% save('ZoutVSI_o_loop.mat','ZoutVSI')
% save('ZoutVSI_o_loop_delay.mat','ZoutVSI')
% save ('ZoutVSI_i_loop_no_decoup.mat','ZoutVSI')
% save ('ZoutVSI_i_loop_no_decoup_delay.mat','ZoutVSI')
% save ('ZoutVSI_v_loop_no_decoup.mat','ZoutVSI')
% save ('ZoutVSI_v_loop_no_decoup_delay.mat','ZoutVSI')

