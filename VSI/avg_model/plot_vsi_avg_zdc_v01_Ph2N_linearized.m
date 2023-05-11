H=tf(Model_sys);

Vavg=H(1,1);

Ilavg=H(2,1);

Isavg=H(3,1);
Zsavg=Vavg/Isavg; Zlavg=Vavg/Ilavg;

fresp=1:1:10000; 
wresp=2*pi*fresp;
Vfr=freqresp(Vavg,wresp); Ilfr=freqresp(Ilavg,wresp); Isfr=freqresp(Isavg,wresp); 
Nfr=length(fresp); Zsfr=zeros(2,2,Nfr); Zlfr=zeros(2,2,Nfr);
for iterator=1:length(Vfr),
    Vfrit=Vfr(1,1,iterator);
    Ilfrit=Ilfr(1,1,iterator);
    Isfrit=Isfr(1,1,iterator) ;
    Zsfr(1,1,iterator)=Vfrit/Isfrit;
    Zlfr(1,1,iterator)=Vfrit/Ilfrit;
end

linewidth=4; fontsize=24;

linestyle1='-b'; linestyle2='-g';
linestyle3='-r'; linestyle4='-m';


figure(5);
clf;
subplot(2,1,1);
semilogx(fresp,20*log10(squeeze(abs(Zsfr(1,1,:)))),linestyle1,'LineWidth',linewidth);
hold on;
semilogx(fresp,20*log10(squeeze(abs(Zlfr(1,1,:)))),linestyle2,'LineWidth',linewidth);
hold off;
grid on;
xlabel('frequency [Hz]','FontSize',fontsize);
ylabel('magnitude [dB]','FontSize',fontsize);
set(gca,'FontSize',fontsize);
subplot(2,1,2);
semilogx(fresp,unwrap(phase(squeeze(Zsfr(1,1,:))))*180/pi,linestyle1,'LineWidth',linewidth);
hold on;
semilogx(fresp,unwrap(phase(squeeze(Zlfr(1,1,:))))*180/pi,linestyle2,'LineWidth',linewidth);
hold off;
grid on;
xlabel('frequency [Hz]','FontSize',fontsize);
ylabel('phase [deg]','FontSize',fontsize);
set(gca,'FontSize',fontsize);

ZinVSI=Zlavg;
% save('ZinVSI_o_loop.mat','ZinVSI')
% save('ZinVSI_o_loop_delay.mat','ZinVSI')
% save('ZinVSI_i_loop_no_decoup.mat','ZinVSI')
% save('ZinVSI_i_loop_no_decoup_delay.mat','ZinVSI')
% save('ZinVSI_v_loop_no_decoup.mat','ZinVSI')
save('ZinVSI_v_loop_no_decoup_delay.mat','ZinVSI')