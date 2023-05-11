% Function gershgorinplot(A,fmin,fmax,res,step,fig,split) plots the
% Nyquist Array of Matrix A superimposing the Gershgorin bands over 
% the Nyquist plots of elements aii(s). It uses figure (fig), frequency
% range given by 10^dec1 to 10^dec2 in Hz, and its resolution by res.
% It plots Gershgorin discs every 'step' samples of the frequency vector.
% If split is chosen the four channels are plotted in separate windows.
%
% Call: gershgorinplot(A,fmin,fmax,res,step,fig,split)

function gershgorinplot(A,Eig_1,Eig_2,k_freq_r,step,fig)

% n=res;
% % Frequency range in Hz
% f=logspace(fmin,fmax,n);
% w=f*2*pi;
% if(1)
%     w=[-fliplr(w) w];
%     n=n*2;
% end
% 
% % Transfer Function Extraction for matrix A
% A=tf(A);
% a11=A(1,1);
% a12=A(1,2);
% a21=A(2,1);
% a22=A(2,2);
% 
% % A(s) frequency response
% Aresp=freqresp(A,w);

% Aii(s) magnitude in dB
a11real(1,:)=real(A(1,1).ResponseData(1:length(A(1,1).ResponseData)*k_freq_r));
a12real(1,:)=real(A(1,2).ResponseData(1:length(A(1,1).ResponseData)*k_freq_r));
a21real(1,:)=real(A(2,1).ResponseData(1:length(A(1,1).ResponseData)*k_freq_r));
a22real(1,:)=real(A(2,2).ResponseData(1:length(A(1,1).ResponseData)*k_freq_r));

% Aii(s) phase in deg
a11imag(1,:)=imag(A(1,1).ResponseData(1:length(A(1,1).ResponseData)*k_freq_r));
a12imag(1,:)=imag(A(1,2).ResponseData(1:length(A(1,1).ResponseData)*k_freq_r));
a21imag(1,:)=imag(A(2,1).ResponseData(1:length(A(1,1).ResponseData)*k_freq_r));
a22imag(1,:)=imag(A(2,2).ResponseData(1:length(A(1,1).ResponseData)*k_freq_r));

% Gershgorin Disc Radii
radii(1,:)=abs(A(1,2).ResponseData(1:length(A(1,1).ResponseData)*k_freq_r));
radii(2,:)=abs(A(2,1).ResponseData(1:length(A(1,1).ResponseData)*k_freq_r));

figHandle=figure(fig)
set(figHandle, 'Position', [50, 10, 1000, 800]);
    clf
    lw=1;
    % Axis determination
    plot(a22real,a22imag,'LineWidth',[lw])
    hold on
    gdiscs(a22real,a22imag,radii,step)
    hold off
    axis22=axis;
    plot(a21real,a21imag,'LineWidth',[lw])
    axis21=axis;
    plot(a12real,a12imag,'LineWidth',[lw])
    axis12=axis;
    % A11 Bode
    subplot(221)
    set(gca,'Position',[0.075 0.515 .415 .405])
    plot(a11real,a11imag,'LineWidth',[lw])
    hold on
    gdiscs(a11real,a11imag,radii(1,:),step)
    plot(-1,0,'r+','LineWidth',[2]);
    plot(Eig_1,'r');
    legend({'{\it\lambda}_{1}(L)\it_{dq}'},'Fontsize',[12],'FontWeight','bold')
    hold off
    axis11=axis;
    taxis=[axis11; axis12; axis21; axis22];
    axis([min(taxis(:,1)) max(taxis(:,2)) min(taxis(:,3)) max(taxis(:,4))]);
    set(gca,'XTickLabel',[])
    grid
    title('From: D','Fontsize',[12])
    ylabel('To: D','Fontsize',[12])
    % A12 Bode
    subplot(222)
    set(gca,'Position',[0.51 0.515 .415 .405])
    plot(a12real,a12imag,'LineWidth',[lw])
    hold on
    legend({'L\it_{dq}'},'Fontsize',[12],'FontWeight','bold')
    plot(-1,0,'r+','LineWidth',[2]);
    hold off
    axis([min(taxis(:,1)) max(taxis(:,2)) min(taxis(:,3)) max(taxis(:,4))]);
    set(gca,'XTickLabel',[])
    set(gca,'YTickLabel',[])
    grid
    title('From: Q','Fontsize',[12])
    % A21 Bode
    subplot(223)
    set(gca,'Position',[0.075 0.085 .415 .405])
    plot(a21real,a21imag,'LineWidth',[lw])
    hold on
    plot(-1,0,'r+','LineWidth',[2]);
    hold off
    axis([min(taxis(:,1)) max(taxis(:,2)) min(taxis(:,3)) max(taxis(:,4))]);
    grid
    ylabel('To: Q','Fontsize',[12])
    % A22 Bode
    subplot(224)
    set(gca,'Position',[0.51 0.085 .415 .405])
    plot(a22real,a22imag,'LineWidth',[lw])
    hold on
    gdiscs(a22real,a22imag,radii(2,:),step)
    plot(-1,0,'r+','LineWidth',[2]);
    plot(Eig_2,'r')
    legend({'{\it\lambda}_{2}(L)\it_{dq}'},'Fontsize',[12],'FontWeight','bold')
    hold off
    axis([min(taxis(:,1)) max(taxis(:,2)) min(taxis(:,3)) max(taxis(:,4))]);
    set(gca,'YTickLabel',[])
    grid