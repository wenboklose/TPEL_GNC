% Function gdiscs(aiireal,aiimag,r,steps) plots the Gershgorin bands
% centered at (aiireal,aiimag) having radius r. It plots discs every steps
% samples of the vectors aiireal, aiimag and r which should have the same
% dimension.
%
% Call: gdiscs(aiireal,aiimag,r,steps)

function gdiscs(aiireal,aiimag,r,steps)

n=length(r);
% Angle for Disc generation
t=linspace(0,2*pi,50);

for k=steps:steps:n
    plot(r(1,k)*cos(t)+aiireal(k),r(1,k)*sin(t)+aiimag(k),'k','LineWidth',[1])
%     plot(r(2,k)*cos(t)+aiireal(k),r(2,k)*sin(t)+aiimag(k),'k','LineWidth',[1])
end



