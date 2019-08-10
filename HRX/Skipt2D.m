%Skript für ein 2D Potential
clear all

xplot = linspace(-2,2);
yplot = linspace(-1,2);

[Xplot,Yplot] = meshgrid(xplot,yplot);

V = @(x,y) 3.*exp(-x.^2-(y-1/3).^2)-3.*exp(-x.^2-(y-5/3).^2) - 5.*exp(-(x-1).^2-y.^2)-5.*exp(-(x+1).^2-y.^2) + ...
    1/5.*x.^4 +1/5*(y-1/3).^4;

dxV =@(x,y) -6.*x.*exp(-x.^2-(y-1/3).^2)+6.*x.*exp(-x.^2-(y-5/3).^2) + 10.*(x-1).*exp(-(x-1).^2-y.^2)+10.*(x+1).*exp(-(x+1).^2-y.^2) + ...
    4/5.*x.^4;

dyV=@(x,y) -6.*(y-1/3).*exp(-x.^2-(y-1/3).^2)+6.*(y-5/3).*exp(-x.^2-(y-5/3).^2) +10.*y.*exp(-(x-1).^2-y.^2)+10.*y.*exp(-(x+1).^2-y.^2) + ...
     +4/5*(y-1/3).^3;

% figure(1)
% surf(Xplot,Yplot,V(Xplot,Yplot))

epsilon = 0.1;

N=1000; Xzero=1; dt=1/N;

dW=sqrt(dt)*randn(1,N);
%Xem = zeros(1,N);
Xzerox = 0;
Xzeroy = 0.5;
Xtempx = Xzerox;
Xtempy = Xzeroy;
j=1;

anzSample=1000;
 %while (Xtemp>b) || (Xtemp<a)
 for j=1:anzSample
     %Winc = sum(dW((j-1)+1:j));
     %Xtemp = Xtemp -(dpoth(Xtemp) + dbasisfun(Xtemp))*dt + sqrt(dt)*epsilon*randn; 
     Xtempx = Xtempx -(dxV(Xtempx,Xtempy))*dt + sqrt(dt)*epsilon*randn;%+1/2*(sqrt(dt)*epsilon*randn)^2-dt;
     Xtempy = Xtempy -(dyV(Xtempx,Xtempy))*dt + sqrt(dt)*epsilon*randn;%+1/2*(sqrt(dt)*epsilon*randn)^2-dt;
     Xemx(j) = Xtempx;
     Xemy(j) = Xtempy;
 %end
 %j=j+1;
 end
 
figure(1)
surf(Xplot,Yplot,V(Xplot,Yplot)), hold on
plot3([Xzerox Xemx],[Xzeroy Xemy],V([Xzerox Xemx],[Xzeroy Xemy]),'r+'),%hold off
hold off
legend('Pot','Pot(X_{SDE}))')
title('asymetisches 3 Well Potential')