%Plot Metastabilit

%close all
% %Range
x=linspace(-1.5,1.5);
%e.g. 
%poth =@(y) 8.*y.^4-44/3.*y.^3 + 2.*y.^2 +11/3.*y+1;
%dpoth =@(y) 32.*y.^3-44.*y.^2+4.*y+11/3;
poth =@(y)  1/2.*(y.^2-1).^2;
dpoth =@(y) 2.*y.*(y.^2-1);
pot = poth(x);
dpot= dpoth(x);

epsilon = 0.5;
cv=1000000;

N=10^6; Xzero=-0.25; dt=1/1000;

Xtemp=zeros(1,cv);
Xtemp(1)=-1;


 for i=1:cv

    Xtemp(i+1) = Xtemp(i)-(dpoth(Xtemp(i))).*dt + sqrt(dt).*epsilon.*randn(1);
                  

 end
 
figure(1)
plot(1:cv, Xtemp(1:cv))
xl= xlabel('time');
yl=ylabel('$X_t$');
set(yl,'Interpreter','Latex');
set(xl,'Interpreter','Latex');

figure(2)
plot(x,poth(x),'Linewidth',3)
xl= xlabel('x');
yl=ylabel('Potential');
set(xl,'Interpreter','Latex');
set(yl,'Interpreter','Latex');
