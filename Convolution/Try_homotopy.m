%Skipt zur Simulation mit Homotopie
clear Xem
close all
% %Range
x=linspace(-0.5,1.5);
%e.g. 
poth =@(y) 8.*y.^4-44/3.*y.^3 + 2.*y.^2 +11/3.*y+1;
dpoth =@(y) 32.*y.^3-44.*y.^2+4.*y.^2+11/3;
pot = poth(x);
dpot= dpoth(x);
%t=0:0.1:1;
t=0.16:-.02:0.0;
%homopot = homopoth(x);





epsilon = 0.3;
% a = 0.5;
% b =-0.5;
% Xzero=1;

a = 1.1;
b = 1.2;


N=1000; Xzero=-0.25; dt=1/N;

dW=sqrt(dt)*randn(1,N);
%Xem = zeros(1,N);



 
 for i=1:length(t)
     clear Xem
     homopoth =@(y) poth(y)+96.*t(i).^2+(4-88.*y+96.*y.^2).*t(i);
     dhomopoth =@(y) dpoth(y)+(-88 +192.*y).*t(i);
     Xtemp=Xzero;
     j=1;
    %for j=1:N
    while (Xtemp < 0.5)
        Xtemp = Xtemp -(dhomopoth(Xtemp))*dt + sqrt(dt)*epsilon*randn;
        Xem(j) = Xtemp;
        if j >=N
              break
        end
        j=j+1;
 %end
        %j=j+1
%         if ((Xtemp>b) || (Xtemp<a))
%             hit = hit+1;
%             break
%         end
    end
%hit 



figure(i)
% subplot(121)
% plot(x,homopoth(x))
% subplot(122)
plot(x,homopoth(x)), hold on
plot([Xzero,Xem],homopoth([Xzero Xem]),'r+'),%hold off
hold off
legend('Homotopy(Pot)','Homotopy(Pot(X_SDE))')
title('asymetisches 2 Well Potential')
xlabel(['Anzahl der Schritte ' num2str(j)])

 end
