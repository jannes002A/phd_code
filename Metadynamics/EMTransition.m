%Skipt zur Simulation mit Metadynamics mit Girsanov
clear Xem
close all
% %Range
x=linspace(-0.5,1.5);
%e.g. 
poth =@(y) 8.*y.^4-44/3.*y.^3 + 2.*y.^2 +11/3.*y+1;
dpoth =@(y) 32.*y.^3-44.*y.^2+4.*y.^2+11/3;
pot = poth(x);
dpot= dpoth(x);





epsilon = 0.3;
% a = 0.5;
% b =-0.5;
% Xzero=1;

a = 1.1;
b = 1.2;


N=100000; Xzero=-0.25; dt=1/N;

dW=sqrt(dt)*randn(1,N);
%Bias Potential

dVbias=0;
sigma=0.2;
tau=1;
omega=1/2000;


 
     Xtemp=Xzero;
     j=1;
    %for j=1:N
    while (Xtemp < 1)
    

        
        Xtemp = Xtemp -(dpoth(Xtemp))*dt + sqrt(dt)*epsilon*randn;
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



% figure(1)
% % subplot(121)
% % plot(x,homopoth(x))
% % subplot(122)
% plot(x,poth(x)), hold on
% plot([Xzero,Xem],poth([Xzero Xem]),'r+'),%hold off
% hold off
% legend('Pot','Pot(X_SDE)')
% title('asymetisches 2 Well Potential')
% xlabel(['Anzahl der Schritte ' num2str(j)])
% 
% figure(2)
% Bx = Basisfunc(x,j-1,Xem(1:j-1),sigma,1); 
% dVbias = tau*omega*sum(Bx,2);
% plot(x,-dpoth(x)'+sqrt(2)*dVbias)
% 
% VbiasP=zeros(1,length(x));
% for k=1:j-1
%     VbiasP = VbiasP - 1/2*erf((Xem(k)-x)./(sqrt(2)*sigma));
% end
% 
% figure(3)
% plot(x,poth(x)+sqrt(2)*tau*omega*VbiasP);hold on
% plot(x,poth(x));hold off
% title('Pot+Vbias')
% 
% figure(4)
% plot(x,poth(x));hold on 
% plot(x,sqrt(2)*tau*omega*VbiasP);hold off


