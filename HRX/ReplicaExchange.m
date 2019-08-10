%Skipt zur Simulation mit replica exchange
clear Xem
close all
% %Range
x=linspace(-0.5,1.5);
%e.g. 
poth =@(y) 8.*y.^4-44/3.*y.^3 + 2.*y.^2 +11/3.*y+1;
dpoth =@(y) 32.*y.^3-44.*y.^2+4.*y.^2+11/3;
pot = poth(x);
dpot= dpoth(x);

beta1=4;
beta2=1;







Xzero=-0.25; dt=1/1000;


S=10000; % Number of Simulations


Xem1=zeros(1,S);
Xem2=zeros(1,S);

X1=Xzero;
X2=Xzero;
count=zeros(1,S);
swap=0;

for k=1:S
     clear Xem
     
  
    %for j=1:N
    
     X1 = X1 -(dpoth(X1))*dt + sqrt(dt*2/beta1)*randn;
     
     X2 = X2 -(dpoth(X2))*dt + sqrt(dt*2/beta2)*randn;
     
     
     if swap == 0
        Xem1(k) = X1;
        Xem2(k) = X2;
     else
        Xem1(k) = X2;
        Xem2(k) = X1;
     end
 
     alpha = min(1,((exp(-beta1 *poth(X1))*exp(-beta2*poth(X2)))/...
         (exp(-beta2 *poth(X1))*exp(-beta1*poth(X2))))) ;
     u=rand;
     
     count(k)= (u< alpha);
     
     if (u<alpha)
         Xh=X1;
         X1=X2;
         X2=Xh;
         swap=1;
     else
         swap=0;
     end
         
end

accrate= sum(count);
fprintf('AccRate %2.4f \n',accrate)
  %für feste Anzhal an schritten schauen wie viele das ziel erreicht haben  


figure(1)


hold on
plot(x,poth(x));
plot([Xzero,Xem1],poth([Xzero,Xem1]),'r+');
hold off

figure(2)

hold on
plot(x,poth(x));
plot([Xzero,Xem2],poth([Xzero,Xem2]),'r+');
hold off

xp=linspace(-1,2);
figure(3)
plot(xp,exp(-beta1*poth(xp)),xp,exp(-beta2*poth(xp)),'Linewidth',3) 
 

 
 
 
 