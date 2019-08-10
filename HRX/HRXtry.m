%Skipt zur Simulation mit replica exchange with homotopy
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
%t=0.16:-.02:0.0;
%homopot = homopoth(x);
beta=4;
l= 0.04;
%l=[0,0];






 Xzero=-0.25; dt=1/1000;




S=10000; % Number of Simulations
homopoth =@(y,l) poth(y)+96.*l^2+(4-88.*y+96.*y.^2).*l;
dhomopoth =@(y,l) dpoth(y)+(-88 +192.*y).*l;
Xem1=zeros(1,S);
Xem2=zeros(1,S);

X1=Xzero;
X2=Xzero;
count=zeros(1,S);
swap=0;

for k=1:S
     clear Xem
     
  
    %for j=1:N
    
     X1 = X1 -(dpoth(X1))*dt + sqrt(dt*2/beta)*randn;
     
     X2 = X2 -(dhomopoth(X2,l))*dt + sqrt(dt*2/beta)*randn;
     
     
     if swap == 0
        Xem1(k) = X1;
        Xem2(k) = X2;
     else
        Xem1(k) = X2;
        Xem2(k) = X1;
     end
 
     alpha = min(1,((exp(-beta *poth(X1))*exp(-beta*homopoth(X2,l)))/...
         (exp(-beta *homopoth(X1,l))*exp(-beta*poth(X2))))) ;
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
 

figure(1)


hold on
plot(x,poth(x));
plot([Xzero,Xem1],poth([Xzero,Xem1]),'r+');
hold off

figure(2)

hold on
plot(x,homopoth(x,l));
plot([Xzero,Xem2],homopoth([Xzero,Xem2],l),'r+');
hold off


xp=linspace(-1,2);
figure(3)
plot(xp,exp(-beta*poth(xp)),xp,exp(-beta*homopoth(xp,l)),'Linewidth',3)

%%
p=[0.0,0.04];

figure(4)
hold on
for i=1:length(p)
    plot(x,(exp(-beta*homopoth(x,p(i)))))
end
hold off
 
 
 
 