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
%t=0.16:-.02:0.0;
%homopot = homopoth(x);
beta=0.05;
l= [0,0.08];




epsilon = 0.3;
% a = 0.5;
% b =-0.5;
% Xzero=1;

a = 1.1;
b = 1.2;


N=1000; Xzero=-0.25; dt=1/1000;

dW=sqrt(dt)*randn(1,N);
%Xem = zeros(1,N);


S=1000; % Number of Simulations
homopoth =@(y,l) poth(y)+96.*l^2+(4-88.*y+96.*y.^2).*l;
dhomopoth =@(y,l) dpoth(y)+(-88 +192.*y).*l;
Xem1=zeros(1,S);
Xem2=zeros(1,S);

X1=Xzero;
X2=Xzero;
 
for k=1:S
     clear Xem
     
  
    %for j=1:N
    
     X1 = X1 -(dhomopoth(X1,l(1)))*dt + sqrt(dt)*epsilon*randn;
     Xem1(k) = X1;
     X2 = X2 -(dhomopoth(X2,l(2)))*dt + sqrt(dt)*epsilon*randn;
     Xem2(k) = X2;
 
     alpha = min(1,((exp(-beta *homopoth(X1,l(1)))*exp(-beta*homopoth(X2,l2)))/...
         (exp(-beta *homopoth(X1,l(2)))*exp(-beta*homopoth(X2,l1)))));
     u=rand;
     
     if (u<alpha)
         lh=l(1);
         l(1)=l(2);
         l(2)=lh;
     end
         
end

  %für feste Anzhal an schritten schauen wie viele das ziel erreicht haben  


% figure(i)
% % subplot(121)
% % plot(x,homopoth(x))
% % subplot(122)
% plot(x,homopoth(x)), hold on
% plot([Xzero,Xem],homopoth([Xzero Xem]),'r+'),%hold off
% hold off
% legend('Homotopy(Pot)','Homotopy(Pot(X_SDE))')
% title('asymetisches 2 Well Potential')
% xlabel(['Anzahl der Schritte ' num2str(j)])

 
 
 hitting_time = hitting_time./S;
 %try what so ever
 % ist einfach nur ein Erwartunsgwert der durch 1000 geteilt wird. 
 % kann man nciht als Wahrscheinlichkeit interpretieren
 ht = hitting_time./N;
 ht = 1-ht;
 
 hitting_time_prob =  hitting_time_prob./S;
 hitting_time_prob2 =  hitting_time_prob2./S;
 hitting_time_prob3 =  hitting_time_prob3./S;
 
 
 
 
 figure(1)
 plot(t,hitting_time,'-o')
 xlabel('Homotopieparameter')
 ylabel('Anzahl der Schritte')
 title(['Mittel der Schrittzahl ueber ' num2str(S) ' Simulationen'])
 
 figure(2)
 plot(t,hitting_time_prob,'-o')
 xlabel('Homotopieparameter')
 ylabel('Wahrscheinlichkeit')
 title(['Exit Prob in ' num2str(N/10) ' Steps '])
 
 figure(3)
 plot(t,hitting_time_prob2,'-o')
 xlabel('Homotopieparameter')
 ylabel('Wahrscheinlichkeit')
 title(['Exit Prob in ' num2str(N/2) ' Steps '])
 
 figure(4)
 plot(t,hitting_time_prob3,'-o')
 xlabel('Homotopieparameter')
 ylabel('Wahrscheinlichkeit')
  title(['Exit Prob in ' num2str(N) ' Steps '])
  
 figure(5) 
 plot(t,ht,'-o')
 xlabel('Homotopieparameter')
 ylabel('Wahrscheinlichkeit')
 title(['Erwartungswert/1000 Vorsicht bei der Interpretation'])