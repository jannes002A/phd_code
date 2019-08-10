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


N=1000; Xzero=-0.25; dt=1/1000;

dW=sqrt(dt)*randn(1,N);
%Xem = zeros(1,N);

hitting_time=zeros(1,length(t));
hitting_time_prob=zeros(1,length(t));
hitting_time_prob2=zeros(1,length(t));
hitting_time_prob3=zeros(1,length(t));
S=1000; % Number of Simulations
 
 for i=1:length(t)
     for k=1:S
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
        
        hitting_time(i)= hitting_time(i)+j; 
        
        if (j <= N/10)
            hitting_time_prob(i) =  hitting_time_prob(i)+1;
        end
        
        if (j <= N/2)
            hitting_time_prob2(i) =  hitting_time_prob2(i)+1;
        end
        if (j<N-1)
            hitting_time_prob3(i) = hitting_time_prob3(i)+1;
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

 end
 
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