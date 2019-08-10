% Skipt zur Simulation mit Metadynamics mit Girsanov
% To do;
% Sampling and Reweigting
% Proof of Novikov condition
% How to put the dVBias maybe there is a clever way
%  

%prevS = rng(0,'v5normal'); % use legacy generator, save previous settings
rng(1,'twister')                 % restore the previous settings


% %Range
x=linspace(-1.5,1.5);
%e.g. 
% poth =@(y) 8.*y.^4-44/3.*y.^3 + 2.*y.^2 +11/3.*y+1;
% dpoth =@(y) 32.*y.^3-44.*y.^2+4.*y.^2+11/3;
% pot = poth(x);
% dpot= dpoth(x);
% Xzero=-0.25;

poth =@(y)  1/2.*(y.^2-1).^2;
dpoth =@(y) 2.*y.*(y.^2-1);
pot = poth(x);
dpot= dpoth(x);
Xzero=-1;





Xtemp=Xzero;

Ns=0000;
Tl=15000;
dt=0.0001;
Xm=zeros(1,1000);
G=zeros(1,Ns);
G2=zeros(1,Ns);
Z2=zeros(1,Ns);
time=ones(1,Ns)*Tl*dt;
ht=zeros(1,Ns);
pt=zeros(1,Ns);
j=0;

sigma=0.8*ones(1,1000);
tau=1;
omega=0.1;
count=0;
beta=3;



% eine Trajektorie um das Bias Potential zu berechnen. 

 %for j=1:Tl
 while( (Xtemp > 0.9 && Xtemp < 1.1)==0 )
 %while( (Xtemp < -0.9 && Xtemp > -1.1)==0 )
    
        if j>1 && mod(j,200)== 0
            count=count+1;
            Xm(count) = Xtemp;
        end
        if count >= 1000
            error('Trajektorie zu lang!')
        end
        
        Bx = Basisfunc(Xtemp,count,Xm(1:count),sigma(1:count),0);
        dVbias = tau*omega*ones(1,count) * Bx';
        
        
        dBt=sqrt(2*dt/beta)*randn;
        Xtemp = Xtemp + (dVbias-dpoth(Xtemp))*dt + dBt;


%         if Xtemp > 1.1 && Xtemp <1.2
%             break;
%         end
        
        j=j+1;
        
 end
 fprintf('Steps needed in the first run %8.2f \n',j)
% Sampling mit dem fixen Bias Potential   
 for i=1:Ns
     
     Xtemp=zeros(1,Tl+1);
     Gs=zeros(1,Tl+1);
     Gd=zeros(1,Tl+1);
     Gh=zeros(1,Tl+1);
     Xtemp(1)=Xzero;
     first=0;
     Zt=zeros(1,Tl+1);
     zt(1)=1;
    
    for j=1:Tl
    
%         if j>1 && mod(j,200)== 0
%             count=count+1;
%             Xm(count) = Xtemp(j);
%         end
        
        Bx = exp(-0.5*(Xtemp(j)-Xm(1:count)).^2./sigma(1:count).^2)./(sigma(1:count).*sqrt(2*pi));
        Bx = sum(Bx);
        dVbias = tau*omega* Bx;
        
        
        dBt=randn;
        Xtemp(j+1) = Xtemp(j) + (dVbias-dpoth(Xtemp(j)))*dt + sqrt(2*dt/beta)*dBt;
        Gs(j+1)= Gs(j) - (dVbias/sqrt(2/beta))*sqrt(dt)*dBt; 
        Gd(j+1)= Gd(j) - 1/2* (dVbias/sqrt(2/beta))^2*dt;
        %Vorzeichen checken!
        Zt(j+1) = Zt(j) + (dVbias/sqrt(2/beta))*sqrt(dt)*dBt;   

         if Xtemp(j+1)>0.9 && Xtemp(j+1)<1.1 && first==0
         %if Xtemp(j+1)<-0.9 && Xtemp(j+1)>-1.1 && first==0
             first=1;
             ht(i)=j+1;
         end

%          if Xtemp(j+1)>1.1 && Xtemp(j+1)<1.2
%                Gh(count2)=(dVbias/sqrt(2/beta))*sqrt(dt)*dBt- 1/2* (dVbias/sqrt(2/beta))^2*dt;
%                count2=count2+1;
%          end
        
        
    end
    
    if ht(i)== 0
         G2(i)= 0;
         Z2(i)= 0;
    else
        time(i) = ht(i)*dt; 
        G2(i)= exp(Gs(ht(i))+Gd(ht(i)));
        Z2(i)= Zt(ht(i));
        
    end
    G(i) = exp(Gs(Tl)+Gd(Tl));
    pt(i)= exp(-beta*ht(i)*dt)*G2(i);
    
 end
%hit 
fprintf('Girsanov Reweighting fixed Parameter \n')
fprintf('Mittelwert G(N): %2.8f \n',mean(G))
fprintf('Mittelwert G(ht): %2.8f \n',mean(G2))
fprintf('Var(E[-beta*tau*Z]): %2.8f \n',var(G2))
fprintf('R(I): %2.8f \n',sqrt(var(G2))/mean(G2))

fprintf('E[-beta*tau*Z]: %2.8f \n',mean(pt))
fprintf('Var(E[-beta*tau*Z]): %2.8f \n',var(pt))
fprintf('R(I): %2.8f \n',sqrt(var(pt))/mean(pt))

fprintf('E[time]: %2.8f \n',mean(time))




x=linspace(-1.5,1.5);
Bx = Basisfunc(x,count,Xm(1:count),sigma,1);
dVbias = tau*omega*ones(1,count) * Bx';

figure(1)
plot(x,dVbias,'g',x, -dpot,'b',x,dVbias-dpot,'r','LineWidth',3); 
l=legend('$\nabla V_{bias}$','$-\nabla V$','$\nabla V_{bias} - \nabla V$');
set(l,'Interpreter','Latex');          


 %%   
% sigma=1.6;
% tau=1;
% omega=0.2;
% count=0;
% beta=3;
% 
%  
%         figure(2)
%         x=linspace(-1.5,1.5);
%         Bx = Basisfunc(x,1,1,sigma(1),1);
%         dVbias = tau*omega*ones(1,1) * Bx';
%         plot(x,dVbias)
%         
% figure(1)
% 
% plot(x,poth(x)), hold on
% plot(Xtemp,poth(Xtemp),'r+'),%hold off
% hold off
% legend('Pot','Pot(X_SDE)')
% title('asymetisches 2 Well Potential')
% xlabel(['Anzahl der Schritte ' num2str(j)])
% 
% figure(2)
% Bx = Basisfunc(x,count,Xm(1:count),sigma,1); 
% dVbias = tau*omega*sum(Bx,2);
% plot(x,dVbias+poth(x)','r','LineWidth',3); hold on
% plot(x,poth(x)','LineWidth',3); hold off
% 
% %plot(x,dVbias-dpoth(x)','r',x,-dpoth(x)','b','LineWidth',3);
% legend('-dPot+Bias','-dPot')
% 
%Integration for Basisfunc
% VbiasP=zeros(1,length(x));
% for k=1:count
%     VbiasP = VbiasP - 1/2*erf((Xm(k)-x)./(sqrt(2)*sigma));
% end
% 

% figure(3)
% sigma=0.3;
% dVbias =@(y) 0.4.* exp(-0.5.*(y+1).^2./sigma.^2);
% plot(x,dVbias(x)'+poth(x)','r','LineWidth',3); hold on
% plot(x,poth(x)','LineWidth',3); hold off
% legend('Pot+Bias','Pot')



% figure(3)
% plot(x,poth(x)-tau*omega*VbiasP);hold on
% plot(x,poth(x));hold off
% title('Pot+Vbias')
% 
% figure(4)
% plot(x,poth(x));hold on 
% plot(x,sqrt(2)*tau*omega*VbiasP);hold off
% 
% %Three well potential
% % x=linspace(-5,5,1000);
% % V=1/200.*(0.5.*x.^6-15.*x.^4+119.*x.^2+28.*x+50);
% % figure(1)
% % plot(x,V)
% x=linspace(-1.5,1.5);
% 
% length_x=length(x);
% Pot_Bias=zeros(1,length_x);
% 
% dx=x(2)-x(1);
% Pot_Bias(1)=omega*sum(Basisfunc(x(1),count,Xm(1:count),sigma,0))*dx;
% figure(5);
% 
% hold on
% xlabel('x')
% title('Biased Potential')
% legend('Original potential','Tilted potential');
% plot(x,poth(x))
% 
%     for j=2:length_x
%         %Compute F(x) from c(x) using the fact that c(x)=-sqrt(2)* dF/dx
%         Pot_Bias(j)=Pot_Bias(j-1)+omega*sum(Basisfunc(x(j),count,Xm(1:count),sigma,0))*dx;
%     end
%     plot(x,+Pot_Bias+poth(x),'r') 
%  hold off

