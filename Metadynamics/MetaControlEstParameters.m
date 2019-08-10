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

Ns=100;
Tl=15000;
dt=1/10000;
Xm=zeros(1,Tl/200);
G=zeros(1,Ns);
G2=zeros(1,Ns);
time=zeros(1,Ns);
ht=zeros(1,Ns);
pt=zeros(1,Ns);
j=0;
k=1;

sigma=zeros(1,Tl/200);
tau=1;
omega=1/10;
count=0;
beta=3;


% eine Trajektorie um das Bias Potential zu berechnen. 

 %for j=1:Tl
 while( (Xtemp > 0.9 && Xtemp < 1.1)==0 )
    
        if j>1 && mod(j,200)== 0
            %if k==200;
            count=count+1;
            Xm(count) = mean(Xhelp);
            sigma(count)= max(abs(Xhelp-Xm(count)))*5;
            k=1;
            Xhelp=zeros(1,200);
            
        end
        
        Bx = Basisfunc(Xtemp,count,Xm(1:count),sigma(1:count),0);
        dVbias = tau*omega*ones(1,count) * Bx';
        
        
        dBt=sqrt(2*dt/beta)*randn;
        Xtemp = Xtemp + (dVbias-dpoth(Xtemp))*dt + dBt;
        Xhelp(k)=Xtemp;

%         if Xtemp > 1.1 && Xtemp <1.2
%             break;
%         end
        
        j=j+1;
        k=k+1;
        
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
    
    for j=1:Tl
    
%         if j>1 && mod(j,200)== 0
%             count=count+1;
%             Xm(count) = Xtemp(j);
%         end
        
        Bx = Basisfunc(Xtemp(j),count,Xm(1:count),sigma,0);
        dVbias = tau*omega*ones(1,count) * Bx';
        
        
        dBt=randn;
        Xtemp(j+1) = Xtemp(j) + (dVbias-dpoth(Xtemp(j)))*dt + sqrt(2*dt/beta)*dBt;
        Gs(j+1)= Gs(j) - (dVbias/sqrt(2/beta))*sqrt(dt)*dBt; 
        Gd(j+1)= Gd(j) - 1/2* (dVbias/sqrt(2/beta))^2*dt; 

         if Xtemp(j+1)>0.9 && Xtemp(j+1)<1.1 && first==0
             first=1;
             ht(i)=j+1;
         end

%          if Xtemp(j+1)>1.1 && Xtemp(j+1)<1.2
%                Gh(count2)=(dVbias/sqrt(2/beta))*sqrt(dt)*dBt- 1/2* (dVbias/sqrt(2/beta))^2*dt;
%                count2=count2+1;
%          end
        
        
    end
    
    if ht(i)== 0
         time(i) = 0;
         G2(i)= 0;
         
    else
        time(i) = j*dt; 
        G2(i)= exp(Gs(ht(i))+Gd(ht(i)));
        
    end
    G(i) = exp(Gs(Tl)+Gd(Tl));
    pt(i)= exp(-beta*ht(i)*dt)*G2(i);
    
 end
%hit 

fprintf('Mittelwert G(N): %2.8f \n',mean(G))
fprintf('Mittelwert G(ht): %2.8f \n',mean(G2(G2 > 0)))
fprintf('E[-beta*tau*Z]: %2.8f \n',mean(pt(pt > 0)))
fprintf('E[-beta*tau*Z]/E[Z]: %2.8f \n',mean(pt(pt > 0))/mean(G2(G2>0)))





% figure(1)
% 
% plot(x,poth(x)), hold on
% plot(Xtemp,poth(Xtemp),'r+'),%hold off
% hold off
% legend('Pot','Pot(X_SDE)')
% title('asymetisches 2 Well Potential')
% xlabel(['Anzahl der Schritte ' num2str(j)])
% 
figure(2)
Bx = Basisfunc(x,count,Xm(1:count),sigma,1); 
dVbias = tau*omega*sum(Bx,2);
plot(x,dVbias-dpoth(x)','r',x,-dpoth(x)','b','LineWidth',3);
legend('Pot+Bias','Pot')
% 
%Integration for Basisfunc
% VbiasP=zeros(1,length(x));
% for k=1:count
%     VbiasP = VbiasP - 1/2*erf((Xm(k)-x)./(sqrt(2)*sigma));
% end
% 
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
