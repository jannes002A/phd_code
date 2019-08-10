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
% For the other Girsanov Formula we use exp as the gradient 
% We need the integral and the 2nd derivative of the ansatz function
% d/dx = w/(sqrt(2 pi)l) exp(-1/2*(x-c)^2/(2l^2))
% int = -1/2 w erf((c-x)/(sqrt(2)l))
% d^2/dx^2 = w(x-c)/(sqrt(2 pi)l^3) exp(-(x-c)^2/(2l^2))

poth =@(y)  1/2.*(y.^2-1).^2;
dpoth =@(y) 2.*y.*(y.^2-1);

pot = poth(x);
dpot= dpoth(x);
Xzero=-1;





Xtemp=Xzero;

Ns=1000;
Tl=15000;
dt=0.0001;
Xm=zeros(1,1000);
G=zeros(1,Ns);
G2=zeros(1,Ns);
time=ones(1,Ns)*Tl*dt;
ht=zeros(1,Ns);
pt=zeros(1,Ns);
jh=0;

sD=0.8*ones(1,1000);
omega=1/10;
count=0;
beta=3;

dv2=@(y,c) omega.*(y-c)./(sqrt(2*pi).*sD(1).^3).*exp(-(y-c).^2./(2.*sD(1).^2));

% eine Trajektorie um das Bias Potential zu berechnen. 

 %for j=1:Tl
 while( (Xtemp > 0.9 && Xtemp < 1.1)==0 )
    
        if jh>1 && mod(jh,200)== 0
            count=count+1;
            Xm(count) = Xtemp;
        end
        if count >= 1000
            error('Trajektorie zu lang!')
        end
        
        Bx = Basisfunc(Xtemp,count,Xm(1:count),sD(1:count),0);
        dVbias = omega*ones(1,count) * Bx';
        
        
        dBt=sqrt(2*dt/beta)*randn;
        Xtemp = Xtemp + (dVbias-dpoth(Xtemp))*dt + dBt;


%         if Xtemp > 1.1 && Xtemp <1.2
%             break;
%         end
        
        jh=jh+1;
        
 end
 fprintf('Steps needed in the first run %8.2f \n',jh)
% Sampling mit dem fixen Bias Potential   
Xhm=Xm(1:count);
G0 = 1/(2*beta^-1)*sum(-1/2*omega*erf((Xhm-Xzero)./(sqrt(2)*sD(1))));

parfor i=1:Ns
     
     Xtemp=zeros(1,Tl+1);
     Gs=zeros(1,Tl+1);
     Gd=zeros(1,Tl+1);
     Gh=zeros(1,Tl+1);
     Xtemp(1)=Xzero;
     first=0;
     Ghit=0;
     
    for j=1:Tl
    
%         if j>1 && mod(j,200)== 0
%             count=count+1;
%             Xm(count) = Xtemp(j);
%         end
        
        Bx = Basisfunc(Xtemp(j),count,Xhm,sD,0);
        dVbias = omega*ones(1,count) * Bx';
        
        
        dBt=randn;
        dpot=dpoth(Xtemp(j+1));
        Xtemp(j+1) = Xtemp(j) + (dVbias-dpot)*dt + sqrt(2*dt/beta)*dBt;
        Gs(j+1)= Gs(j) + 1/(2*beta^-1).* ( -dVbias*dpot  + 1/2.*dVbias.*dVbias + (beta^-1)*sum(dv2(Xtemp(j),Xhm)) )*dt; 
        

         if Xtemp(j+1)>0.9 && Xtemp(j+1)<1.1 && first==0
             first=1;
             ht(i)=j+1;
             Ghit = -1/(2*beta^-1)*sum(-1/2*omega*erf((Xhm-Xtemp(j))./(sqrt(2)*sD(1))));
         end

%          if Xtemp(j+1)>1.1 && Xtemp(j+1)<1.2
%                Gh(count2)=(dVbias/sqrt(2/beta))*sqrt(dt)*dBt- 1/2* (dVbias/sqrt(2/beta))^2*dt;
%                count2=count2+1;
%          end
        
        
    end
    
    if ht(i)== 0
         G2(i)= 0;
         
    else
        time(i) = ht(i)*dt; 
        G2(i)= exp(Ghit+Gs(ht(i))+G0);
        
    end
    Glast = -1/(2*beta^-1)*sum(-1/2*omega*erf((Xhm-Xtemp(j))./(sqrt(2)*sD(1))));
    G(i) = exp(Gs(j)+Glast+G0);
    pt(i)= exp(-beta*ht(i)*dt)*G2(i);
    
 end
%hit 

fprintf('Girsanov Reweighting fixed Parameter and different Girsanov \n')
fprintf('Mittelwert G(N): %2.8f \n',mean(G))
fprintf('Mittelwert G(ht): %2.8f \n',mean(G2(G2 > 0)))
fprintf('Var(E[-beta*tau*Z]): %2.8f \n',var(G2(G2 > 0)))
fprintf('R(I): %2.8f \n',sqrt(var(G2(G2 > 0)))/mean(G2(G2 > 0)))

fprintf('E[-beta*tau*Z]: %2.8f \n',mean(pt(pt > 0)))
fprintf('Var(E[-beta*tau*Z]): %2.8f \n',var(pt(pt > 0)))
fprintf('R(I): %2.8f \n',sqrt(var(pt(pt > 0)))/mean(pt(pt > 0)))

fprintf('E[time]: %2.8f \n',mean(time))





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
Bx = Basisfunc(x,count,Xm(1:count),sD,1); 
dVbias = omega*sum(Bx,2);
plot(x,dVbias-dpoth(x)','r',x,-dpoth(x)','b',x,dVbias,'LineWidth',3);
legend('Pot+Bias','Pot')
% 
%Integration for Basisfunc
% VbiasP=zeros(1,length(x));
% for k=1:count
%     VbiasP = VbiasP - 1/2*erf((Xm(k)-x)./(sqrt(2)*sD));
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
