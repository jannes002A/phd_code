%Skipt zur Simulation mit Homotopie
clear Xem

rng(3,'twister')                 % restore the previous settings


%close all
% %Range
x=linspace(-0.5,1.5);
%e.g. 
poth =@(y) 8.*y.^4-44/3.*y.^3 + 2.*y.^2 +11/3.*y+1;
dpoth =@(y) 32.*y.^3-44.*y.^2+4.*y+11/3;
%homopoth =@(y,t) poth(y)+96.*t.^2+(4-88.*y+96.*y.^2).*t;
homopoth =@(y,t) poth(y)+96.*(t^2/2).^2+(4-88.*y+96.*y.^2).*(t^2/2);

pot = poth(x);
dpot= dpoth(x);


% figure(1)
% plot(x,homopoth(x,0),'LineWidth',4); hold on
% plot(x,homopoth(x,0.15),'LineWidth',4); 
% plot(x,homopoth(x,0.2),'LineWidth',4);
% plot(x,homopoth(x,0.25),'LineWidth',4); 
% hold off
% xl=xlabel('$x$');
% yl=ylabel('$V(x)$');
% l= legend('$\lambda = 0 $','$\lambda = 0.15 $','$\lambda = 0.2 $','$\lambda = 0.25 $','Location','southwest');
% set(l,'Interpreter','Latex');
% set(xl,'Interpreter','Latex');
% set(yl,'Interpreter','Latex');


a=0.5;

beta = 0.5;
lambda=0.15;


Xzero=-0.25; 
dt=1/1000;
Ns=100;

times=zeros(1,Ns);

S=zeros(1,Ns);



test=1;
e=zeros(1,test);

nMax=10000;
G=zeros(length(lambda),Ns);

dhomopoth =@(y) dpoth(y)+(-88 +192.*y).*lambda^2/2; 
diffdpoth=@(y) (-88 +192.*y).*lambda^2/2;

    
 for k=1:Ns

            
            
            Xtemp=zeros(1,nMax);
            Gs=zeros(1,Ns);
            Gd=zeros(1,Ns);
            Xtemp(1)=Xzero;
            
            for j=1:nMax
                dBt=sqrt(2*dt).*randn(1);
                Xtemp(j+1) = Xtemp(j) -(dhomopoth(Xtemp(j))).*dt + beta.*dBt;
                Gs(j+1) = Gs(j) - diffdpoth(Xtemp(j+1))/beta .*dBt; 
                Gd(j+1) = Gd(j) - 1/2*(diffdpoth(Xtemp(j+1))/beta).^2.*dt;
                %SI = SI + dBt;
                
                if Xtemp(j+1)> a
                    break;
                end
                
                
            end
        times(k)=j*dt;
        G(k) = exp(Gs(j+1)+Gd(j+1));
        %S(i,k)=SI;
end
   

fprintf('Mittelwert: %2.16f \n',mean(G))
fprintf('Varianz: %2.16f \n',var(G))
fprintf('rel. Error: %2.16f\n',std(G)/mean(G))

%gtimes=times.*(G-repmat(mean(G,2),1,100));
%gtimes2=times.*G;
%mgt=mean(gtimes2,2);
%vartimes=var(gtimes2,0,2);
 
% par=linspace(eps,lambda(end));
% % figure(4)
% % semilogy(l,exp(m*l))
% % 
% %  p=polyfit(log(lambda(1:2))',(rates(1:2)),1);
% % p1=p(1);
% % p2=p(2);
% % %  
% % figure(5)
% % plot(l,p(2)+p(1).*log(l));hold on
% % plot(lambda',rates,'rx')
% % % 
% % f2=@(x) p(2)+p(1).*log(x);
% 
% % p=[lambda'.^2,[1 1 1]']\log(rates);
% % % 
% % f2=@(x) exp(p(1).*x.^2+p(2));
% % figure(5)
% % plot(par,f2(par));hold on
% % plot(lambda,rates,'rx');hold off
% 
% p=[lambda'.^4,lambda'.^2,[1 1 1]']\log(rates);
% % % 
% f3=@(x) exp(p(1).*x.^4+p(2).*x.^2+p(3));
% figure(5)
% plot(par,f3(par));hold on
% plot(lambda,rates,'rx');hold off
% 
% 
% f=@(l) rates(1).*exp((l-lambda(1)).*(log(rates(2))-log(rates(1)))./(lambda(2)-lambda(1)));
% f1= rates(1).*exp((0-lambda(1)).*(log(rates(2))-log(rates(1)))./(lambda(2)-lambda(1)));
% figure(6)
% plot(l,f(l));hold on
% plot(lambda,rates,'rx');hold off

%e(g)=abs(f3(0)-pdeVal);



