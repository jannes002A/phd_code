%Skipt zur Simulation mit Homotopie
clear Xem
load PDETimes.mat
load MCTimes.mat
mct= mtimes;
mcvar = vartimes;
%rng(128,'twister')
rng(7,'twister')

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
h2d=@ (y,t) 96.*y^2-88.*y+4+ 392.*(t.^2/2);


figure(1)
plot(x,homopoth(x,0),'LineWidth',4); hold on
plot(x,homopoth(x,0.15),'LineWidth',4); 
plot(x,homopoth(x,0.2),'LineWidth',4);
plot(x,homopoth(x,0.25),'LineWidth',4); 
hold off
xl=xlabel('$x$');
yl=ylabel('$V(x)$');
l= legend('$\lambda = 0 $','$\lambda = 0.15 $','$\lambda = 0.2 $','$\lambda = 0.25 $','Location','southwest');
set(l,'Interpreter','Latex');
set(xl,'Interpreter','Latex');
set(yl,'Interpreter','Latex');
title('Convoluted Potentials')

figure(2)
plot(x,exp(-0.3*homopoth(x,0))/(sum(exp(-0.3*homopoth(x,0)))),'LineWidth',4); hold on
plot(x,exp(-0.3*homopoth(x,0.15))/(sum(exp(-0.3*homopoth(x,0.15)))),'LineWidth',4); 
plot(x,exp(-0.3*homopoth(x,0.2))/(sum(exp(-0.3*homopoth(x,0.2)))),'LineWidth',4);
plot(x,exp(-0.3*homopoth(x,0.25))/(sum(exp(-0.3*homopoth(x,0.25)))),'LineWidth',4); 
hold off
xl=xlabel('$x$');
yl=ylabel('$\exp(-\beta V(x))$');
l= legend('$\lambda = 0 $','$\lambda = 0.15 $','$\lambda = 0.2 $','$\lambda = 0.25 $','Location','southwest');
set(l,'Interpreter','Latex');
set(xl,'Interpreter','Latex');
set(yl,'Interpreter','Latex');
title('Boltzmann Distributions')


%%
a=0.5;

epsilon = 0.3;
%lambda=[0.01,0.02];
lambda=[0.15,0.2,0.25];


Xzero=-0.25; 
dt=1/1000;
Ns=1000;

pdeVal=1/ext(1);
test=1;
e=zeros(1,test);

for g=1:test
    times=zeros(length(lambda),Ns);
    for i=1:length(lambda)
        for k=1:Ns
            %homopoth =@(y) poth(y)+96.*t(i).^2+(4-88.*y+96.*y.^2).*t(i);
            dhpoth =@(y) dpoth(y)+(-88 +192.*y).*lambda(i)^2/2; 
    
            Xtemp=Xzero;
    
     
            j=1;
            while (Xtemp<=a)
                Xtemp = Xtemp-(dhpoth(Xtemp)).*dt + sqrt(2*dt*epsilon).*randn(1);
                j=j+1;
            end
        times(i,k)=j*dt;
        end
    end

mtimes=mean(times,2);
vartimes=var(times,0,2);
rates= (1./mtimes);

 
par=linspace(eps,lambda(end));
% figure(4)
% semilogy(l,exp(m*l))
% 
%  p=polyfit(log(lambda(1:2))',(rates(1:2)),1);
% p1=p(1);
% p2=p(2);
% %  
% figure(5)
% plot(l,p(2)+p(1).*log(l));hold on
% plot(lambda',rates,'rx')
% % 
% f2=@(x) p(2)+p(1).*log(x);

% p1=[lambda'.^3,lambda'.^2,lambda'.^1,[1 1 1 1]']\log(rates);
% % 
% f3=@(x) exp(p1(1).*x.^3+p1(2).*x.^2+p1(3).*x+p1(4));
% figure(5)
%plot(par,f2(par));hold on
%plot(lambda,rates,'rx');hold off

p=[lambda'.^4,lambda'.^2,[1 1 1]']\log(rates);
rates2=1./(mtimes);
rh= log(rates2) + 1/2*log(-h2d(0.4175,lambda).*h2d(-0.23101,lambda))'-log(2*pi);
p2=1/epsilon.*[lambda'.^4,lambda'.^2,[1 1 1]']\rh;

% Rate as a function of lambda

f4=@(l)  2*pi./(sqrt(-h2d(0.4175,l).*h2d(-0.23101,l))).*exp((p2(1).*l.^4+p2(2).*l.^2+p2(3))/epsilon);
f3=@(l)  exp((p(1).*l.^4+p(2).*l.^2+p(3)));
% figure(5)
% plot(par,f3(par));hold on
% plot(lambda,rates,'rx');hold off

%extrapolated Stopping time
e(g)=1/f3(0);

% f=@(l) rates(1).*exp((l-lambda(1)).*(log(rates(2))-log(rates(1)))./(lambda(2)-lambda(1)));
% f1= rates(1).*exp((0-lambda(1)).*(log(rates(2))-log(rates(1)))./(lambda(2)-lambda(1)));
% figure(6)
% plot(l,f(l));hold on
% plot(lambda,rates,'rx');hold off

%e(g)=abs(f3(0)-pdeVal);

end

1/mean(e) - pdeVal
%2.0835e+03
var(e);
%1.0967e+06
%Konfidenzintervall
xu= mean(e)+1.96*sqrt(var(e)/Ns);
xd= mean(e)-1.96*sqrt(var(e)/Ns);
% 2.2888e+03
% 1.8783e+03

%ExtraRate1 Diff = 9.9834e-05
%ExtraRate2 Diff = 4.9108e-05
%ExtraTime2 Diff = 123.1427
t=0.0:.01:0.25;

% figure(1)
% plot(t(1:end),ext(1:end),'LineWidth',4); hold on
% plot(t,1./f3(t),'LineWidth',4)
% %plot(t,1./f2(t),'LineWidth',4)
% hold off
% xl=xlabel('Lambda');
% yl=ylabel('Average Exit Time for x=-0.25');
% title('Stopping times')
% l= legend('Exact','Extra','Location','northeast');
% set(xl,'Interpreter','Latex');
% set(yl,'Interpreter','Latex');
% set(l,'Interpreter','Latex');
% 
% figure(2)
% semilogy(t(1:end),1./ext,'LineWidth',4);hold on
% semilogy(t(1:end),f3(t(1:end)),'g','LineWidth',4);
% %semilogy(t(1:end),f2(t(1:end)),'k','LineWidth',4);
% semilogy(lambda,rates','rx','LineWidth',4);
% %semilogy(t(1:end),abs(f2(t)),'m');
% hold off
% title('Rates')
% xlabel('Lambda')
% ylabel('Average Exit Rate for x=-0.25')
% l= legend('Exact','Extra','Sampling','Location','northwest');
% set(l,'Interpreter','Latex');
% 
% figure(3)
% plot(t(1:end),1./ext,'LineWidth',4);hold on
% plot(t,f3(t),'LineWidth',4)
% %plot(t,f2(t),'LineWidth',4)
% plot(lambda,rates','rx','LineWidth',4);
% hold off
% title('Rates')
% xl=xlabel('Lambda');
% yl=ylabel('Average Exit Rate for x=-0.25');
% set(xl,'Interpreter','Latex');
% set(yl,'Interpreter','Latex');
% l= legend('Exact','Extra','Sampling','Location','northwest');
% set(l,'Interpreter','Latex');
% 
% 
% figure(4)
% plot(t(1:end),1./(ext),'LineWidth',4);hold on
% plot(t,f4(t),'LineWidth',4)
% %plot(t,f2(t),'LineWidth',4)
% plot(lambda,1./(mtimes),'rx','LineWidth',4);
% hold off
% title('Rates')
% xl=xlabel('Lambda');
% yl=ylabel('Average Exit Rate for x=-0.25');
% set(xl,'Interpreter','Latex');
% set(yl,'Interpreter','Latex');
% l= legend('Exact','Extra','Sampling','Location','northwest');
% set(l,'Interpreter','Latex');

%1/(2*ext(1))-f4(0)

figure(5)
semilogy(t(1:end),1./(ext),'b','LineWidth',4);hold on
semilogy(t(1:end),f4(t(1:end)),'g','LineWidth',4);
%semilogy(t(1:end),f2(t(1:end)),'k','LineWidth',4);
semilogy(lambda,1./(mtimes),'rx','LineWidth',4);
semilogy(0,1./(mct),'mx','LineWidth',4);
%semilogy(t(1:end),abs(f2(t)),'m');
hold off
title('Rates')
xlabel('Lambda')
ylabel('Average Exit Rate for x=-0.25')
l= legend('Exact','Extra','Sampling','MC','Location','northwest');
set(l,'Interpreter','Latex');

figure(6)
plot(t(1:end),ext(1:end),'b','LineWidth',4); hold on
plot(t,1./(f4(t)),'g','LineWidth',4)
plot(lambda,mtimes,'rx','LineWidth',4)
plot(0,mct,'mx','LineWidth',4)
%plot(t,1./f2(t),'LineWidth',4)
hold off
xl=xlabel('Lambda');
yl=ylabel('Average Exit Time for x=-0.25');
title('Stopping times')
l= legend('Exact','Extra','Sampling','MC','Location','northeast');
set(xl,'Interpreter','Latex');
set(yl,'Interpreter','Latex');
set(l,'Interpreter','Latex');
