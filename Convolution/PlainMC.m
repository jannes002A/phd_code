%Skipt zur Simulation Monte Carlo

%close all
% %Range
x=linspace(-0.5,1.5);
%e.g. 
poth =@(y) 8.*y.^4-44/3.*y.^3 + 2.*y.^2 +11/3.*y+1;
dpoth =@(y) 32.*y.^3-44.*y.^2+4.*y+11/3;
pot = poth(x);
dpot= dpoth(x);

epsilon = 0.3;

N=10^6; Xzero=-0.25; dt=1/1000;
a=0.5;
Ns=1000;

cv=1;

times=zeros(1,Ns);


 for i=1:cv
    Xtemp = Xzero;
    
           for k=1:Ns    
            Xtemp=Xzero;
            fprintf('Sampling %3i \n', k )
     
            j=1;
                while (Xtemp<=a)
                    Xtemp = Xtemp-(dpoth(Xtemp)).*dt + sqrt(2*dt.*epsilon).*randn(1);
                    j=j+1;
                end
                times(i,k)=j*dt;
           end
           
 end

mcmtimes=mean(times);
mcvartimes=var(times);
%  mtimes = 1.6958e+03
%  1/mtimes = 5.8968e-04 
%  vtimes = 2.9471e+06

save('MCTimes.mat','mt','mtimes','vart','vartimes','-v7.3')

xu= mean(times)+1.96*sqrt(var(times)/Ns);
xd= mean(times)-1.96*sqrt(var(times)/Ns);

% mrate=zeros(length(t),1);
% cc=zeros(length(t),1);
% sigma=zeros(length(t),1);
% for i=1:length(t)
%     mrate(i) = mean(Xtime(i,:));
%     sigma(i)= sum((Xtime(i,:)-mrate(i)).^2);
%     cc(length(t)+1-i) = 1.96*sqrt(sigma(i)/S);
% end
% 
% %pfit=p(1).*t.^4+p(2).*t.^3+p(3).*t+p(4);
% %pfit=p(1).*t.^2+p(2).*t+p(3);
% l=linspace(0,0.16);
% pfit=p(1).*t+p(2);
% I=[ones(1,length(hr)); th.^2 ;th.^4]';
% ex=I\rates(hr);
% fit=ex(3).*l.^4+ex(2).*l.^2+ex(1);
% figure(3)
% plot(l,fit,'Linewidth',2)
% hold on
% plot(t(hr),1./mtimes(hr),'.-r','Linewidt',2)
% %plot(t,1./ext,'g','Linewidth',2)
% hold off
% legend('extr. rate','samp. rate','Location','NorthWest')
% xlabel('Lambda')
% ylabel('Exit rate for x=-.25')
% % figure(1)
% % plot(t,pfit,'g');
% 
% figure(1)
% errorbar(t(hr),rates(hr),cc(hr),'r.','Linewidth',2)
% xlabel('Lambda')
%ylabel('Variance of the estimator')

