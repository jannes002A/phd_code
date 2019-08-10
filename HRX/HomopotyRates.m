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
%t=0.08:-.001:0.02;
t=0.16:-0.04:0;
%homopot = homopoth(x);
% t=0.16;
% homopoth =@(y) poth(y)+96.*t.^2+(4-88.*y+96.*y.^2).*t;
% figure(1)
% plot(x,pot,'k','LineWidth',2)


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

% matrix for saving the trajectories
ma = zeros(N,S);
ra = randn(N,S);
 

num=zeros(length(t),1000);
for i=1:length(t)
    homopoth =@(y) poth(y)+96.*t(i).^2+(4-88.*y+96.*y.^2).*t(i);
    dhomopoth =@(y) dpoth(y)+(-88 +192.*y).*t(i); 
    
%     homopoth =@(y) poth(y)+24.*t(i).^4+(2-44.*y+48.*y.^2).*t(i).^2;
%     dhomopoth =@(y) dpoth(y)+(-44 +48.*y).*t(i).^2; 
    
     Xtemp=ones(1,S)*Xzero;
     ma(1,:)=Xtemp;
     
     
     %j=1;
    for j=2:N
       % while (Xtemp < 0.5)
            Xtemp = Xtemp -(dhomopoth(Xtemp)).*dt + sqrt(dt).*epsilon.*ra(j,:);
            ma(j,:)= Xtemp;
%         if j >=N 
%               break
%         end
        %j=j+1;
 %end
        %j=j+1
%         if ((Xtemp>b) || (Xtemp<a))
%             hit = hit+1;
%             break
%         end
    end
     ma_help=(ma>0.5);

    pos=zeros(1,N);
    for m=1:N
        if (sum(ma_help(:,m))== 0)
            pos(m)=1000;
        else  
          pos(m) = find(ma_help(:,m),1,'first');
        end
    end
    
    for l=1:max(pos)
        num(i,l)= length(find(pos>=l,1000));
    end
end


pnum=mean(num,2);
%  figure()
%  plot(t,pnum,'ko-','LineWidth',2)
%  l1=xlabel('parameter $\lambda$ ');
%  l=ylabel('Average number of steps  that the process stays in the well');
% set(l,'Interpreter','Latex');
%  set(l1,'Interpreter','Latex');

 %%
 %t=[0,0.04,0.08,0.12,0.16];
 figure(1)
 hold on
 for k=1:10:length(t)
     %t(k)
 plot(num(k,:),'LineWidth',4)
 
 xlabel('time')
 ylabel('Number of Processes in Well')
 %title('Prob Exit Rate')
 end
 hold off
 

l= legend('$\lambda = 0.08 $','$\lambda = 0.07 $','$\lambda = 0.06 $','$\lambda = 0.05 $','$\lambda = 0.04 $','$\lambda = 0.03 $','$\lambda = 0.02 $','Location','southwest');

set(l,'Interpreter','Latex');

% lgened('')
 
 
  figure(2)
 %hold on
 for k=1:length(t)

 
 semilogy(num(k,:))
 hold on
 xlabel('time')
 ylabel('Number of Processes in Well')
 title('Prob Exit Rate')
 end 
 %hold off

 
 
 rate=zeros(1,length(t));

 for l=1:length(t)
     
     pos_help=find(num(l,:)> 0,1,'last');
     time=(1:pos_help);
     prob=log(num(l,1:pos_help)/1000);
     p=polyfit(time,prob,1);
     rate(l)=p(1);
   
end
 

 figure(3)
 plot(t,log(-rate),'ko--','LineWidth',2)
 l1=ylabel('log(rate)');
 l=xlabel('parameter $\lambda$');
 set(l,'Interpreter','Latex');
 set(l1,'Interpreter','Latex');
 
 figure(6)
 plot(t,rate)
 
 pos_help=find(-rate>0,1,'last');
 rate_help=log(-rate(1:pos_help));
 
 
 p1=polyfit(t(1:pos_help),rate_help,1);
 p4=polyfit(t(1:pos_help),rate_help,4);
 
%  figure(4)
%  yfit1 = p1(1).*t+p1(2);
%  plot(t,yfit1)
 yfit4 =p4(1).*t.^4+ p4(2).*t.^3+ p4(3).*t.^2+p4(4).*t+p4(5);
 %yfit4 =p4(1).*t.^4+ p4(3).*t.^2+p4(5);
 figure(4)
 plot(t,yfit4);
 hold on
 plot(t,log(-rate));
 hold off
 %semilogy(t,yfit4);
 figure(5)
 semilogy(t,yfit4)
 
 % extrapolation ohne auasreißer
 
 pos_help=find(t>0.03,1,'last');
 rate_help=log(-rate(1:pos_help));
 
 p4=polyfit(t(1:pos_help),rate_help,4);
 %%
 t_help=linspace(0,0.08);
 
 yfit4 =p4(1).*t_help.^4+ p4(2).*t_help.^3+ p4(3).*t_help.^2+p4(4).*t_help+p4(5);
 %yfit4 =p4(1).*t_help.^4+ p4(3).*t_help.^2+p4(5);
 
 figure(7)
 plot(t_help,yfit4,'LineWidth',4);
 hold on
 plot(t(1:5:pos_help),log(-rate(1:5:pos_help)),'o','LineWidth',4);
 hold off
 xlabel('Homotopy Parameter')
 ylabel('log(Rates)')
 legend('Extrapolated Rates','Real Rates','Location','southeast')
 
 
 t=[0,0.04,0.08,0.12,0.16];
 homopoth =@(y,h) poth(y)+96.*h.^2+(4-88.*y+96.*y.^2).*h;
 figure(8)
 
 for i=1:5
     V = homopoth(x,t(i));
     plot(x,V,'LineWidth',4)
     hold on
 end
 xlabel('position space')
 ylabel('potential energy')
 hold off
 ylim([-0.5 17 ])
l=legend('$\lambda = 0 $','$\lambda = 0.04 $','$\lambda = 0.08 $','$\lambda = 0.12 $','$\lambda = 0.16 $','Location','northwest');
 set(l,'Interpreter','Latex');

 
 t=[0,0.04,0.08,0.12,0.16];
 homopoth =@(y,h) poth(y)+96.*h.^2+(4-88.*y+96.*y.^2).*h;
 figure(9)
 for i=1:5
     V = homopoth(x,t(i));
     B = exp(-1 *V);
     plot(x,B,'LineWidth',4)
     hold on
 end
 xlabel('position space')
 ylabel('Boltzman distribution')
 hold off
l=legend('$\lambda = 0 $','$\lambda = 0.04 $','$\lambda = 0.08 $','$\lambda = 0.12 $','$\lambda = 0.16 $','Location','northwest');


%%
% Kurve mit geraden Exponenten
t=0.08:-.001:0.02;
T=zeros(length(t(1:pos_help)),3);
T(:,1)=1;
T(:,2)=t(1:pos_help).^2;
T(:,3)=t(1:pos_help).^4;

par=T\rate_help';
yfit4 =par(3).*t_help.^4+ par(2).*t_help.^2+par(1);

plot(t_help,yfit4,'k','LineWidth',2)
hold on
plot(t(1:5:pos_help),log(-rate(1:5:pos_help)),'ko','LineWidth',4);
hold off
l=xlabel('parameter $\lambda$');
l1=ylabel('log(Rates)');
legend('Extrapolated Rates','Real Rates','Location','southeast')
set(l,'Interpreter','Latex');
set(l1,'Interpreter','Latex');




