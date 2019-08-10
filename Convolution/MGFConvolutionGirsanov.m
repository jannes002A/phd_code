% Convolution Girsanov reweighting
% Moment Generating function


clear all

%% DRAW POTENTIAL

%close all
% %Range
x=linspace(-0.8,1.8);
%e.g. 
V =@(y) 8.*y.^4-44/3.*y.^3 + 2.*y.^2 +11/3.*y+1;
gradV =@(y) 32.*y.^3-44.*y.^2+4.*y+11/3;
%homopoth =@(y,t) poth(y)+96.*t.^2+(4-88.*y+96.*y.^2).*t;
Vc =@(y,t) V(y)+96.*(t^2/2).^2+(4-88.*y+96.*y.^2).*(t^2/2);
gradVc =@(y,t) gradV(y)+(-88 +192.*y).*t^2/2; 

%convolution parameter
l=0.2;

figure(1)
plot(x,V(x),'LineWidth',3); hold on
plot(x,Vc(x,l),'LineWidth',3); hold off

figure(2)
plot(x,gradV(x),'LineWidth',3); hold on
plot(x,gradVc(x,l),'LineWidth',3); hold off



%% TRAJECTORY
dt = 0.0001;
sdt = sqrt(dt);
beta = 3;
sigma = sqrt( 2/beta);

nvs = 1;
ntrjs = 100; %number of trajectories








for v=1:nvs
    %rng(v,'twister')



    %Estimators
 
    time = zeros(ntrjs,1);
    mgf = zeros(ntrjs,1);

    for n = 1:ntrjs

        
        I1 = 0;
        I2 = 0;
        
        x = -0.25;
        t=0;
        while (x < 0.41)
            t=t+1;
            eta=randn(1);
            hgrad = - (gradVc(x,l));
            
            x = x + (hgrad) * dt + eta * sigma*sdt;
            
            
            diffDrift = gradV(x) +hgrad; 
            
            I1 = I1 - (diffDrift) * eta / sigma * sdt; %- nabla U(x)
            I2 = I2 - (diffDrift).^2 / sigma^2 *dt;

            
        end
        time(n)=t*dt;
        mgf(n) = exp(-beta*t)*exp(I1+0.5*I2);
    end

end

fprintf('E[MGF] %2.8f \n',mean(mgf))
fprintf('Var[MGF]: %2.8f \n',var(mgf))
fprintf('E[time] %2.8f \n',mean(time))
fprintf('Var[time] %2.8f \n',var(time))


%%
%Plots
% x = linspace(-1.6,1.6,1000);
% V = @(x) 0.5*(x.^2-1).^2;
% 
% 
% gradV = @(x) 2*x.*(x.^2-1);
% Vx = V(x);
% gV= gradV(x);
% dvbias=zeros(1,length(x));
% vbias=zeros(1,length(x));
% 
% for i=1:length(x0)
%     dvbias = dvbias -1/s^2*(x - x0(i)).*w .* exp( - 0.5.*(x - x0(i)).^2/s^2);
%     vbias = vbias + w .* exp( - 0.5.*(x - x0(i)).^2/s^2);
% end
% 
% 
% figure(1)
% plot(x,dvbias,'g',x,-(gV+dvbias),'r',x,-gV,'b','LineWidth',3)
% 
% figure(2)
% plot(x,vbias,'g',x,Vx+vbias,'r',x,Vx,'b','LineWidth',3)
% 
% figure(3)
% plot(x,V_meta,'g',x,Vx+V_meta,'r',x,Vx,'b','LineWidth',3)
% legend('Bias','Pot+Bias','Pot')
% 
% figure(4)
% plot(x,V_meta,'k',x,Vx+V_meta,':k',x,Vx,'-.k','LineWidth',3)
% legend('Bias','Pot+Bias','Pot')

