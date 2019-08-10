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
ntrjs = 10000; %number of trajectories
finalp = zeros(nvs,1);
finalpw = zeros(nvs,1);
finalpt = zeros(nvs,1);


nsteps = 15000;





for v=1:nvs
    %rng(v,'twister')

    %Girsanov weight
    W = zeros(1,ntrjs);

    %Estimators
    p = zeros(1,ntrjs);    %P(X_t in T|t<T)
    pt = zeros(1,ntrjs);   %E[exp(-beta*tau)Z]
    time = nsteps*ones(1,ntrjs); %time of each trajectory

    for n = 1:ntrjs
%         if mod(n,10)==0
%             %clc
%             disp([num2str(n/ntrjs*100), '%'])
%         end

        %Girsanov integrals
        eta = randn(nsteps,1);
        I1 = 0;
        I2 = 0;
        X = zeros(nsteps,1);

        x = -0.25;
        first = 0;
        for t = 1:nsteps-1 
            
            hgrad = - (gradVc(x,l));
            
            x = x + (hgrad) * dt + eta(t) * sigma*sdt;
            X(t) = x;
            
            diffDrift = gradV(x)+hgrad; 
            
            I1 = I1 - (diffDrift) * eta(t) / sigma * sdt; %- nabla U(x)
            I2 = I2 - (diffDrift).^2 / sigma^2 *dt;

            if  x > 1.15 && x < 1.25 && first == 0
            %if x > -1.1 && x < -0.9 && first == 0

                p(n) = 1;
                pt(n) = exp(-beta*(t)*dt);
                time(n) = t - 1;
                first = 1;
                break;
            end
        end


        if p(n)== 0
            W(n) = 0;         
        else
            W(n) = exp( I1 + 0.5*I2); 
        end
        
        pt(n) = pt(n) * W(n);
    end
    finalp(v) = mean(p);

    finalpw(v) = mean(W);
    finalpt(v) = sum(pt(pt>0))/ntrjs;
    if isnan(finalpt(v))==1
         finalpt(v) = 0;
    end
end

fprintf('E[P(X_t in B| t< T)]: %2.8f \n',mean(finalpw))
fprintf('Var[P(X_t in B| t< T)]: %2.8f \n',var(W))
fprintf('R(I): %2.8f \n',sqrt(var(W))/mean(finalpw))
fprintf('\n')
fprintf('E[exp(-beta*tau)]: %2.8f \n',mean(finalpt))
fprintf('Var(E[exp(-beta*tau)]): %2.8f \n',var(pt))
fprintf('R(I): %2.8f \n',sqrt(var(pt))/mean(finalpt))
fprintf('\n')
fprintf('Trajectories in T: %2i \n', sum(time<15000)/ntrjs)

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

