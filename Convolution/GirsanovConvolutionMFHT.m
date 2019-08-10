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
legend('Org. Pot', 'Conv. Pot')


figure(2)
plot(x,gradV(x),'LineWidth',3); hold on
plot(x,gradVc(x,l),'LineWidth',3); hold off



%% TRAJECTORY
dt = 0.0001;
sdt = sqrt(dt);
beta = 3;
sigma = sqrt( 2/beta);
ntrjs = 1000; %number of trajectories




    %rng(v,'twister')



    %Estimators

    pt = zeros(1,ntrjs);   %E[exp(-beta*tau)Z]
    time = zeros(1,ntrjs); %time of each trajectory

    for n = 1:ntrjs

      %Girsanov integrals
        
        I1 = 0;
        I2 = 0;
        x = -0.25;
        t=0;
        while  x < 0.5
            t=t+1;
            eta = randn(1);
            hgrad = - (gradVc(x,l));
            
            x = x + (hgrad) * dt + eta * sigma*sdt;
      
            
            diffDrift = gradV(x)+hgrad; 
            
            I1 = I1 - (diffDrift) * eta / sigma * sdt; %- nabla U(x)
            I2 = I2 - (diffDrift).^2 / sigma^2 *dt;
        end

    pt(n) = exp(-beta*(t)*dt)*exp(I1+0.5*I2);
    time(n) = t*dt;
    end


fprintf('E[exp(-beta*tau)]: %2.8f \n',mean(pt))
fprintf('Var(E[exp(-beta*tau)]): %2.8f \n',var(pt))
fprintf('R(I): %2.8f \n',sqrt(var(pt))/mean(pt))
fprintf('\n')
fprintf('E[time]: %2.8f \n',mean(time))
fprintf('Var[time]: %2.8f \n',var(time))




