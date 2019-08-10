clear all

%% DRAW POTENTIAL

%close all
% %Range
x=linspace(-0.8,1.8);
%e.g. 
V =@(y) 8.*y.^4-44/3.*y.^3 + 2.*y.^2 +11/3.*y+1;
gradV =@(y) 32.*y.^3-44.*y.^2+4.*y+11/3;

% TRAJECTORY
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


        
 
        x = -0.25;
        t=0;
        while  x < 0.5
            t=t+1;
            eta = randn(1);
            
            
            x = x - gradV(x) * dt + eta * sigma*sdt;
        end

    pt(n) = exp(-beta*(t)*dt);
    time(n) = t*dt;
    end


fprintf('E[exp(-beta*tau)]: %2.8f \n',mean(pt))
fprintf('Var(E[exp(-beta*tau)]): %2.8f \n',var(pt))
fprintf('R(I): %2.8f \n',sqrt(var(pt))/mean(pt))
fprintf('\n')
fprintf('E[time]: %2.8f \n',mean(time))
fprintf('Var[time]: %2.8f \n',var(time))




