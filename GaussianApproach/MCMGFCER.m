% Cross Entropy Gaussian Process Model Double Well
% Corrected Nonparametric Estimator with different Girsanov

clear 
% Sampling 
% a bias can be included here
V=@(x)  1/2.*(x.^2-1).^2;
gradV = @(x) 2*x.*(x.^2-1);
dt = 0.001;
sdt = sqrt(dt);
beta = 2;
sigma = sqrt(2/beta);

nvs = 1;
ntrjs = 10000; %number of trajectories



% the considers path functional is the moment generating function of the
% stopping time



        
        %Eta= randn(ntrjs,nsteps-1);
time=zeros(1,ntrjs);
pathfunc = zeros(1,ntrjs);

        
    for i = 1:ntrjs
        
        t=0;
        x = -1;



        
        while (x<0) 
            t=t+1;
            eta=randn(1);
            x = x - (gradV(x)) * dt + eta * sigma*sdt;
         

 
        end
        time(i)=t;
        pathfunc(i)=exp(-1/beta*t*dt);
    end
    


  

    
    
    p=mean(time*dt);
    fprintf('Mean average hitting time %f \n', p )
    fprintf('Var(time) %f \n', var(time*dt))
    fprintf('Mean average pathfunctional	 %f \n', mean(pathfunc) )
    fprintf('Variance pathfunctional %f \n', var(pathfunc))

    
