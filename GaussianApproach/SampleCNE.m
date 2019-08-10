% Cross Entropy Gaussian Process Model Double Well
% 1. Sampling of biased trajectories and evaluation of path functional
% 2. Build Matrix for solving regularized linear equation 

clear 
% Sampling 
% a bias can be included here
V=@(x)  1/2.*(x.^2-1).^2;
gradV = @(x) 2*x.*(x.^2-1);
dt = 0.01;
sdt = sqrt(dt);
beta = 3;
sigma = sqrt(2/beta);

nvs = 1;
ntrjs = 100; %number of trajectories
strjs = 1000;
opt_steps=1;
nsteps = 150;
n_pred=100;
sk=8;
l=0.2;
pathfunc =  ones(ntrjs,1);
spathfunc = ones(strjs,1);
stime = zeros(strjs,1);

c_old = zeros(opt_steps,n_pred);
c_pred=0;
girsanov = zeros(1,ntrjs);


% the considers path functional is the moment generating function of the
% stopping time
        
        bias=0;
        X = zeros(ntrjs,nsteps);
        X(:,1)=-1;
        X_nonbias = zeros(ntrjs,nsteps);
        X_nonbias(:,1) = -1;

    for i = 1:ntrjs
       x = -1;
        
        for j = 2:nsteps
            eta=randn(1);
            
            x = x + (- gradV(x)) * dt + eta * sigma*sdt;
            
            X(i,j)=x;
            X_nonbias(i,j) = x;

             if  x > 0.9 && x < 1.1 
                 
                 pathfunc(i) = 1; %weighted path functional 
                 
                 break;
             else
                 pathfunc(i) = 0;
             end
        end

    end
    
    data = X';
    data= data(:);
    data_nonbias = X_nonbias';
    data_nonbias = data_nonbias(:);
    %pathweight = repmat(pathfunc./ntrjs,nsteps,1);
    
    
    K=zeros(nsteps*ntrjs,nsteps*ntrjs);

    pathweight = repmat(pathfunc,nsteps,1);
    pathweight= pathweight(:);
    
   
    for j=1:length(data)
            K(j,:) = (pathweight'/ntrjs).*sk/sqrt(2*pi)*l.^2.*exp(-0.5*(data'-data(j)).^2/l.^2);
    end
    
    A = K*dt + 2/beta*speye(nsteps*ntrjs,nsteps*ntrjs); 
    b = K*data_nonbias;
    
    c = A\-b;
    
    pathweight = pathweight(:)./((2/beta)*(ntrjs));
    
    for i = 1:strjs

        Is=0;
        Id=0;

        x = -1;
        
        for j = 2:nsteps
            eta=randn(1);
           
                
            K_pred = pathweight'.* sk/sqrt(2*pi)*l.^2.*exp(-0.5*(data'-x).^2/l.^2);
            bias = - K_pred*c*dt - K_pred*(data_nonbias);  
            x = x + (-gradV(x)+bias)*dt + eta * sigma*sdt;
            
            Is = Is - bias * eta/ sigma * sdt;
            Id = Id - bias.^2 / sigma^2 *dt;
            

             if  x > 0.9 && x < 1.1 
                 stime(i) = j;
                 spathfunc(i) = exp(Is+0.5*Id); %weighted path functional 
                 girsanov(i) = exp(Is+0.5*Id);
                 break;
             else
                 spathfunc(i) = 0;
             end
        end

    end
    p=[sum(stime>0),strjs];
    fprintf('Trajectories in T %d / %d \n', p )
   
    
    % Estimators
    mgf = spathfunc;
    
    


fprintf('E[P(X_t in B| t< T)]: %2.8f \n',mean(girsanov))
fprintf('Var[P(X_t in B| t< T)]: %2.8f \n',var(girsanov))
fprintf('R(I): %2.8f \n',sqrt(var(girsanov))/mean(girsanov))
fprintf('\n')
fprintf('E[exp(-beta*tau)]: %2.8f \n',mean(mgf))
fprintf('Var(E[exp(-beta*tau)]): %2.8f \n',var(mgf))
fprintf('R(I): %2.8f \n',sqrt(var(mgf))/mean(mgf))


