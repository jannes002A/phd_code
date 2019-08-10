% Cross Entropy Gaussian Process Model Double Well
% 1. Sampling of biased trajectories and evaluation of path functional
% 2. Build Matrix for solving regularized linear equation 


% Sampling 
% a bias can be included here
V=@(x) 1/2.*x.^4-x.^2;
gradV = @(x) 2*x.*(x.^2-1);
dt = 0.01;
sdt = sqrt(dt);
beta = 3;
sigma = sqrt(2/beta);

nvs = 1;
ntrjs = 10; %number of trajectories


nsteps = 150;
time=zeros(1,ntrjs);
X = zeros(ntrjs,nsteps);
X(:,1)=-1;
X_nonbias = zeros(ntrjs,nsteps);
X_nonbias(:,1) = -1;
% the considers path functional is stopping plus some stopping condition
% F(Z) = exp(-1/beta int_0 (tau & T) G(z_s) ds - 1/beta H(z_(tau & T)) )

pathfunc = exp(0.1)*ones(ntrjs,1); % H = -beta*0.1 
bias=0;

for i = 1:ntrjs

        Is=0;
        Id=0;

        x = -1;
        
        for j = 2:nsteps 
            eta=randn(1);
            x = x + (bias - gradV(x) ) * dt + eta * sigma*sdt;
            X(i,j) = x;
            X_nonbias(i,j) = x-bias*dt;
            
            Is = Is - bias * eta/ sigma * sdt;
            Id = Id - bias.^2 / sigma^2 *dt;
            

             if  x > 0.9 && x < 1.1 
             %if x > -1.1 && x < -0.9 && first == 0
                 time(i) = j;
                 pathfunc(i) = exp(-beta*j*dt)*exp(Is+0.5*Id); %weighted path functional 
                 X(i,j:end)= x;
                 X_nonbias(i,j:end) = x-bias*dt;
                
                 break;
             end
        end

end


%%
% Building the matrix to solve the optimization problem 
l=1;
sk=1;
k = @(a,b) sk/sqrt(2*pi)*l.^2.*exp(-0.5*(a-b).^2/l.^2);
%k = @(a,b) 1*exp(-0.5*(a-b).^2/l);


%vector with observed data
% X = X';
% data_obs = X(:);

%K=zeros(length(data_obs),length(data_obs));
K=zeros(nsteps,nsteps);
A=zeros(nsteps,nsteps);
b=zeros(nsteps,1);

for t=1:ntrjs
    for i=1:length(X(t,:))
        for j=1:length(X(t,:))
            K(i,j) = sk/sqrt(2*pi)*l.^2.*exp(-0.5*(X(t,i)-X(t,j)).^2/l.^2);
        end
    end
    A = A + pathfunc(t)*K*dt + 2*beta*eye(nsteps,nsteps); 
    b = b + pathfunc(t)*(K*(X_nonbias(t,:)'));
end
% A= A/ntrjs;
% b= b/ntrjs;
c = A\-b;

figure(1)
plot(c); 
title('Prediction')

% xtest = linspace(-2,2,150);
% 
% figure(2)
% plot(xtest, -gradV(xtest)); hold on
% plot(xtest, (c - gradV(xtest)'))
% title('Gradients')
% hold off
% % 
% % 
% dx=xtest(2)-xtest(1);
% per_pot = zeros(1,nsteps+1);
%  
% for i=2:nsteps+1
%     per_pot(i) = per_pot(i-1) + ( - c(i-1) + gradV(xtest(i-1)) )*dx;
% end
% 
% figure(3)
% plot(per_pot(2:end)+4); hold on
% plot(V(xtest)); hold off
% title('Perturbed Potential')


% c Prediction
n_pred=100;
x_pred = linspace(-2,2,n_pred);
K_pred = zeros(length(x_pred),nsteps);
c_pred=0;

for t=1:ntrjs
    for i=1:length(x_pred)
        for j=1:length(X(t,:))
            K_pred(i,j) = sk/sqrt(2*pi)*l.^2.*exp(-0.5*(x_pred(i)-X(t,j)).^2/l.^2);
        end
    end
    c_pred =  c_pred - pathfunc(t)/(2*beta*ntrjs)*( K_pred*c*dt + K_pred*(X_nonbias(t,:)'));  
end


figure(4)
plot(x_pred, c_pred)


figure(5)
plot(x_pred,  -gradV(x_pred)); hold on
plot(x_pred,   c_pred -gradV(x_pred)'  )
legend('-gradV','-gradV+cPred')
title('Gradients')
hold off

dx=x_pred(2)-x_pred(1);
per_pot = zeros(1,n_pred+1);
 
for i=2:n_pred+1
    per_pot(i) = per_pot(i-1) + (-c_pred(i-1) + gradV(x_pred(i-1)) )*dx;
end

figure(6)
plot(x_pred,V(x_pred)); hold on
plot(x_pred,per_pot(2:end)+4); hold off
title('Perturbed Potential')