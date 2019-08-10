% Cross Entropy Gaussian Process Model Double Well
% Corrected Nonparametric Estimator with different Girsanov

clear 
% Sampling 
% a bias can be included here
V=@(x)  1/2.*(x.^2-1).^2;
gradV = @(x) 2*x.*(x.^2-1);
dt = 0.001;
sdt = sqrt(dt);
beta = 4;
sigma = sqrt(2/beta);

nvs = 1;
ntrjs = 100; %number of trajectories
opt_steps=2;
nsteps = 150;
npoints = 100;
points = linspace(-2,2,npoints);
n_pred=100;
sk=1;
l=.1;
pathfunc =  ones(ntrjs,1);
girsanov =  ones(ntrjs,1);
c_old = zeros(opt_steps,n_pred);
time = zeros(ntrjs,1);
b = zeros(npoints,1);

A = zeros(npoints,npoints);


% the considers path functional is the moment generating function of the
% stopping time

bias=0;

for opt = 1:opt_steps

        
        %Eta= randn(ntrjs,nsteps-1);
        time=zeros(1,ntrjs);
        costfunc = zeros(1,ntrjs);
        trajlength = ones(ntrjs+1,1);
        
        
    for i = 1:ntrjs
        
        Is = 0;
        Id = 0;
        x = -1;
        t=0;
        cost=0;
        j=1;
        X = [];
        X(1)=-1;
        B = [];
        B(1)=0;
        Ki = zeros(npoints,npoints);
        Kii= zeros(npoints,1);
        Kiii = zeros(npoints,npoints);
        
        while (x<0) 
            j=j+1;
            t=t+1;
            eta=randn(1);
            
            
            if opt==1
                x = x - (gradV(x)+bias) * dt + eta * sigma*sdt;
            else
                K_pred =  sk*exp(-0.5*(x-points).^2/l.^2);
                bias =  K_pred*par; 
                x = x - (gradV(x)+bias)*dt + eta * sigma*sdt;
            end
            X(j) = x;
            B(j)= eta*sdt;
            cost=cost+bias*bias*dt;
            
            Is = Is + bias * eta/ sigma * sdt;
            Id = Id + (bias).^2 / sigma^2 *dt;

 
        end
        time(i)=t;
        girsanov(i)=exp(Is-0.5*Id);
        pathfunc(i)=exp(-1/beta*t*dt)*exp(Is-0.5*Id);
        costfunc(i)=1/4 * cost;
        
        for k=1:npoints
            for l=1:npoints
                Ki(k,l)= Ki(k,l) + sk.*exp(-0.5*(X-points(k)).^2/l.^2)*(sk.*exp(-0.5*(X-points(l)).^2/l.^2)')*dt;
            end
        end
        A=A+pathfunc(i)/sigma^2*Ki;
        
        
        
        for k=1:length(X)
            Kii =  Kii + sk.*exp(-0.5*(X(k)-points').^2/l.^2)*B(k);
        end
        
        b = b+ pathfunc(i).*(1/sigma).*(Kii);
        

    end
    
 
    
    
    for k=1:npoints
        Kiii(k,:) = sk.*exp(-0.5*(points(k)-points).^2/l.^2);
    end
    A=(A./ntrjs)+0.1*Kiii+0.001*speye(npoints,npoints);
    b=b./ntrjs;

    
    par = A\-b;
    
    n_pred=100;
    x_pred=linspace(-2,2,n_pred);
    K_pred=zeros(1,npoints);
    for k=1:n_pred
        K_pred(k) =  sk.*exp(-0.5*(x_pred(k)-points).^2/l.^2)*par;
    end
    
    figure(opt)
    plot(x_pred,-K_pred);
    
    
    p=mean(time*dt);
    fprintf('Mean average hitting time %f \n', p )
    fprintf('Var(time) %f \n', var(time*dt))
    fprintf('Mean average pathfunctional	 %f \n', mean(pathfunc) )
    fprintf('Variance pathfunctional %f \n', var(pathfunc))
    fprintf('Mean average Girsanov	 %f \n', mean(girsanov) )
    fprintf('Variance Girsanov %f \n', var(girsanov))
    
end
%%
% figure(1)
% plot(x_pred, c_old(end,:),'LineWidth',3)
% title('Control')

figure(6)
plot(x_pred,  -gradV(x_pred),x_pred, -(gradV(x_pred)'+ K_pred),'LineWidth',3 ); 
legend('-gradV','-gradV+cPred')
title('Gradients')

dx=x_pred(2)-x_pred(1);
per_pot = zeros(1,n_pred+1);
control = zeros(1,n_pred+1); 
c_pred = K_pred;
 
for i=2:n_pred+1
    per_pot(i) = per_pot(i-1) + (c_pred(i-1) + gradV(x_pred(i-1)) )*dx;
    control(i) = control(i-1) + c_pred(i-1)*dx;
end

% figure(7)
% plot(x_pred,V(x_pred),x_pred,per_pot(2:end)+14.5,'LineWidth',3); 
% legend('Potential','Perturbed Potential')
% title('Perturbed Potential')

figure(8)
plot(x_pred,V(x_pred),x_pred,V(x_pred)+control(2:end)+0,'LineWidth',3); 
legend('Potential','Perturbed Potential')

% figure(9)
% plot(x_pred,control(2:end)+11,'LineWidth',3); 
% legend('Predicted Control (Integral)')

