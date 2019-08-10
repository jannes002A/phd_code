% Cross Entropy Gaussian Process Model Double Well
% Corrected Nonparametric Estimator with different Girsanov

clear 
% Sampling 
% a bias can be included here
V=@(x) 1/2.*(x.^2-1).^2;
gradV = @(x) 2*x.*(x.^2-1);
dt = 0.01;
sdt = sqrt(dt);
beta = 3;
sigma = sqrt(2/beta);

nvs = 1;
ntrjs = 50; %number of trajectories
opt_steps=2;
nsteps = 150;
n_pred=100;
sk=20;
l=0.1;
pathfunc =  ones(ntrjs,1);
c_old = zeros(opt_steps,n_pred);
c_pred=0;

% the considers path functional is the moment generating function of the
% stopping time

bias=0;

for opt = 1:opt_steps

        
        %Eta= randn(ntrjs,nsteps-1);
        time=zeros(1,ntrjs);
        X = zeros(ntrjs,nsteps);
        X(:,1)=1;
        X_nonbias = zeros(ntrjs,nsteps);
        X_nonbias(:,1) = 1;
        pathweight = repmat(pathfunc,nsteps,1);
        pathweight = pathweight(:)./((2/beta)*(ntrjs));
    for i = 1:ntrjs

        Is=0;
        Id=0;

        x = 1;
        
        for j = 2:nsteps
            eta=randn(1);
            
            
            if opt==1
                x = x + (- gradV(x)+bias) * dt + eta * sigma*sdt;
            else
                
                K_pred = pathweight'.* sk/sqrt(2*pi)*l.^2.*exp(-0.5*(data'-x).^2/l.^2);
                bias =  K_pred*c*dt - K_pred*(data);  
                x = x + (-gradV(x)+bias)*dt + eta * sigma*sdt;
            end
            X(i,j) = x;
          
            
            Is = Is - bias * eta/ sigma * sdt;
            Id = Id - bias.^2 / sigma^2 *dt;
            

             if  x < -0.9 && x > -1.1 
                 time(i) = j;
                 pathfunc(i) = exp(-beta*j*dt)*exp(Is+0.5*Id); %weighted path functional 
                 X(i,j:end)=x;
                 break;
             else
                 pathfunc(i) = exp(1)*exp(Is+0.5*Id);
                 
             end
        end

    end
    
    data = X';
    data= data(:);

    K=zeros(nsteps*ntrjs,nsteps*ntrjs);

    pathweight = repmat(pathfunc,nsteps,1);
    pathweight= pathweight(:);
    
   
    for j=1:length(data)
            K(j,:) = (pathweight'/ntrjs).*sk/sqrt(2*pi)*l.^2.*exp(-0.5*(data'-data(j)).^2/l.^2);
    end
    
    A = -K*dt + 2/beta*speye(nsteps*ntrjs,nsteps*ntrjs); 
    b = K*sqrt(2/beta)*data;
    
    c = A\-b;
    
    c_old(opt,:) = c_pred;
    % plot zur Kontolle
    x_pred = linspace(-2,2,n_pred);
    K_pred = zeros(length(x_pred),nsteps*ntrjs);

   
   for i=1:length(x_pred)
       K_pred(i,:) = (pathweight'/((2/beta)*ntrjs)).* sk/sqrt(2*pi)*l.^2.*exp(-0.5*(data'-x_pred(i)).^2/l.^2);
   end
   c_pred = K_pred*c*dt - K_pred*(data);  
    


    figure(opt)
    plot(x_pred, c_pred)
    
    p=[sum(time>0),ntrjs];
    fprintf('Trajectories in T %d / %d \n', p )
    fprintf('|c_new-c_old|_2 = %f \n', norm(c_pred'-c_old(opt,:)))
    
end
%%
figure(1)
plot(x_pred, c_old(end,:),'LineWidth',3)
title('Control')

figure(6)
plot(x_pred,  -gradV(x_pred),x_pred, -gradV(x_pred) + c_old(end,:),'LineWidth',3 ); 
legend('-gradV','-gradV+cPred')
title('Gradients')

dx=x_pred(2)-x_pred(1);
per_pot = zeros(1,n_pred+1);
control = zeros(1,n_pred+1); 
 
for i=2:n_pred+1
    per_pot(i) = per_pot(i-1) + (-c_old(end,i-1) + gradV(x_pred(i-1)) )*dx;
    control(i) = control(i-1) -c_old(end,i-1)*dx;
end

% figure(7)
% plot(x_pred,V(x_pred),x_pred,per_pot(2:end)+14.5,'LineWidth',3); 
% legend('Potential','Perturbed Potential')
% title('Perturbed Potential')

figure(8)
plot(x_pred,V(x_pred),x_pred,V(x_pred)+control(2:end)+11,'LineWidth',3); 
legend('Potential','Perturbed Potential')

figure(9)
plot(x_pred,control(2:end),'LineWidth',3); 
legend('Predicted Control (Integral)')

