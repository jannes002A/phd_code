clear all

%% DRAW POTENTIAL

f = linspace(-1.6,1.6,1000);
V = @(x) 0.5*(x.^2-1).^2;

% figure
% plot(f,y)

%% DERIVATIVE
gradV = @(x) 2*x.*(x.^2-1);
y = V(f);

%% TRAJECTORY
dt = 0.0001;
sdt = sqrt(dt);
beta = 7;
sigma = sqrt(2/beta);

nvs = 1;
ntrjs = 1000; %number of trajectories
finalp = zeros(nvs,1);
finalpw = zeros(nvs,1);
finalpt = zeros(nvs,1);


nsteps = 15000;



%metadynamics
s = 0.8;
w =0.05/sqrt(2*pi*s^2);
% s = 0.2;
% w = 0.01/sqrt(2*pi*s^2);
delta_t = 200; %200;




    k = 1;
    x = -1; %starting point (left well)


    t = 1;    
    stop = 0;
X = zeros(60*nsteps,1);

    while x < 0 

        if mod(t,delta_t) == 0 
            k = k + 1;
        end

        x0 = X(1:delta_t:k*delta_t);
        exp_meta = w * exp( - 0.5*(x - x0).^2/s^2);
        dVx_meta = - sum(1/s^2*(x - x0).*exp_meta);%- sum(exp_meta);%
        x = x - gradV(x) * dt - dVx_meta*dt + sigma* randn*sdt;
        t = t + 1; 
        X(t) = x;

        if    x > 0.9 && x < 1.1 %x>-1.1 && x<-0.9%
            stop = 1;
        end
        
    end
    
    % Draw potential

    V_meta = zeros(1,1000);
    gV_meta = zeros(1,1000);
    
    
    for i=1:k
        for r=1:1000
            %x0 = X(1+i*delta_t);   
            V_meta(r) = V_meta(r) + w * exp(-0.5*(f(r) - x0(i))^2/s^2);
            gV_meta(r) = - sum(1/s^2*(f(r) - x0).*w.* exp( - 0.5*(f(r) - x0).^2/s^2));
        end  
    end
    %figure
    %plot(f,y,'b',f,y+V_meta,'r')
    %pause(1)

    t
    k

    %Girsanov weight
    W = zeros(1,ntrjs);

    %Estimators
   
    pt = zeros(1,ntrjs);   %E[exp(-beta*tau)Z]
    time = nsteps*ones(1,ntrjs); %time of each trajectory

    for n = 1:ntrjs
%         if mod(n,10)==0
%             %clc
%             disp([num2str(n/ntrjs*100), '%'])
%         end

        %Girsanov integrals
        
        I1 = 0;
        I2 = 0;
        t=0;
        x = -1;
        while x < 0
            eta = randn;
            exp_meta = w * exp( - 0.5*(x - x0).^2/s^2);
            dVx_meta = - sum(1/s^2*(x - x0).*exp_meta);%- sum(exp_meta);%


            x = x - (gradV(x) + dVx_meta) * dt + eta * sigma*sdt;
            I1 = I1 + dVx_meta * eta / sigma * sdt; %- nabla U(x)
            I2 = I2 + (dVx_meta).^2 / sigma^2 *dt;
            t=t+1;

        end
           pt(n) = exp(-beta*(t)*dt);
           time(n) = t;
           W(n) = exp( I1 - 0.5*I2); 
        
        
        pt(n) = pt(n) * W(n);
    end

   





fprintf('E[exp(-beta*tau)]: %2.8f \n',mean(pt))
fprintf('Var(E[exp(-beta*tau)]): %2.8f \n',var(pt))
fprintf('R(I): %2.8f \n',sqrt(var(pt))/mean(pt))


%%
%Plots
x = linspace(-1.6,1.6,1000);
V = @(x) 0.5*(x.^2-1).^2;


gradV = @(x) 2*x.*(x.^2-1);
Vx = V(x);
gV= gradV(x);
dvbias=zeros(1,length(x));
vbias=zeros(1,length(x));

for i=1:length(x0)
    dvbias = dvbias -1/s^2*(x - x0(i)).*w .* exp( - 0.5.*(x - x0(i)).^2/s^2);
    vbias = vbias + w .* exp( - 0.5.*(x - x0(i)).^2/s^2);
end


figure(1)
plot(x,dvbias,'g',x,-(gV+dvbias),'r',x,-gV,'b','LineWidth',3)

figure(2)
plot(x,vbias,'g',x,Vx+vbias,'r',x,Vx,'b','LineWidth',3)

figure(3)
plot(x,V_meta,'g',x,Vx+V_meta,'r',x,Vx,'b','LineWidth',3)
legend('Bias','Pot+Bias','Pot')

figure(4)
plot(x,V_meta,'k',x,Vx+V_meta,':k',x,Vx,'-.k','LineWidth',3)
legend('Bias','Pot+Bias','Pot')

