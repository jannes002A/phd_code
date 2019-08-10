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
beta = 3;
sigma = sqrt(2/beta);

nvs = 1;
ntrjs = 1000; %number of trajectories
finalp = zeros(nvs,1);
finalpw = zeros(nvs,1);
finalpt = zeros(nvs,1);


nsteps = 15000;



%metadynamics
delta_t = 200; %200;
dVx_meta=0;

for v=1:nvs
    rng(v,'twister')
    k = 0;
    x = -1; %starting point (left well)
    
    x0=zeros(1,nsteps/delta_t);
    s = zeros(1,nsteps/delta_t);

    Xhelp=zeros(1,delta_t);
    count=1;

    t = 1;    
    stop = 0;
    X = zeros(60*nsteps,1);

    while stop==0 

        if mod(t,delta_t) == 0 
            k = k + 1;
            count=1;
            x0(k) = mean(Xhelp);
            s(k)= max(abs(Xhelp-x0(k)))*3;
            Xhelp=zeros(1,delta_t);
        end
        if k>0
        %x0 = X(1:delta_t:k*delta_t);
        w =0.05./sqrt(2*pi.*s(1:k).^2);
        exp_meta = w .* exp( - 0.5*(x - x0(1:k)).^2./s(1:k).^2);
        dVx_meta = - sum(1./s(1:k).^2.*(x - x0(1:k)).*exp_meta);%- sum(exp_meta);%
        end
        x = x - gradV(x) * dt - dVx_meta*dt + sigma* randn*sdt;
        Xhelp(count)=x;
        t = t + 1; 
        %X(t) = x;
        count=count+1;
        
        if    x > 0.9 && x < 1.1 %x>-1.1 && x<-0.9%
            stop = 1;
        end
        
    end
    
    % Draw potential

    V_meta = zeros(1,1000);

    for i=1:k
        for r=1:1000
            %x0 = X(1+i*delta_t);   
            V_meta(r) = V_meta(r) + 0.05./sqrt(2*pi.*s(i).^2)*exp(-0.5*(f(r) - x0(i)).^2./s(i).^2);
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

        x = -1;
        first = 0;
        for t = 1:nsteps-1 
        w =0.05./sqrt(2*pi.*s(1:k).^2);
        exp_meta = w .* exp( - 0.5*(x - x0(1:k)).^2./s(1:k).^2);
        dVx_meta = - sum(1./s(1:k).^2.*(x - x0(1:k)).*exp_meta);%- sum(exp_meta);%


            x = x - (gradV(x) + dVx_meta) * dt + eta(t) * sigma*sdt;
            X(t) = x;
            I1 = I1 + dVx_meta * eta(t) / sigma * sdt; %- nabla U(x)
            I2 = I2 + (dVx_meta).^2 / sigma^2 *dt;

            if  x > 0.9 && x < 1.1 && first == 0
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
            W(n) = exp( I1 - 0.5*I2); 
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
% plot(x,-dvbias,'b',x,-gV,'r',x,-gV-dvbias,'g','LineWidth',3)
% 
% figure(2)
% plot(x,-vbias,'b',x,Vx,'r')



