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
ntrjs = 100000; %number of trajectories
finalp = zeros(nvs,1);
finalpw = zeros(nvs,1);
finalpt = zeros(nvs,1);


nsteps = 15000;



%metadynamics
% s = 0.8;
% w = 0.05/sqrt(2*pi*s^2);
s = 0.2;
w =0.01/sqrt(2*pi*s^2);
delta_t = 200; %200;


for v=1:nvs
    rng(v,'twister')
    k = 1;
    x = -1; %starting point (left well)


    t = 1;    
    stop = 0;
X = zeros(60*nsteps,1);

    while stop==0 

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

    for i=1:k
        %for r=1:1000
            %x0 = X(1+i*delta_t);   
            %V_meta(r) = V_meta(r) + w * exp(-0.5*(f(r) - x0(i))^2/s^2);
            V_meta = V_meta + w.* exp(-0.5.*(f - x0(i)).^2/s.^2);
        %end  
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
        I = 0;
        X = zeros(nsteps,1);

        x = -1;
        first = 0;
        Vx0 =sum( w * exp( - 0.5*(x - x0).^2/s^2));
        
        for t = 1:nsteps-1 
            exp_meta = w * exp( - 0.5*(x - x0).^2/s^2);
            dVx_meta = - sum(1/s^2*(x - x0).*exp_meta);%- sum(exp_meta);%
            dVxx_meta = sum(1/s^4.*(x.^2 + x0.^2-2.*x0.*x-s.^2).*exp_meta);

            x = x - (gradV(x) + dVx_meta) * dt + eta(t) * sigma*sdt;
            X(t) = x;
            I = I + 1/(2*beta^-1).* (dVx_meta*gradV(x)+ 1/2*dVx_meta*dVx_meta - beta^-1*dVxx_meta)*dt;

            if  x > 0.9 && x < 1.1 && first == 0
            %if x > -1.1 && x < -0.9 && first == 0

                p(n) = 1;
                pt(n) = exp(-beta*(t)*dt);
                time(n) = t - 1;
                first = 1;
                Vxt = sum(w * exp( - 0.5*(x - x0).^2/s^2));
                break;
            end
        end


        if p(n)== 0
            W(n) = 0;         
        else
            W(n) = exp(I+ 1/(2*beta^-1)*(Vxt-Vx0)); 
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
plot(x,dvbias,'b',x,-gV,'r',x,-(gV+dvbias),'g','LineWidth',3)

figure(2)
%plot(x,-vbias,'b',x,Vx,'r')
plot(x,vbias)
figure(3)
plot(x,V_meta)

