clear all

%% DRAW POTENTIAL


xp = linspace(-2.6,2.6,1000);
V = @(x) 0.5*(x.^2-1).^2;

c=0;
fac=1;
l=1;
%Plot Mollifier
lower = (xp < c-0.99);
lower = find(lower,1,'last');
upper = (xp > c+0.99);
upper = find(upper,1,'first');
xm = linspace((c-0.99)/l^2,(c+0.99)/l^2,upper-lower);
m= @(x) exp(-1./(1-abs((x-c)/l).^2));
mp=zeros(1,1000);
mp(lower:upper-1) = fac*m(xm)./sum(m(xm));

s = 0.2;
w = 0.01/sqrt(2*pi*s^2);
exp_meta = w .* exp( - 0.5.*(xp - 0).^2/s^2);

figure(1)
plot(xp,mp,xp,exp_meta)

% figure(2)
 gradV = @(x) 2*x.*(x.^2-1);
% y = V(xp);

% plot(xp,-gradV(xp),xp,mp,xp,mp-gradV(xp))


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

xm = linspace(-0.99, 0.99,1000);
const =  sum(exp(-1./(1-abs(xm).^2)));

%metadynamics

delta_t = 200; %200;

for v=1:nvs
    rng(1,'twister')
    k = 0;
    x = -1; %starting point (left well)


    t = 1;    
    stop = 0;
    X = zeros(60*nsteps,1);
    mh= zeros(60*nsteps,1);

    while stop==0 

        if mod(t,delta_t) == 0 
            k = k + 1;
        end
        
        dVx_meta=0;
       
        if k > 0
            c = X(delta_t:delta_t:k*delta_t);
            for l=1:k
                if (abs(x-c(l)) < 1) 
                    mh(l) =   w .* exp( - 0.5.*(x - 0).^2/s^2);
                else
                    mh(l) =0;
                end
            end
            dVx_meta = sum(fac*mh(1:k)./(const));
        end
  
        x = x - gradV(x) * dt + dVx_meta*dt + sigma* randn*sdt;
        t = t + 1; 
        X(t) = x;
        


        if    x > 0.9 && x < 1.1 %x>-1.1 && x<-0.9%
            stop = 1;
        end
        
    end
    
    % Draw potential

%     V_meta = zeros(1,1000);
% 
%     for i=1:k
%         for r=1:1000
%             %x0 = X(1+i*delta_t);   
%             V_meta(r) = V_meta(r) + w * exp(-0.5*(f(r) - x0(i))^2/s^2);
%         end  
%     end
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
           
            for l=1:k
                if (abs(x-c(l)) < 1) 
                    mh(l) =   w .* exp( - 0.5.*(x - 0).^2/s^2);
                else
                    mh(l) =0;
                end
            end
            dVx_meta = sum(fac*mh(1:k)./(const));

            x = x - gradV(x) * dt + dVx_meta*dt + sigma* randn*sdt;
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
            W(n) = exp(-I1 - 0.5*I2); 
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
vbias=zeros(1,length(x));
mh=zeros(1000,k);

for j=1:1000 
    for l=1:k
        if (abs(x(j)-c(l)) < 1) 
            mh(j,l) =   w .* exp( - 0.5.*(x(j) - 0).^2/s^2);
        else
            mh(j,l) = 0;
        end
    end
end
dvbias = sum(fac*mh./(const),2);


figure(1)
plot(x,dvbias,'g',x,-gV+dvbias','r',x,-gV,'b','LineWidth',3)

figure(2)
plot(x,dvbias)

% figure(2)
% plot(x,vbias,'g',x,Vx+vbias,'r',x,Vx,'b','LineWidth',3)
% 
% figure(3)
% plot(x,V_meta,'g',x,Vx+V_meta,'r',x,Vx,'b','LineWidth',3)
% legend('Bias','Pot+Bias','Pot')

