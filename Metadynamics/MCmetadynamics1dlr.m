%Monte Carlo Sampling for P(X_t in A |t<T)

%rng(2,'twister')                 % restore the previous settings


% %Range
x=linspace(-1.5,1.5);
%e.g. 
% poth =@(y) 8.*y.^4-44/3.*y.^3 + 2.*y.^2 +11/3.*y+1;
% dpoth =@(y) 32.*y.^3-44.*y.^2+4.*y.^2+11/3;
% pot = poth(x);
% dpot= dpoth(x);
% Xzero=-0.25; 

poth =@(y)  1/2.*(y.^2-1).^2;
dpoth =@(y) 2.*y.*(y.^2-1);
pot = poth(x);
dpot= dpoth(x);
Xzero=1;


beta=3;
Xtemp=Xzero;

Ns=1000;
%dt=1e-4;
%Tl=1/dt;
%Tl=1.5/dt;
%dt=1/10000;
Tl=15000;
dt=0.0001;

time=ones(1,Ns)*Tl*dt;
p=zeros(1,Ns);
pt=zeros(1,Ns);
Xhelp=zeros(1,Tl);

% Sampling um die Varianz des Schätzers zu betrachten 

 rng(3,'twister')

 for i=1:Ns
     
     
     Xtemp=Xzero;
    
    for j=1:Tl
    
        
        % Sampling of the trajectory 
        dBt=sqrt(2*dt/beta)*randn;
        Xtemp = Xtemp - (dpoth(Xtemp))*dt + dBt;
        
        % #of trajectories that hit the set within the time
        if Xtemp > -1.1 && Xtemp < -0.9
            p(i)=1;
            pt(i)=exp(-beta*j*dt);
            time(i) = j*dt; 
            %Xhelp(1:j+1)=Xtemp(1:j+1);
            break;
        end
        
        
    end
    
     
   
    
 end
 
 
 %hits
 
 


fprintf('MCMetadynamics \n')
fprintf('P[X_t in B| t< T)]: %2.8f \n',mean(p))
fprintf('Var[P(X_t in B| t< T)]: %2.8f \n',var(p))
fprintf('Relative Error [P(X_t in B| t< T)]: %2.8f \n',sqrt(var(p))/mean(p))
fprintf('E[exp(-beta*tau)]: %2.8f \n',mean(pt))
fprintf('Var(E[exp(-beta*tau)]): %2.16f \n',var(pt))
fprintf('Relative Error E[exp(-beta*tau)]: %2.8f \n',sqrt(var(pt))/mean(pt))
fprintf('E[time]: %2.8f \n',mean(time))




figure(1)
plot(x,poth(x),'b','LineWidth',3), hold on
plot(Xhelp,poth(Xhelp),'r+'),%hold off
hold off
legend('Pot','X_t')
