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
Xzero=-1;


beta=7;
Xtemp=Xzero;

Ns=1000;
dt=0.0001;
VL=1;








% Sampling um die Varianz des Schätzers zu betrachten 

   
 time=zeros(1,Ns);  
 pt=zeros(1,Ns);
 for i=1:Ns
     
     
    Xtemp=Xzero;
    t=0; 
    while Xtemp<0
    
        
        % Sampling of the trajectory 
        dBt=sqrt(2*dt/beta)*randn;
        Xtemp = Xtemp - (dpoth(Xtemp))*dt + dBt;
        
        % #of trajectories that hit the set within the time
        
        t=t+1;    
        
            
        
        
        
    end
    pt(i)=exp(-beta*t*dt);
    time(i) = t*dt; 
     
   
    
 end
 
 
fprintf('E[exp(-beta*tau)]: %2.8f \n',mean(pt))
fprintf('Var(E[exp(-beta*tau)]): %2.8f \n',var(pt))
fprintf('Relative Error E[exp(-beta*tau)]: %2.8f \n',sqrt(var(pt))/mean(pt))
fprintf('E[time]: %2.8f \n',mean(time))


