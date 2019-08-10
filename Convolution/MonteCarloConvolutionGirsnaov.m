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

%V=@(x) 1/2.*x.^4-x.^2 - 0.2*x+0.3;
%gradV = @(x) 2*x.*(x.^2-1)-0.2;
V =@(y) 8.*y.^4-44/3.*y.^3 + 2.*y.^2 +11/3.*y+1;
gradV =@(y) 32.*y.^3-44.*y.^2+4.*y+11/3;


dt = 0.0001;
sdt = sqrt(dt);
beta = 3;
sigma = sqrt( 2/beta);



Ns=1000;
%dt=1e-4;
%Tl=1/dt;
%Tl=1.5/dt;
%dt=1/10000;
Tl=15000;


time=ones(1,Ns)*Tl*dt;
p=zeros(1,Ns);
pt=zeros(1,Ns);



% Sampling um die Varianz des Schätzers zu betrachten 

 parfor i=1:Ns
     
     
     Xtemp=-0.25;
    
    for j=1:Tl
    
        
        % Sampling of the trajectory 
        dBt=sqrt(2*dt/beta)*randn;
        Xtemp = Xtemp - (gradV(Xtemp))*dt + dBt;
        
        % #of trajectories that hit the set within the time
        if Xtemp > 1.15 && Xtemp < 1.25 
            p(i)=1;
            pt(i)=exp(-beta*j*dt);
            time(i) = j*dt; 
            %Xhelp(1:j+1)=Xtemp(1:j+1);
            break;
        end
        
        
    end
    
     
   
    
 end
 
 

fprintf('Mone Carlo Estimator \n')
fprintf('E[P(X_t in B| t< T)]: %2.8f \n',mean(p))
fprintf('Var[P(X_t in B| t< T)]: %2.8f \n',var(p))
fprintf('Relative Error [P(X_t in B| t< T)]: %2.8f \n',sqrt(var(p))/mean(p))
fprintf('E[exp(-beta*tau)]: %2.8f \n',mean(pt))
fprintf('Var(E[exp(-beta*tau)]): %2.8f \n',var(pt))
fprintf('Relative Error E[exp(-beta*tau)]: %2.8f \n',sqrt(var(pt))/mean(pt))
fprintf('E[time]: %2.8f \n',mean(time))


