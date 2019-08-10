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


beta=3;
%dt=1/10000;
dt=1e-4;
Xtemp=Xzero;

Ns=100;
%Tl=15000;
Tl=1/dt;
VL=1;

time=zeros(1,Ns);
p=zeros(1,Ns);
pt=zeros(1,Ns);
Xhelp=zeros(1,Tl);
prob=zeros(1,VL);
probtime=zeros(1,VL);


% Sampling um die Varianz des Schätzers zu betrachten 
for v=1:VL
    rng(v,'twister')
    
    p=zeros(1,Ns);
    pt=zeros(1,Ns);
 for i=1:Ns
     
     Xtemp=zeros(1,Tl+1);
     Xtemp(1)=Xzero;
    
    for j=1:Tl
    
        
        % Sampling of the trajectory 
        dBt=sqrt(2*dt/beta)*randn;
        Xtemp(j+1) = Xtemp(j) - (dpoth(Xtemp(j)))*dt + dBt;
        
        % #of trajectories that hit the set within the time
        if Xtemp(j+1)>0.9 && Xtemp(j+1)<1.1
            p(i)=1;
            pt(i)=exp(-beta*j*dt);
            time(i) = j*dt; 
            Xhelp(1:j+1)=Xtemp(1:j+1);
            break;
        end
        
        
    end
    
     
   
    
 end
 
 
 %hits
 prob(v)=sum(p)/Ns;
 probtime(v)= sum(pt(pt>0))/Ns;
 
end

%fprintf('P(X_t in B| t< T): %2.8f \n',prob)
%fprintf('E[exp(-beta*tau)]: %2.8f \n',probtime)

re = sqrt(var(pt(pt>0)))/(mean(pt(pt>0)))


% Varinaz des Schätzers
mprobtime = mean(probtime);
varprobtime = var(probtime);

fprintf('E[P(X_t in B| t< T)]: %2.8f \n',mean(prob))
fprintf('Var[P(X_t in B| t< T)]: %2.8f \n',var(prob))

fprintf('E[exp(-beta*tau)]: %2.8f \n',mprobtime)
fprintf('Var(E[exp(-beta*tau)]): %2.8f \n',varprobtime)

figure(1)
plot(x,poth(x),'b','LineWidth',3), hold on
plot(Xhelp,poth(Xhelp),'r+'),%hold off
hold off
legend('Pot','X_t')

figure(2)
plot(x,poth(x),'-.k','LineWidth',3), hold on
plot(Xhelp,poth(Xhelp),'k+'),%hold off
hold off
legend('Pot','X_t')
