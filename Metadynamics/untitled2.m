% Skipt zur Simulation mit Metadynamics mit Girsanov
% To do;
% Sampling and Reweigting
% Proof of Novikov condition
% How to put the dVBias maybe there is a clever way
%  

%prevS = rng(0,'v5normal'); % use legacy generator, save previous settings
%rng(1,'twister')                 % restore the previous settings


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

beta=5;
Tl=1000000;
dt=1/1000;

Xtemp=zeros(1,1000);
Xtemp(1)=Xzero; 
 
 for i=1:Tl
        dBt=sqrt(2*dt/beta)*randn;
        Xtemp(i+1) = Xtemp(i) - (dpoth(Xtemp(i)))*dt + dBt;
 end

 figure(1)
 plot(1:Tl,Xtemp(1:Tl));
 xlabel('Iteration')
 ylabel('Position')
 
 figure(2)
 plot(x,poth(x),'LineWidth',2)
 xlabel('x')
 ylabel('Potential')
 
 figure(3)
 plot(x,exp(-beta*poth(x))/sum(exp(-beta*poth(x))),'LineWidth',2)
 xlabel('x')
 ylabel('Boltzmann Distribution')
 

