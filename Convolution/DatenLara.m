%Skipt zur Simulation mit Homotopie



%close all
% %Range
x=linspace(-0.5,1.5);
%e.g. 
poth =@(y) 8.*y.^4-44/3.*y.^3 + 2.*y.^2 +11/3.*y+1;
dpoth =@(y) 32.*y.^3-44.*y.^2+4.*y+11/3;
%homopoth =@(y,t) poth(y)+96.*t.^2+(4-88.*y+96.*y.^2).*t;
homopoth =@(y,t) poth(y)+96.*(t^2/2).^2+(4-88.*y+96.*y.^2).*(t^2/2);
dhomopoth =@(y,lambda) dpoth(y)+(-88 +192.*y).*lambda^2/2; 


epsilon = 0.4;

lambda=[0.1,0.15,0.2,0.25,0.3,0.4,0.5];
l=length(lambda);


Xzero=-0.25; 
dt=1/1000;
Ns=10000; % trajectory length
N=1000;  % number of trajectories


data= zeros(Ns,N,l);
   
 for i=1:l
     Xtemp=zeros(N,Ns);
        Xtemp(:,1)=Xzero;
        for k=1:Ns-1
            Xtemp(:,k+1) = Xtemp(:,k)-(dhomopoth(Xtemp(:,k),lambda(i))).*dt + sqrt(2*dt*epsilon).*randn(N,1);
        end
    data(:,:,i)=Xtemp';
 end

% Die Daten sind als Matrix gespeichert. Mit data(:,:,i) nimmst du die
% Daten für das enstprechende lambda(i). In der Matrix data(:,:,i) sind die
% Trajektorien als Spaltenvektor gespeichert. Die Trajektorien haben eine
% Länge von Ns*dt Zeitschritten. Für jedes lambda(i) habe ich N Trajektoren
% berechent. 
save('ConvData','data','-v7.3')

