%V=@(x) 1/2.*x.^4-x.^2 - 0.2*x+0.3;
%gradV = @(x) 2*x.*(x.^2-1)-0.2;
poth =@(y) 8.*y.^4-44/3.*y.^3 + 2.*y.^2 +11/3.*y+1;
dpoth =@(y) 32.*y.^3-44.*y.^2+4.*y+11/3;
homopoth =@(y,t) poth(y)+96.*(t^2/2).^2+(4-88.*y+96.*y.^2).*(t^2/2);
gradV =@(y,lambda) dpoth(y)+(-88 +192.*y).*lambda^2/2; 

%pot = V(x);
%dpot= gradV(x);
Xzero=-0.25;


beta=3;
Tl=1000000;
dt=0.0001;
sdt=sqrt(dt);

x = zeros(Tl,1);
x2 = zeros(Tl,1);
x(1)= Xzero;
x2(1)=Xzero;


    
    for j=2:Tl
    
        
        % Sampling of the trajectory 
        dBt=sqrt(2*dt/beta)*randn;
        x(j) = x(j-1) - (gradV(x(j-1),0.3))*dt + dBt;
        
        
        %Monte Carlo
        dBt=sqrt(2*dt/beta)*randn;
        x2(j) = x2(j-1) - (dpoth(x2(j-1)))*dt + dBt;
        
        
        
        
        
    end
    
    figure(1)
    hist(x,500)
    xlabel('Position')
    ylabel('Number of Visits')
    
    figure(2)
    hist(x2,500)
    xlabel('Position')
    ylabel('Number of Visits')
    
    figure(3)
    plot(1:length(x),x)
    xlabel('time')
    ylabel('V(x)')
    
    figure(4)
    plot(1:length(x2),x2)
    xlabel('time')
    ylabel('V(x)')
   
    

