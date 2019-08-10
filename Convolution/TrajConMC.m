%V=@(x) 1/2.*x.^4-x.^2 - 0.2*x+0.3;
%gradV = @(x) 2*x.*(x.^2-1)-0.2;
V=@(x) 1/2.*(x.^2-1).^2;
gradV = @(x) 2*x.*(x.^2-1);

%pot = V(x);
%dpot= gradV(x);
Xzero=-1;


beta=5;
Tl=1000000;
dt=0.0001;
sdt=sqrt(dt);

lambda=0.3;
x = zeros(Tl,1);
x2 = zeros(Tl,1);
x(1)= Xzero;
x2(1)=Xzero;


    
    for j=2:Tl
    
        
        % Sampling of the trajectory 
        dBt=sqrt(2*dt/beta)*randn;
        
        xapp = x(j-1)+lambda.*randn(20,1);
        congradV = sum(gradV(xapp))/20;
        
        x(j) = x(j-1) - (congradV)*dt + dBt;
        
        
        %Monte Carlo
        %dBt=sqrt(2*dt/beta)*randn;
        x2(j) = x2(j-1) - (gradV(x2(j-1)))*dt + dBt;
        
        
        
        
        
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
 