%V=@(x) 1/2.*x.^4-x.^2 - 0.2*x+0.3;
%gradV = @(x) 2*x.*(x.^2-1)-0.2;
V=@(x) 1/2.*(x.^2-1).^2;
gradV = @(x) 2*x.*(x.^2-1);
%pot = V(x);
%dpot= gradV(x);
Xzero=-1;


beta=5;
epsilon=1;
gamma = 1;
beta2 = 50;
sdt=sqrt(dt);


Tl=1000000;
dt=0.0001;

x = zeros(Tl,1);
x2 = zeros(Tl,1);
x(1)= Xzero;
x2(1)=Xzero;

y=0;

    
    for j=2:Tl
    
        
        % Sampling of the trajectory 
        dBt=sqrt(2*dt/beta)*randn;
        x(j) = x(j-1) - (gradV(x(j-1)-y))*dt + dBt;
        
        y= y - 1/(epsilon*gamma)*y*dt + 1/sqrt(epsilon*beta2)*sdt*randn(1);
        
        dBt=sqrt(2*dt/beta)*randn;
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
    

