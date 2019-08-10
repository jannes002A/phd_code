% Likelihood approach to SOC

beta=4; % noise parameter
%Potential and Derivative
V =@(x) (x.*x-1).^2 - 0.2*x + 0.3;
gradV =@(x)  4*x.*(x.*x-1) - 0.2;

nB = 50; %number of Basis functions
leftlim = -1.1; % boundary of metastable set
rightlim = 1.9;

% boundary of the hitting set
a=-1.1;
b=-1.0;
num_iter = 10;
num_traj = 10;

%start values
x0=1; % starting point
sl=0.003; % steeplength
% get best approximation of the HJB solution
xh = linspace(-1,2,1000);
c=getcoeffs(nB, leftlim,rightlim,xh',pot);

% perturb the initial guess a bit
c = c+ 0.05*ones(size(c));

objfunc= zeros(1,num_traj);
derivative1 = zeros(1,num_traj);
derivative2 = zeros(1,num_traj);
hitting_time = zeros(1,num_traj);

% loop for optimization

for k=1:num_iter
    
    disp(k)
    
    % Compute expected value of
    % E[ \tau + 1/4 \int |control|^2 ds]
    
    
    
    % Loop for computing the mean first hitting time
    for i=1:num_traj
    
    dt= 0.001;
    sigma = sqrt(2/beta);
    x=x0;
    cost=0;
    W1=0;
    W2=0;
    t=0;
    
        while (x>b) || (x<a)
            
        % basis function evaluated at x
        Bx = Basisfunc(x,nB,leftlim,rightlim);
        % control
        control = c'*Bx';
            
        % cost    
        cost = cost + 0.25 * control.^2 *dt;
        
        %random number
        rN =randn;
        % Brownian motion
        dB = sigma * sqrt(dt)* rN;   
        % Euler- Maruyama for the SDE
        x= x + (- gradV(x) + sqrt(2)* control)* dt * dB;
        
        % Derivative
        W1 = W1 + 0.5 * (control.*Bx')*dt;
        W2 = W2 -(t+cost)* ((sqrt(2)*Bx')/sigma * sqrt(dt)*rN);
        
        % time 
        t=t+dt;
        
        end
    % the objective function
    hitting_time(i) = t;
    objfunc(i) = t+ cost;
    derivative1(i) = W1;
    derivative2(i) = W2;
    
    end
    
    %average of the objective function and the derivatives
    Eof = mean(objfunc);
    Ed1 = mean(W1);
    Ed2 = mean(W2);
    
    %gradient descent
    c = c- sl*(Ed1+Ed2);
    

end