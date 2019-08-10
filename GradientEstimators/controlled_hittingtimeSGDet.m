function [cost,hitting_time,grad,gir] = controlled_hittingtimeSGDet(x0,BETA,a,b,POT,c,leftlim,rightlim)
%This function computes the hitting time of a controlled random process
%which lives on a potential energy landscape given by POT of a set (a,b)
%and returns the cost of control, the gradient of the free energy F and the
%gradient of the Onsager-Machlup function Sh

%Input: 
% POT - potential,
% temperature - BETA, 
% starting point - x0, 
% target set - (a,b)
% control - c

%Output: exit time 'extime', cost of control, 'cost'

%Uses: driver.m
%Used by: 

dt = 0.001;
noise_prefactor=sqrt(2*dt/BETA);
prefactor= sqrt(2/BETA);
x=x0;
t=0;
nB=length(c);
sig=1;
cost=0;
dcc=0;
%


[~,~,Estart]= BasisfuncX(x,nB,leftlim,rightlim);
W1= Estart;

while (x>b)||(x<a)
    % Bx = Basisfunc(x,nB) gives the vector (b_1,\ldots,b_nB)(x)
    %      of basis function values at x
    [Bx,dB]= BasisfuncX(x,nB,leftlim,rightlim); 
    
    
    
    % cc gives the control policy c(x) 
    cc = c*Bx';
    cch= Bx;
    
    
    
    % CostControl=\int_0^\tau \frac{1}{2} c(X_t)dt
    cost = cost + 0.5*cc^2*dt;
    % Derivative of the cost wrt the parameter
    % d/da 1/2|a*f(x)|^2 = a*f(x)^2;
    dcc =  dcc  + cc.*cch*dt;

    %random_num is a standard normal random variable
    random_num=randn;
    %eta is normally-distributed scalar random variable with mean zero and
    %variance = 2 h \varepsilon
    eta = noise_prefactor*random_num;
    
    % 1/sqrt(2*beta)^2*(int -nabla V(x_t)*nabla C(x_t)-|nabla V(x_t)|^2 + 1/beta delta(x_t)
    % dt-C(x_t)+C(x_0)) 
    % Laplace Basisfunction
    lB =  Bx.*dB;
    dPx=dPot(x,POT);
    %Vorzeichen!
    W1 = W1 + ( (dPx.*Bx) - 2 *(cc).*Bx - 1/BETA.*lB).*dt;
    

    
    %Euler-Maruyama scheme for SDE
    x = x+( - dPx + sqrt(2)*cc)*dt + eta;
    t=t+dt;
    
    %If set still not yet reached, extend the terminal time
    if t>sig*100
        disp(['t=' num2str(t), ', x=' num2str(x)]);
        sig=sig+1;
    end
end
%Store the first hitting time 
hitting_time = t;

%The gradient is expressed in terms of Girsanov then you differentiate wrt
%to the pertuabtion which is t*W1 for the first function for the second
%function you have to apply a chain rule this is why you have
%(cost+dcost)*W1

[~,~,Eend]= BasisfuncX(x,nB,leftlim,rightlim); 
W1=W1-(Eend);
W1=W1*sqrt(2)/prefactor^2;


grad = dcc' + (t+cost).*W1';

