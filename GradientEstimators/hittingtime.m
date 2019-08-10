function hitting_time = hittingtime(x0,BETA,a,b,POT)
%This function computes the hitting time of an uncontrolled random process
%which lives on a potential energy landscape given by POT of a set (a,b)

%Input: 
% POT - potential,
% temperature - BETA, 
% starting point - x0, 
% target set - (a,b)
% control - c

dt = 0.001;
x=x0;
t=0;
sig=1;

while (x>b)||(x<a)
    %eta is normally-distributed scalar random variable with mean zero and
    %variance = 2 h \varepsilon
    eta = sqrt(2*dt/BETA)*randn;
    
    %Euler-Maruyama scheme for SDE
    x = x - dt*dPot(x,POT) + eta;
    t=t+dt;
    
    %If set still not yet reached, extend the terminal time
    if t>sig*100
        disp(['t=' num2str(t), ', x=' num2str(x)]);
        sig=sig+1;
    end
end
%Store the first hitting time 
hitting_time = t;
