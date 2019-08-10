function [cost,hitting_time,gradJ,grad_Sh] = controlled_hittingtime(x0,BETA,a,b,POT,c,leftlim,rightlim)
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
x=x0;
t=0;
nB=length(c);
sig=1;
cost=0;
%
sum_partial_tilde_gh = zeros(size(c));
grad_Sh = zeros(size(c));

while (x>b)||(x<a)
    % Bx = Basisfunc(x,nB) gives the vector (b_1,\ldots,b_nB)(x)
    %      of basis function values at x
    Bx = Basisfunc(x,nB,leftlim,rightlim); 
    
    % cc gives the control policy c(x) 
    cc = c*Bx'; 
    
    % At the end of the while loop, the term below is equal to
    % \sum^{N_\tau-1}_{k=0} h c(\tilde{X}_k)\cdot b_j(\tilde{X}_k)
    %That is, it corresponds to the running cost
    % \sum^{N_\tau-1}_{k=0}\frac{\partial \tilde{g}_h}{\partial a_j} in the
    % Efficient Rare Event Simulation (ERES) paper by Schuette and Hartmann
    sum_partial_tilde_gh = sum_partial_tilde_gh + cc*Bx*dt;
    
    % CostControl=\int_0^\tau \frac{1}{2} c(X_t)dt
    cost = cost + 0.5*cc^2*dt;

    %random_num is a standard normal random variable
    random_num=randn;
    %eta is normally-distributed scalar random variable with mean zero and
    %variance = 2 h \varepsilon
    eta = noise_prefactor*random_num;
    
    %At the end of the while loop, the j-th component of the term below is 
    %\frac{\partial S_h}{\partial a_j} from the ERES paper.
    % Some prefactors may be incorrect.
    % Notes:  Int1 is a row vector the same size as Bx (i.e. length nB)
    grad_Sh = grad_Sh + random_num*Bx;
    
    %Euler-Maruyama scheme for SDE
    x = x+( - dPot(x,POT) + sqrt(2)*cc)*dt + eta;
    t=t+dt;
    
    %If set still not yet reached, extend the terminal time
    if t>sig*100
        disp(['t=' num2str(t), ', x=' num2str(x)]);
        sig=sig+1;
    end
end
%Store the first hitting time 
hitting_time = t;

%The term below is the gradient of S_h wrt the control coefficients (a_j)_j
% i.e. the j-th entry of nabla_Sh is equal to
% \frac{\partial S_h}{\partial a_j}
grad_Sh=grad_Sh*(-sqrt(dt*BETA));

%The term below is the sum 
%  \sum^{N_\tau-1}_{k=0}\frac{\partial \tilde{g}_h}{\partial a_j}
%    + 
% [ \sum^{N_\tau-1}_{k=0} \tilde{g}_h ] * [ \nabla S_h ]
%
% Note that gradJ is a vector with length equal to the length of
% the number of basis functions in the finite-dimensional ansatz space
gradJ = sum_partial_tilde_gh -(cost+hitting_time)*grad_Sh;
