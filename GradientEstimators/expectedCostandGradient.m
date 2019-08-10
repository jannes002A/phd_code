function grad_tildeI=expectedCostandGradient(Ecost,Ehit,EgradJ,EgradSh)
%This function computes the expected cost and gradient by taking
%averages over 'n' random trials.

%Output: expected cost, 'costC', expected exit time, 'expExtime', expected
%gradient of potential, 'expdPhi'

%Depends on: controlled_hittingtime (function), 
%   BETA (parameter), POT (data), x0 (grid), c (vector of coefficients of
%   control vector given Gaussian basis functions b_j. See Basisfunctions.m

avg_cost=0.;
avg_hitting_time=0.;
avg_gradJ=0.;
avg_gradSh=0.;

for k=1:n
    [cost_add,hitting_time_add,gradJ_add,gradSh_add] = controlled_hittingtime(x0,BETA,a,b,POT,c,leftlim,rightlim);
    avg_cost = avg_cost + cost_add;
    avg_hitting_time = avg_hitting_time + hitting_time_add;
    avg_gradJ = avg_gradJ+gradJ_add;
    avg_gradSh = avg_gradSh + gradSh_add;
end

EgradJ = avg_gradJ/n;

EgradSh = avg_gradSh/n;
Ecost = avg_cost/n;
Ehit = avg_hitting_time/n;

grad_tildeI = EgradJ + (Ecost+Ehit)*EgradSh;

