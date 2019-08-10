function [Ecost,Ehit,EgradJ,EgradSh] = mean_controlled_hittingtime( x0,BETA,a,b,POT,c,leftlim,rightlim,n)
%Compute the mean hitting time over n samples

avg_cost=0.;
avg_hit=0.;
avg_gradJ=0.;
avg_gradSh=0.;

for k=1:n
    [cost_add,hit_add,gradJ_add,gradSh_add] = controlled_hittingtime(x0,BETA,a,b,POT,c,leftlim,rightlim);
    avg_cost = avg_cost + cost_add;
    avg_hit = avg_hit + hit_add;
    avg_gradJ = avg_gradJ+gradJ_add;
    avg_gradSh = avg_gradSh + gradSh_add;
end

EgradJ = avg_gradJ/n;
EgradSh = avg_gradSh/n;
Ecost = avg_cost/n;
Ehit = avg_hit/n;

end

