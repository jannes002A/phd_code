function [Ecost,Ehit,stdCost,stdHit,EgradJ,EgradSh] = mean_std_controlled_hittingtime( x0,BETA,a,b,POT,c,leftlim,rightlim,n)
%Compute the mean hitting time over n samples

cost=zeros(n,1);
hitting_time=zeros(n,1);
avg_gradJ=0.;
avg_gradSh=0.;
for k=1:n
    [cost(k),hitting_time(k),gradJ_add,gradSh_add] = controlled_hittingtime(x0,BETA,a,b,POT,c,leftlim,rightlim);
    avg_gradJ = avg_gradJ+gradJ_add;
    avg_gradSh = avg_gradSh + gradSh_add;
end

EgradJ = avg_gradJ/n;
EgradSh = avg_gradSh/n;

Ecost = mean(cost);
stdCost = std(cost);

Ehit=mean(hitting_time);
stdHit=std(hitting_time);

end

