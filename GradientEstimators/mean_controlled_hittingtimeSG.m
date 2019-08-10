function [Ecost,Ehit,Egrad,varHit,varCost,varObjf,varGrad] = mean_controlled_hittingtimeSG( x0,BETA,a,b,POT,c,leftlim,rightlim,n)
%Compute the mean hitting time over n samples

avg_cost=0;
avg_hit=0;
avg_grad= zeros(20,n);
varhit=zeros(1,n);
cost = zeros(1,n);
objfun= zeros(1,n);


for k=1:n
    [cost_add,hit_add,grad_add] = controlled_hittingtimeSG(x0,BETA,a,b,POT,c,leftlim,rightlim);
    avg_cost = avg_cost + cost_add;
    avg_hit = avg_hit + hit_add;
    avg_grad(:,k) = grad_add;
    varhit(k)= hit_add;
    cost(k)= cost_add;
    objfun(k) = hit_add+cost_add; 
    
end

Egrad = sum(avg_grad,2)/n;
Ecost = avg_cost/n;
Ehit = avg_hit/n;
varHit=var(varhit);
varCost= var(cost);
varObjf= var(objfun);
varGrad = var(avg_grad,0,2);



end