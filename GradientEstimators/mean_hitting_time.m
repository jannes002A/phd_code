function mean_hit = mean_hitting_time( x0,BETA,a,b,POT,n )
%Compute the mean hitting time over n samples

mean_hit=0.0;
n=round(n);
for j=1:n
   mean_hit = mean_hit + hittingtime(x0,BETA,a,b,POT);
end

mean_hit=mean_hit/(1.0*n);

end

