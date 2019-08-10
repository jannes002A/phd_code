%Compute Mean First Hitting Time of the set [-1.1,-1] from
%different starting points on the grid 'x' stored in DWPOT.mat

%Method: sampling for sample size of 10 per starting point

clear all
load DWPOT.mat
N=length(x);
MFHT=zeros(size(x));

POT='asym2wellPot';
a=-1.1;
b=-1;
n=10;
for i=1:N
    MFHT(i)=mean_hitting_time(x(i),BETA,a,b,POT,n);
end
    
save MFHT_sampling_10samples.mat

clear all
load DWPOT.mat
K=find(abs(x)<2);
N=length(K);
MFHT=zeros(length(K),1);
POT='asym2wellPot';
a=-1.1;
b=-1;
n=50;
for i=1:N
    MFHT(i)=mean_hitting_time(x(K(i)),BETA,a,b,POT,n);
end

save MFHT_sampling_50samples.mat 

    

