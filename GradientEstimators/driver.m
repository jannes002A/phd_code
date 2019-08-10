function d=driver(c,x,leftlim,rightlim)
%Input: row vector of control coefficients 'c' and position 'x'
%Output: a scalar d, which gives the control c(x) at point 'x'

%Uses: Basisfunc

nB=length(c);
B=Basisfunc(x,nB,leftlim,rightlim);
d=c*B';

    
