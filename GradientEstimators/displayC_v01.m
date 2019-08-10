function d=displayC_v01(c)
%c needs to be a row vector

%Me: set nB to be the length of input vector c 
nB=length(c);
%Me: define a grid of step size 0.05 from -2 to 2 (why this interval and
%why this step size)?
x=-2:0.05:2;

%Me: d was iteratively defined so I added the line below
d=zeros(length(x),1);
for k=1:length(x)
    %Bx gives a row vector the same length as c. Bx(k) gives the value of
    %the k-th (Gaussian) basis function at point x(k)
    Bx=Basisfunc(x(k),nB);
    
    %d(k) gives the control c(x(k)) at the point x(k)
    d(k)=c*Bx';
    
end

%Plot the control policy c(x) as a function of x
figure(12)
clf
plot(x,d);