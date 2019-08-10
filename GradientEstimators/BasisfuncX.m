function [B,deltaB,Estart]=BasisfuncX(x,nB,leftlim,rightlim)

%This function computes a row vector B, such that B(k) is the density
%of a Gaussian with mean given by `means(k)' and variance 'sigma^2'

%Input: 'x' is a scalar. The position of some point in 1D. 'nB' is the
%number of basis functions.

%row vector of values of basis functions at point x
B=zeros(1,nB);
deltaB=zeros(1,nB);
Estart=zeros(1,nB);

%step size
h=(rightlim-leftlim)/(nB-1);
%Me: means are the centers of the Gaussians, evenly spaced on (a,b) with
%step size h
means = leftlim:h:rightlim;
%Sigma is the standard deviation, same for all the Gaussians
sigma = h;
for k=1:nB
    %Me: B(k) is a number. It is the evaluation of a Gaussian density with
    %mean given by means(k) and variance given by sigma^2 at the point x.
    B(k) = exp(-0.5*(x-means(k))^2/sigma^2)/(sigma*sqrt(2*pi));
    %Laplace of basisfunction
    deltaB(k) = -((x-means(k))/sigma^2); 
    Estart(k) = -0.5*erf(0.707107*(means(k)-x)/sigma);
end

 