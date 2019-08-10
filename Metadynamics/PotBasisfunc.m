function B=PotBasisfunc(x,nB,Nmeans, sigma,varplot)

%This function computes a row vector B, such that B(k) is the density
%of a Gaussian with mean given by `means(k)' and variance 'sigma^2'

%Input: 'x' is a scalar. The position of some point in 1D. 'nB' is the
%number of basis functions.

% This basisfunction makes sense because it does not split when
% differentiating

%row vector of values of basis functions at point x
if varplot==0
    B=zeros(1,nB);
    %step size

    for k=1:nB
        %Me: B(k) is a number. It is the evaluation of a Gaussian density with
        %mean given by means(k) and variance given by sigma^2 at the point x.
        B(k) = exp(-0.5*(x-Nmeans(k)).^2/sigma(k).^2)./(sigma(k).*sqrt(2*pi));
    end

else
    B=zeros(length(x),nB);

    for k=1:nB
    %Me: B(k) is a number. It is the evaluation of a Gaussian density with
    %mean given by means(k) and variance given by sigma^2 at the point x.
    B(:,k) = exp(-0.5.*(x-Nmeans(k)).^2./sigma(k).^2)./(sigma(k).*sqrt(2*pi));
    end
end

%plot(x,omega*exp(-0.5.*(x-0).^2./sigma^2)./(sigma*sqrt(2*pi)))
