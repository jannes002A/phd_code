%Gaussian prioir
n=50;
Xtest= linspace(-5,5,n);
l=0.1;
sf=1;

K=zeros(n,n);

for i=1:length(Xtest)
    for j=1:length(Xtest)
        K(i,j)= sf^2*exp(-(Xtest(i)-Xtest(j))^2/(2*l));
    end
end

L=chol(K+1e-15*eye(n));

data = randn(n,3);

f_prior = L*data;

% figure(1)
% plot(Xtest,f_prior);


% Predicting

Xtrain = [-4,-3,-2,-1, 1];
ytrain = sin(Xtrain);

K=zeros(length(Xtrain),length(Xtrain));
sf=1;
l=0.1;

for i=1:length(Xtrain)
    for j=1:length(Xtrain)
        K(i,j)= sf^2*exp(-(Xtrain(i)-Xtrain(j))^2/(2*l));
    end
end
L=chol(K+0.00005*eye(length(Xtrain)));

Ks=zeros(length(Xtrain),n);

for i=1:length(Xtrain)
    for j=1:length(Xtest)
        Ks(i,j)= sf^2*exp(-(Xtrain(i)-Xtest(j))^2/(2*l));
    end
end

Lk=K\Ks;
mu= Lk'* (L\ytrain');

figure(1)
plot(Xtest,mu); hold on
plot(Xtrain,ytrain,'x')
hold off