%This script computes the estimate of the Mean First Hitting Time using
%2000 sampled trajectories of the controlled SDE for a grid of starting
%points. 

clear all
load MFHT_PDEFiniteElement.mat c x KK KK_ POT BETA epsilon 

%Define the number of basis functions in the finite-dimensional ansatz
%space
numBasisFun=50;
%Find expansion coefficients of c wrt 50 Gaussian basis functions
% x_ and c_ are the grid and control for points in (-1,2)
leftlim=-1.1;
rightlim=1.9;
[opt_coeffs,~]=getcoeffs(numBasisFun,leftlim,rightlim,x(KK),c(KK_));
n=2000;

a=-1.1;
b=-1;
x0=-1:0.25:1.5;

Ecost=zeros(length(x0),1);
Ehit=zeros(length(x0),1);
stdCost=zeros(length(x0),1);
stdHit=zeros(length(x0),1);
EgradJ=zeros(length(x0),numBasisFun);
EgradSh=zeros(length(x0),numBasisFun);

for i=1:length(x0)
   [Ecost(i),Ehit(i),stdCost(i),stdHit(i),EgradJ(i,:),EgradSh(i,:)]= mean_std_controlled_hittingtime( x0(i),BETA,a,b,POT,opt_coeffs',leftlim,rightlim,n); 
end

save MFHT_stats numBasisFun leftlim rightlim opt_coeffs n a b Ecost Ehit stdCost stdHit EgradJ EgradSh x0

%% Do plots

%Plot psi
h1=figure(1);
clf
plot(x0,Ehit)
xlabel('x')
ylabel('E_x[\tau] with optimal control')
title('Plot of Mean First Hitting Time E_x[\tau] with optimal control, sample size=2000')
title_fig1='Figures/MFHT_optcontrol_2000samples';
print(h1,'-dpng',title_fig1)
saveas(gcf,strcat(title_fig1,'.fig'));
print(h1,'-depsc',title_fig1)

hold on
for k=1:length(x0)
    %95 % confidence interval - [mean - 1.96 * stdev, mean + 1.96 * stdev]
    plot([x0(k) x0(k)],[Ehit(k)-1.96*stdHit(k) Ehit(k)+1.96*stdHit(k)],'r');
end
title('Plot of Mean First Hitting Time E_x[\tau] with optimal control, sample size=2000, \newline 95% confidence intervals')
title_fig2='Figures/MFHT_optcontrol_2000samples_95confinterval';
print(h1,'-dpng',title_fig2)
saveas(gcf,strcat(title_fig2,'.fig'));
print(h1,'-depsc',title_fig2)

