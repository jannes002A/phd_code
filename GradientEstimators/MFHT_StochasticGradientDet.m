%clear all
%Load:
% - the vector of optimal control values 'c', calculated using the Finite
%Element method as shown in the FMPT_01.m script. The optimal control is
%defined on (-1,2)
% - the vector of grid points 'x'
% - the vector of indices of grid points belonging to (-1,2), 'KK'
% - the string which denotes which potential energy function we use, 'POT'
% - the inverse temperature, 'BETA'
% - the inverse of BETA, 'epsilon'
% - the cumulant generating function, 'F(x)' on [-3,3]
load MFHT_PDEFiniteElement.mat c x KK KK_ POT BETA V F

%Define the number of basis functions in the finite-dimensional ansatz
%space
numBasisFun=20;
%Find expansion coefficients of c wrt 50 Gaussian basis functions
% x_ and c_ are the grid and control for points in (-1,2)
leftlim=-1.1;
rightlim=1.9;
[opt_coeffs,basisvalues]=getcoeffs(numBasisFun,leftlim,rightlim,x(KK),c(KK_));

%Define number of desired iterations of gradient descent method
num_iter=15;
%Set the default stepsize for gradient descent
graddesc_stepsize = 0.01; %0.007

%Define matrix of control coefficient iterates
coeffs=zeros(num_iter+1,numBasisFun);

%As our initial guess of a control policy, perturb opt_coeffs by
%translation
coeffs(1,:)=opt_coeffs + 0.2*ones(size(opt_coeffs));
%Cheat a little by making sure the control vector is all negative
% for i=1:numBasisFun
%     coeffs(1,i)=coeffs(1,i)*(-sign(coeffs(1,i)));
% end

%Define a vector for containing the (monotonically decreasing) value
%function / free energy function, determined by the control policy
F_at_x0=zeros(num_iter+1,1);

%Define the number of samples for the sampling method
n=1;
%Define the set S=[-1.1,-1] which we use in computing the hitting time
a=-1.1;
b=-1.0;
%Specify the initial starting point for all trajectories
x0=1;
VarHit=zeros(1,num_iter+1);
VarCost = zeros(1,num_iter+1);
VarObjf = zeros(1,num_iter+1);
VarGrad = zeros(numBasisFun,num_iter+1);
AvHit = zeros(1,num_iter+1);

%For loop corresponding to gradient descent
for k=1:num_iter
    
    %Compute:    
    % Ecost - expected cost of control,
    % E[ \int_0^\tau \frac{1}{2} \vert c(X_s) \vert^2 ds ]
    
    % Ehit - expected hitting time 
    % E[ \tau ( x0, [-1.1,-1] ) ] 
    
    % grad_tildeI - gradient of discretized cmlt. gen. fun \tilde{I}
    %[~,~,Egrad,~,~,~] = mean_controlled_hittingtimeSGDet( x0,BETA,a,b,POT,coeffs(k,:),leftlim,rightlim,n);
    [Ecost,Ehit,Egrad,vHit,VCost,VObjf,VGrad] = mean_controlled_hittingtimeSGDet( x0,BETA,a,b,POT,coeffs(k,:),leftlim,rightlim,200);
    grad_tildeI = Egrad;
    VarHit(k)=vHit;
    VarCost(k)= VCost;
    VarObjf(k) = VObjf;
    VarGrad(:,k) = VGrad;
    AvHit(k)=Ehit;
    % F_approx(k) - is the value function / cumulant generating function
    % F_approx(k) = E_{Q_x} [ \int_0^\tau 1 + 0.5|c(X_s)|^2 ds ]
    F_at_x0(k)=Ecost+Ehit; 
 
    %Obtain newest iterate of the control expansion coefficients 
    % by performing gradient descent. Store.
    coeffs(k+1,:) = coeffs(k,:) - (graddesc_stepsize/(k+10))*grad_tildeI';
    disp(k)    
end

%Get last value for F_approx
[Ecost,Ehit,Egrad,vHit,VCost,VObjf,VGrad] = mean_controlled_hittingtimeSGDet( x0,BETA,a,b,POT,coeffs(num_iter+1,:),leftlim,rightlim,200);
F_at_x0(k+1)=Ecost+Ehit; 
VarHit(k+1)=vHit;
VarCost(k+1)= VCost;
VarObjf(k+1) = VObjf;
VarGrad(:,k+1) = VGrad;
AvHit(k+1)=Ehit;

%% Plots

%Figure 1
h1=figure(1);
plot(1:1:(num_iter+1),F_at_x0,'bo')
xlabel('iterate i')
ylabel('F^i(1)')
title('Free energy iterates F^i(x) at x=1')
%title_fig1='Figures/MFHT_GradientDescent_Fiterates';
%print(h1,'-dpng',title_fig1)
%saveas(gcf,strcat(title_fig1,'.fig'));
%print(h1,'-depsc',title_fig1)

%Figure 2
err_coeffs=(coeffs-ones(num_iter+1,1)*opt_coeffs');

h2=figure(2);
plot(1:1:length(opt_coeffs),err_coeffs);
xlabel('j-th coefficient')
ylabel('(c^i-c)_j')
title('Plot of error in coefficients (c^i-c)_j at x=1' )
%title_fig2='Figures/MFHT_GradientDescent_errcontrolcoeffs';
%print(h2,'-dpng',title_fig2)
%saveas(gcf,strcat(title_fig2,'.fig'));
%print(h2,'-depsc',title_fig2)

%Figure 3
abserr_coeffs=abs(err_coeffs);

h3=figure(3);
plot(1:1:length(opt_coeffs),abserr_coeffs);
xlabel('j-th coefficient')
ylabel('|(c^i-c)_j|')
title('Plot of absolute error in coefficients | (c^i-c)_j | at x=1')
%title_fig3='Figures/MFHT_GradientDescent_abserr_controlcoeffs';
%print(h3,'-dpng',title_fig3)
%saveas(gcf,strcat(title_fig3,'.fig'));
%print(h3,'-depsc',title_fig3)

%Figure 4
%zmat1(i,j)=value of i-th iterate of control at j-th grid point on 
%           [-1.1,1.9] (grid has 50 points)
zmat1=(basisvalues*coeffs')';
%zmat2(i,j)=difference between value of optimal control and i-th iterate at
%           j-th grid point on [-1.1,1.9]
zmat2=ones(num_iter+1,1)*(basisvalues*opt_coeffs)'-zmat1;
abs_zmat2=abs(zmat2);
x_=x(KK_);

h4=figure(4);
clf
plot(x_,zmat1(1,:),'b',x_,(basisvalues*opt_coeffs),'r',x_,zmat1(11,:),'k','LineWidth',2)
legend('Initial iterate','Optimal control','Last iterate');
xlabel('x')
ylabel('Control c(x)')
title('Plot of control')
%title_fig4='Figures/MFHT_GradientDescent_comparecontrols';
%print(h4,'-dpng',title_fig4)
%saveas(gcf,strcat(title_fig4,'.fig'));
%print(h4,'-depsc',title_fig4)

%Figure 5

%F_array(i,j)=value of i-th iterate of free energy at j-th grid point of
%the original grid x on [-3,3] 
length_x=length(x);
F_array=zeros(num_iter+1,length_x);
dx=x(2)-x(1);
[~,pos]=min(abs(x+1)); % find position of x==-1
h5=figure(5);
clf
plot(x(KK),V(KK),'b',x(KK),V(KK)+2*F(KK),'k','LineWidth',2)
hold on
xlabel('x')

for i=1:(num_iter+1)
    for j=(pos+1):length_x
        %Compute F(x) from c(x) using the fact that c(x)=-sqrt(2)* dF/dx
        F_array(i,j)=F_array(i,j-1)-(1/sqrt(2.0))*basis(coeffs(i,:),x(j),leftlim,rightlim)*dx;
    end
    plot(x(KK),V(KK)+2*F_array(i,KK)','r') 
    hold on
end
title('Plot of potential energy functions')
legend('Original potential','Optimally tilted potential');
%title_fig5='Figures/MFHT_GradientDescent_approximated_potentials';
%print(h5,'-dpng',title_fig5)
%saveas(gcf,strcat(title_fig5,'.fig'));
%print(h5,'-depsc',title_fig5)

%save MFHT_GradientDescent

figure(6)
plot(1:num_iter+1,VarHit,1:num_iter+1,VarCost,1:num_iter+1,VarObjf,'LineWidth',2)
legend('Hitting Time','Cost','Obj. Fun.')
title('Variances')

figure(7)
plot(1:num_iter,VarGrad(:,1:num_iter));
xlabel('Iteration');
ylabel('Variance')
title('Variance of gradient estimator for each component')