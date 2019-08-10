%This script tests the following functions:
% getcoeffs.m
% hittingtime.m
% controlled_hittingtime.m
% mean_hitting_time.m

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

a=-1.1;
b=-1;

x0=1;

[cost,hitting_time_c,gradJ,grad_Sh] = controlled_hittingtime(x0,BETA,a,b,POT,opt_coeffs',leftlim,rightlim);
hitting_time = hittingtime(x0,BETA,a,b,POT);
disp('First hitting time \tau (1, [-1.1,1] )')
disp(['  SDE, optim. control = ' num2str(hitting_time_c)]);
disp(['  SDE without control = ' num2str(hitting_time)]);

n=100;
x0=1;
[Ecost,Ehit,EgradJ,EgradSh]= mean_controlled_hittingtime( x0,BETA,a,b,POT,opt_coeffs',leftlim,rightlim,n);
mean_hit = mean_hitting_time( x0,BETA,a,b,POT,n );

disp('Mean first hitting time E_x[ \tau (1, [-1.1,1] ) ]')
disp(['  SDE, optim. control = ' num2str(Ehit)]);
disp(['  SDE without control = ' num2str(mean_hit)]);

save MFHTdata_100samples_x0eq1

x0=0.75;
[Ecost,Ehit,EgradJ,EgradSh]= mean_controlled_hittingtime( x0,BETA,a,b,POT,opt_coeffs',leftlim,rightlim,n);
mean_hit = mean_hitting_time( x0,BETA,a,b,POT,n );

disp('Mean first hitting time E_x[ \tau (0.75, [-1.1,1] ) ]')
disp(['  SDE, optim. control = ' num2str(Ehit)]);
disp(['  SDE without control = ' num2str(mean_hit)]);

save MFHTdata_100samples_x0eq_pt75

x0=1.25;
[Ecost,Ehit,EgradJ,EgradSh]= mean_controlled_hittingtime( x0,BETA,a,b,POT,opt_coeffs',leftlim,rightlim,n);
mean_hit = mean_hitting_time( x0,BETA,a,b,POT,n );

disp('Mean first hitting time E_x[ \tau (1, [-1.1,1] ) ]')
disp(['  SDE, optim. control = ' num2str(Ehit)]);
disp(['  SDE without control = ' num2str(mean_hit)]);

save MFHTdata_100samples_x0eq1pt25.mat 

x0=0.5;
[Ecost,Ehit,EgradJ,EgradSh]= mean_controlled_hittingtime( x0,BETA,a,b,POT,opt_coeffs',leftlim,rightlim,n);
mean_hit = mean_hitting_time( x0,BETA,a,b,POT,n );

disp('Mean first hitting time E_x[ \tau (1, [-1.1,1] ) ]')
disp(['  SDE, optim. control = ' num2str(Ehit)]);
disp(['  SDE without control = ' num2str(mean_hit)]);

save MFHTdata_100samples_x0eq0pt5

x0=1.5;
[Ecost,Ehit,EgradJ,EgradSh]= mean_controlled_hittingtime( x0,BETA,a,b,POT,opt_coeffs',leftlim,rightlim,n);
mean_hit = mean_hitting_time( x0,BETA,a,b,POT,n );

disp('Mean first hitting time E_x[ \tau (0.75, [-1.1,1] ) ]')
disp(['  SDE, optim. control = ' num2str(Ehit)]);
disp(['  SDE without control = ' num2str(mean_hit)]);

save MFHTdata_100samples_x0eq_1pt5

x0=1.75;
[Ecost,Ehit,EgradJ,EgradSh]= mean_controlled_hittingtime( x0,BETA,a,b,POT,opt_coeffs',leftlim,rightlim,n);
mean_hit = mean_hitting_time( x0,BETA,a,b,POT,n );

disp('Mean first hitting time E_x[ \tau (1, [-1.1,1] ) ]')
disp(['  SDE, optim. control = ' num2str(Ehit)]);
disp(['  SDE without control = ' num2str(mean_hit)]);

save MFHTdata_100samples_x0eq1pt75.mat 