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
numBasisFun=50;
%Find expansion coefficients of c wrt 50 Gaussian basis functions
% x_ and c_ are the grid and control for points in (-1,2)
leftlim=-1.1;
rightlim=1.9;
[opt_coeffs,basisvalues]=getcoeffs(numBasisFun,leftlim,rightlim,x(KK),c(KK_));

%Define number of desired iterations of gradient descent method
num_iter=10;
%Set the default stepsize for gradient descent
graddesc_stepsize =  0.004;%0.003

%Define matrix of control coefficient iterates
coeffs=zeros(num_iter+1,numBasisFun);

%As our initial guess of a control policy, perturb opt_coeffs by
%translation
coeffs(1,:)=opt_coeffs+0.05*ones(size(opt_coeffs));
%Cheat a little by making sure the control vector is all negative
% for i=1:numBasisFun
%     coeffs(1,i)=coeffs(1,i)*(-sign(coeffs(1,i)));
% end

%Define a vector for containing the (monotonically decreasing) value
%function / free energy function, determined by the control policy
F_at_x0=zeros(num_iter+1,1);

%Define the number of samples for the sampling method
n=200;
%Define the set S=[-1.1,-1] which we use in computing the hitting time
a=-1.1;
b=-1.0;
%Specify the initial starting point for all trajectories
x0=1;
VarHit=zeros(1,num_iter);

%For loop corresponding to gradient descent
for k=1:num_iter
    
    %Compute:    
    % Ecost - expected cost of control,
    % E[ \int_0^\tau \frac{1}{2} \vert c(X_s) \vert^2 ds ]
    
    % Ehit - expected hitting time 
    % E[ \tau ( x0, [-1.1,-1] ) ] 
    
    % grad_tildeI - gradient of discretized cmlt. gen. fun \tilde{I}
    [Ecost,Ehit,Egrad,VHit,EstIn] = mean_controlled_hittingtimeSG( x0,BETA,a,b,POT,coeffs(k,:),leftlim,rightlim,n);
    grad_tildeI = Egrad;
    VarHit(k)=VHit;
    
    % F_approx(k) - is the value function / cumulant generating function
    % F_approx(k) = E_{Q_x} [ \int_0^\tau 1 + 0.5|c(X_s)|^2 ds ]
    F_at_x0(k)=Ecost+Ehit+EstIn; 
 
    %Obtain newest iterate of the control expansion coefficients 
    % by performing gradient descent. Store.
    coeffs(k+1,:) = coeffs(k,:) - graddesc_stepsize*grad_tildeI';
    disp(k)    
end
%%
strjs=1000;
dt = 0.001;
nsteps= 15000;
sdt=sqrt(dt);
time = zeros(strjs,1);
pathfunc =zeros(strjs,1);
girsanov = zeros(strjs,1);

for i = 1:strjs

        Is=0;
        Id=0;

        x = x0;
        
        for j = 2:nsteps
            eta=randn(1);
            Bx= Basisfunc(x,numBasisFun,leftlim,rightlim);
            cc = coeffs(end,:)*Bx';    
            
             
            x = x+( - dPot(x,POT) + sqrt(2)*cc)*dt + sqrt(2*dt/BETA)*eta;
            
            Is = Is - cc * eta/ sqrt(2/BETA) * sdt;
            Id = Id - cc.^2 / sqrt(2/BETA)^2 *dt;
            

             if  x > a && x < b 
                 time(i) = j;
                 pathfunc(i) = exp(-beta*j*dt)*exp(Is+0.5*Id); %weighted path functional 
                 girsanov(i) = exp(Is+0.5*Id);
                 break;
             else
                 pathfunc(i) = 0;
             end
        end

end

fprintf('E[P(X_t in B| t< T)]: %2.8f \n',mean(girsanov))
fprintf('Var[P(X_t in B| t< T)]: %2.8f \n',var(girsanov))
fprintf('R(I): %2.8f \n',sqrt(var(girsanov))/mean(girsanov))
fprintf('\n')
fprintf('E[exp(-beta*tau)]: %2.8f \n',mean(pathfunc))
fprintf('Var(E[exp(-beta*tau)]): %2.8f \n',var(pathfunc))
fprintf('R(I): %2.8f \n',sqrt(var(pathfunc))/mean(pathfunc))    
