function [control_coeffs,basisvalues] = getcoeffs(numBasisFun,leftlim,rightlim,x,c )
%Obtain the expansion coefficients of a control c defined over grid x with
%respect to Gaussian basis functions

%Input
% numBasisFun - number of basis functions
% leftlim     - left limit, grid of uniformly spaced means of b. functions
% rightlim    - right limi, ''             ''                  ''
% x           - grid on which c is defined
% c           - control vector defined on x, such that length(c)=length(x)

%Output
% control_coeffs - expansion coefficients
% basisvalues    - matrix such that the i-th row is the vector of values 
%                   of the basis functions evaluated at x(i)

basisvalues=zeros(length(x),numBasisFun);

for i=1:length(x)
    basisvalues(i,:)=Basisfunc(x(i),numBasisFun,leftlim,rightlim);
end

control_coeffs=basisvalues\c;

end

