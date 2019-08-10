
%function [A]=infgenerator2D (n,x,y,sigma)

% function for calculating the inf gen in 2D
% n is the number of grid points
% dV is the potential evaluated at this grid
% x is a verctor with start and end point in x direction
% y is a verctor with start and end point in y direction
% sigma is the factor of the Brownian motion
x=[-2,0];
y=[-2,1];
%x=[-2,2];
%y=[-2,2];
sigma=0.3;
n=100-1;
%n=3;

V1 = @(x,y)  3.*exp(-x.^2-(y-1/3).^2)...
            -3.*exp(-x.^2-(y-5/3).^2)...
            -5.*exp(-(x-1).^2-y.^2)...
            -5.*exp(-(x+1).^2-y.^2) + ...
            1/5.*x.^4 +1/5*(y-1/3).^4;

dxV1 =@(x,y) - 6.*x.*exp(-x.^2-(y-1/3).^2)...
             + 6.*x.*exp(-x.^2-(y-5/3).^2)...
             +10.*(x-1).*exp(-(x-1).^2-y.^2) ...
             +10.*(x+1).*exp(-(x+1).^2-y.^2) + ...
             + 4/5.*x.^3;

dyV1=@(x,y) - 6.*(y-1/3).*exp(-x.^2-(y-1/3).^2)+ ...
            + 6.*(y-5/3).*exp(-x.^2-(y-5/3).^2) + ...
            +10.*y.*exp(-(x-1).^2-y.^2) + ...
            +10.*y.*exp(-(x+1).^2-y.^2) + ...
            + 4/5*(y-1/3).^3;        

x=linspace(x(1),x(2),n);
y=linspace(y(1),y(2),n);

[X,Y]=meshgrid(x,y);

V=V1(X,Y);
dVx=dxV1(X,Y);
dVx=dVx(:);
dVy=dyV1(X,Y);
dVy=dVy(:);
figure(1)
surf(X,Y,V)



%Laplacian
%derivative in x direction
cell1 =1/dx^2*sparse(diag(-2*ones(n,1))+diag(ones(n-1,1),-1)+diag(ones(n-1,1),1));
tridiag=kron(eye(n),cell1);
%derivative in y direction
cell2 =1/(dy^2)*sparse(diag(-2*ones(n^2,1))+diag(-1*ones(n^2-n,1),-n)+diag(-1*ones(n^2-n,1),n));
A=(sigma^2/2)*sparse(tridiag + cell2);

%first derivatives
cellx = spdiags(ones(n,1)*[-1 0 1],[1 0 -1],n,n);
tridiagx=kron(eye(n),cellx);
cellx = 1/(2*dx)*sparse(-diag(dVx)*tridiagx);
celly = 1/(2*dy)*sparse(-diag(dVy)*sparse(diag(-ones(n^2-n,1),-n)+ diag(ones(n^2-n,1),n)));
A=A+cellx+celly;
%A(1,1)=-A(1,2);
%A(end,end)=-A(end,end-1);

f=-ones(n^2,1);
u=A\f;
u=reshape(u,[n n]);
%u=[zeros(n+2,1) [zeros(1,n); u; zeros(1,n)] zeros(n+2,1)];
figure(2)
surf(X,Y,u)



