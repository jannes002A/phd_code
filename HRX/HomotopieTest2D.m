% Homotopie für das 2D Potential mit einer Gauß Hermit Formel

function [fx,dfx,dfy]=HomotopieTest2D(x,y,lambda)

% wie müssen die Stützstellen gewählt werden?

w = [0.01995324,0.39361932,0.94530872,0.39361932,0.01995324];
xh= (x +[ -2.02018287,-0.95857246,0,0.95857246,2.02018287]);
yh= (y +[ -2.02018287,-0.95857246,0,0.95857246,2.02018287]);

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



fx=0;
dfx=0;
dfy=0;

for j=1:5
    wj=w(j);
    for i=1:5
       wi=w(i); 
        fx  = fx  + wj*wi*V1(lambda*xh(i),lambda*yh(j));
        dfx = dfx + wj*wi*(lambda*dxV1(lambda*xh(i),lambda*yh(j)));
        dfy = dfy + wj*wi*(lambda*dyV1(lambda*xh(i),lambda*yh(j)));
    end
end

fx  = 1/pi * fx;
dfx = 1/pi *dfx;
dfy = 1/pi *dfy;

%dfx = dfx + wj*wi*(-2*xh(i)*V1(lambda*xh(i),lambda*yh(j)))+lambda*dxV1(lambda*xh(i),lambda*yh(j));
%dfy = dfy + wj*wi*(-2*yh(i)*V1(lambda*xh(i),lambda*yh(j)))+lambda*dyV1(lambda*xh(i),lambda*yh(j));