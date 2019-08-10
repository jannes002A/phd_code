function [cV,cdx,cdy]=ConBrownianMotion(x,y,lambda,V,dx,dy)

n=20;

xapp = x + lambda.*randn(n,1);
yapp = y + lambda.*randn(n,1);

[xh,yh]=meshgrid(xapp,yapp);

cV =  sum(sum(V(xh,yh)))/n^2;
cdx = sum(sum(dx(xh,yh)))/n^2;
cdy = sum(sum(dy(xh,yh)))/n^2;