% skript for solving the exit rate for the 1d double well potential

poth =@(y) 8.*y.^4-44/3.*y.^3 + 2.*y.^2 +11/3.*y+1;
%pothm =@(y,h) poth(y)+96.*h.^2+(4-88.*y+96.*y.^2).*h;
hpoth =@(y,t) poth(y)+96.*(t^2/2).^2+(4-88.*y+96.*y.^2).*(t^2/2);
dpoth =@(y) 32.*y.^3-44.*y.^2+4.*y+11/3;
dpothm =@(y,t) dpoth(y)+(-88 +192.*y).*t^2/2;



t=0:0.1:0.4;

dt=1/1000;
ext=zeros(length(t),1);

ei=zeros(6,length(t));

for i=1:length(t)

epsilon = sqrt(2*0.4);%2*gamma/sigma^2;


L_bound=-1;
R_bound=2;
n_boxes=500;
dx=(R_bound-L_bound)/(n_boxes-1);
x=L_bound:dx:R_bound;
e=ones(n_boxes,1);

A = spdiags(epsilon^2/2*e*[1 -2 1]/(dx^2),[-1 0 1],n_boxes,n_boxes) + ...
    spdiags(-1*dpothm(x,t(i))'*[-1 0 1]/(2*dx),[1 0 -1],n_boxes,n_boxes)';
    % diag(-1*dpothm(x,t))*spdiags(e*[-1 0 1]/(2*dx),[-1 0 1],n_boxes,n_boxes)
%Neumann Randbedingungen auf der rechten Seite und Dirichlet auf der linken Seite     
%A(1,1) = -A(1,2); 
%A(n_boxes,n_boxes) = -A(n_boxes,n_boxes-1);


%  A1 = spdiags(sigma^2/(2*gamma^2)*e*[1 -2 1]/(dx^2),[-1 0 1],n_boxes,n_boxes) + ...
%       sparse(diag(-1*dpothm(x,t))*spdiags(e*[-1 0 1]/(2*dx),[-1 0 1],n_boxes,n_boxes));
     % so implementiert wird falsch mit nabla V multipliziert
     % spdiags(-1*dpothm(x,t)'*[-1 0 1]/(2*dx),[-1 0 1],n_boxes,n_boxes);

b=-ones(n_boxes,1);
%b(end)=-1-1/dx^2;
%b(1)=b(1)-1/dx^2;
%b(end)=b(end)-1/dx^2;


Pdt=expm(A*dt);
ei(:,i)=eigs(Pdt);

% figure(i)
% subplot(1,2,1)
% plot(1:length(ei(:,i)),ei(:,i),'ro');
% subplot(1,2,2)
% plot(x,hpoth(x,t(i)))
% axis([-1 2 -1 30])
end

figure()
plot(ei,'o-','MarkerSize',7)
xlabel('Component')
ylabel('Eigenvalue')
title('Eigenvalues Transferoperator')
l=legend('$\lambda=0$','$\lambda=0.1$','$\lambda=0.2$','$\lambda=0.3$','$\lambda=0.4$');
set(l,'Interpreter','Latex');  