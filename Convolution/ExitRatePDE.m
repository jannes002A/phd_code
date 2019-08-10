% skript for solving the exit rate for the 1d double well potential

poth =@(y) 8.*y.^4-44/3.*y.^3 + 2.*y.^2 +11/3.*y+1;
pothm =@(y,h) poth(y)+96.*h.^2+(4-88.*y+96.*y.^2).*h;
dpoth =@(y) 32.*y.^3-44.*y.^2+4.*y+11/3;
dpothm =@(y,t) dpoth(y)+(-88 +192.*y).*t^2/2;

t=0.0:.01:0.25;
%t=0;
dt=1/1000;
ext=zeros(length(t),1);

for i=1:length(t)



epsilon = sqrt(2*0.3);%2*gamma/sigma^2;


L_bound=-1;
R_bound=0.5;
n_boxes=1000;
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

%exit1=A1\b;
exit=A\b;

%exit time at -0.25
ext(i)=(exit(750)+exit(751))/2;

%Pdt=expm(A*dt);
%e=eigs(Pdt);
%figure(1)
%plot(1:length(e),e,'ro');
%figure(2)
%plot(x,poth(x))

end


save('PDETimes.mat','exit','ext','-v7.3')

%%
% figure(1)
% plot(x,pothm(x,t))
% 
% % figure(2)
% % plot(x,exit)
% 
% figure(3)
% plot(x,exit)

figure(1)
plot(t(1:end),ext(1:end),'LineWidth',4); hold on
plot(t,1./f3(t),'LineWidth',4)
%plot(t,1./f2(t),'LineWidth',4)
hold off
xl=xlabel('Lambda');
yl=ylabel('Average Exit Time for x=-0.25');
title('Stopping times')
l= legend('Exact','Extra','Location','northeast');
set(xl,'Interpreter','Latex');
set(yl,'Interpreter','Latex');
set(l,'Interpreter','Latex');

figure(2)
semilogy(t(1:end),1./ext,'LineWidth',4);hold on
semilogy(t(1:end),f3(t(1:end)),'g','LineWidth',4);
%semilogy(t(1:end),f2(t(1:end)),'k','LineWidth',4);
semilogy(lambda,rates','rx','LineWidth',4);
%semilogy(t(1:end),abs(f2(t)),'m');
hold off
title('Rates')
xlabel('Lambda')
ylabel('Average Exit Rate for x=-0.25')
l= legend('Exact','Extra','Sampling','Location','northwest');
set(l,'Interpreter','Latex');

figure(3)
plot(t(1:end),1./ext,'LineWidth',4);hold on
plot(t,f3(t),'LineWidth',4)
%plot(t,f2(t),'LineWidth',4)
plot(lambda,rates','rx','LineWidth',4);
hold off
title('Rates')
xl=xlabel('Lambda');
yl=ylabel('Average Exit Rate for x=-0.25');
set(xl,'Interpreter','Latex');
set(yl,'Interpreter','Latex');
l= legend('Exact','Extra','Sampling','Location','northwest');
set(l,'Interpreter','Latex');
% figure(1)
% plot(x,exit)
% 
% figure(2)
% plot(x(100:end-100),exit(100:end-100))
% 
% figure()
% plot(x(100:end-100),1./exit(100:end-100))