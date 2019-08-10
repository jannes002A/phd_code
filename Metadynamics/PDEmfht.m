% skript for solving the exit rate for the 1d double well potential

poth =@(y)  1/2.*(y.^2-1).^2;
dpoth =@(y) 2.*y.*(y.^2-1);



beta=3;


epsilon = sqrt(2/beta);%2*gamma/sigma^2;


L_bound=-1.5;
R_bound=1;
n_boxes=1000;
dx=(R_bound-L_bound)/(n_boxes-1);
x=L_bound:dx:R_bound;
e=ones(n_boxes,1);

A = spdiags(epsilon^2/2*e*[1 -2 1]/(dx^2),[-1 0 1],n_boxes,n_boxes) + ...
    spdiags(-1*dpoth(x)'*[-1 0 1]/(2*dx),[1 0 -1],n_boxes,n_boxes)';
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
%ext(i)=(exit(750)+exit(751))/2;


plot(x,exit)


%save('PDETimes.mat','exit','ext','-v7.3')

%%
% figure(1)
% plot(x,pothm(x,t))
% 
% % figure(2)
% % plot(x,exit)
% 
% figure(3)
% plot(x,exit)

% figure(1)
% plot(t(1:end),ext(1:end),'LineWidth',4); hold on
% plot(t,1./f3(t),'LineWidth',4)
% %plot(t,1./f2(t),'LineWidth',4)
% hold off
% xl=xlabel('Lambda');
% yl=ylabel('Average Exit Time for x=-0.25');
% title('Stopping times')
% l= legend('Exact','Extra','Location','northeast');
% set(xl,'Interpreter','Latex');
% set(yl,'Interpreter','Latex');
% set(l,'Interpreter','Latex');
% 
% figure(2)
% semilogy(t(1:end),1./ext,'LineWidth',4);hold on
% semilogy(t(1:end),f3(t(1:end)),'g','LineWidth',4);
% %semilogy(t(1:end),f2(t(1:end)),'k','LineWidth',4);
% semilogy(lambda,rates','rx','LineWidth',4);
% %semilogy(t(1:end),abs(f2(t)),'m');
% hold off
% title('Rates')
% xlabel('Lambda')
% ylabel('Average Exit Rate for x=-0.25')
% l= legend('Exact','Extra','Sampling','Location','northwest');
% set(l,'Interpreter','Latex');
% 
% figure(3)
% plot(t(1:end),1./ext,'LineWidth',4);hold on
% plot(t,f3(t),'LineWidth',4)
% %plot(t,f2(t),'LineWidth',4)
% plot(lambda,rates','rx','LineWidth',4);
% hold off
% title('Rates')
% xl=xlabel('Lambda');
% yl=ylabel('Average Exit Rate for x=-0.25');
% set(xl,'Interpreter','Latex');
% set(yl,'Interpreter','Latex');
% l= legend('Exact','Extra','Sampling','Location','northwest');
% set(l,'Interpreter','Latex');
% figure(1)
% plot(x,exit)
% 
% figure(2)
% plot(x(100:end-100),exit(100:end-100))
% 
% figure()
% plot(x(100:end-100),1./exit(100:end-100))