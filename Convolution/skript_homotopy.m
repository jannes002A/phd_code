clear;

%Temperatur
beta=3; beta2=0.1;

P=[];
Y=[];
V2=[]; U1=[]; VV=[];S=[];
%Homotopieparameter
r=0:0.0075:0.15;
%%r=0.15;
%r=1;
x=-0.5:0.02:1.5;
m=length(x);
%Pentential
f=8*x.^4-44/3*x.^3+2*x.^2+11/3*x+1;

% Ft=f + 96*r^2 + r*(4-88*x+96*x.^2);
% %p = exp(-beta*Ft); p=p/sum(p);
% 
% figure(1)
% plot(x,Ft)
% hold on

for t=r
   % Potential in der Homotopie
   Ft=f + 96*t^2 + t*(4-88*x+96*x.^2);
   % Bolzmann Verteilung des Potentials
   p = exp(-beta*Ft); p=p/sum(p);
   P=[P;p];
   % Entropie
   s=-sum(p.*log(p));
   S=[S,s];
   % Infinitissimaler Generator in einer Gitterdiskretisierung 
   % basierend auf den Bolzmanngewichten 
%    Q=zeros(m,m);
%    Q(1,2)=1;
%    for i=2:m-1
%        Q(i,i+1)=1;
%        Q(i,i-1)=1;
%    end
   q=ones(1,m-1);
   Q=diag(q,1)+diag(q,-1);
%  Q(m,m-1)=1;
   Q=diag(1./sqrt(p))*Q*diag(sqrt(p));
   Q=Q-diag(sum(Q,2));
   % Linke Eigenwerte des infinitissimalen approximierten Generators
   [evec,eval]=eigs(Q,3,'lr');
   Y=[Y,diag(eval)];
   % Weiterarbeiten mit den Eigenvektoren
   v2=evec(:,2)*sign(evec(1,2));
   % Comittorfunktionen ? oder f(X_t)??
   chi=(v2-min(v2))/(max(v2)-min(v2));
   V2=[V2,chi];
   % zu lösende Eigenwertgleichung ?
   Qm=Q-beta2*diag(chi);
   [psi,v]=eigs(Qm,1,'lr');
   % 
   psi=psi*sign(psi(floor(m*2/3)));
   psi=psi/max(psi);
   U1=[U1,psi];
   VV=[VV,v];
end
P=P';

ss=[9,12,15,18,21];
% plot des Potentials mit unterscheidlichen Homotopieparametern
figure(1); plot(x,P(:,[1,ss]));
% Speed up
figure(2); plot(r([1,ss]),Y(2,[1,ss])/Y(2,1),'-*');
% Zugehörigkeitsfunktion
figure(3); plot(x,V2);
% Zugehörigkeit versus Verweildauer (?)
figure(4); plot(x,[V2(:,[1,ss]),U1(:,[1,ss])]);
% Homotopie versus Entropie
figure(5); plot(r,S,'-*');

y=log(-VV)'; A=[ones(21,1),r',(r.^2)', (r.^3)', (r.^4)'];
%Monroe Inverse
ww=pinv(A(ss,:))*y(ss);
% Umgehen der Inversen
wwJ=A(ss,:)\y(ss);
% Approximation der Übergangswahrscheinlichkeit
figure(6); plot(r,[A*ww,y],'-*');
figure(7); plot(r,[A*wwJ,y],'-*');