clear all 
%dataset below contains:
% -inf. gen A
% -grid x with x(2)-x(1)=0.0040
% -inv. temp. BETA=4 and
% -Pot type 'POT'='asym2wellPot'
% -Potential energy function values 'V', defined on x
load asym2wellPot.mat A BETA POT 
% The inf. gen A is sparse, so make A full for later computations
AA=full(A);
% Store \varepsilon=BETA^{-1}
epsilon=1.0/BETA;
% J: x=-3:3:1500 V=V(x) 
load DWPOT.mat V x 

% I contains all the indices of the points in x
I=1:1:length(x);

% I1 contains the indices of points in x belonging to S:=[-1.1,-1]
I1  = find((x>=-1.1)&(x<=-1));

% J ctns indices of points in x not in S
J=zeros(1,length(I)-length(I1));
j=0;
for i=1:length(x)
    if (x(i)<-1.1)||(x(i)>-1)
            j=j+1;
            J(j)=i;
    end
end

% K contains the indices of points in x within the disk of radius 2
K=find(abs(x)<2);

% L is the inf. gen restricted to points not in [-1.1,-1]
L = AA(J,J);

%% Compute and Display Psi(x)=E_x[exp(-\tau / \varepsilon)]

%Solve the PDE \varepsilon L \psi (x) - \sigma \psi (x) = 0
% using the finite element method. Do this in a few stages.

%Stage 0: Define the discretised solution \psi to be a vector of 1's
%         This follows since for x in S, \tau(x)=0 and hence \psi(x)=1
Psi = ones(size(I))'; 

%Stage 1: Compute \psi(x) for points x which are not in S
%         Store these values in PsiJ
tempmat=AA(J,I1)';
zx=sum(tempmat)';
%zx=zeros(length(L),1);
PsiJ = -(epsilon*L-eye(size(L)))\(epsilon*zx);
%Note: it is not clear why the right hand side is not just a vector of 0's

%Stage 2: Update the discretised solution \psi with PsiJ
Psi(J) = PsiJ;

%% Compute the mean first hitting time E_x[\tau] 
dPsi_dsigma=zeros(length(K),1);

%Apply centered finite differences to psi
sigma=0.0000001;
PsiJ_sigmaplus=-(epsilon*L-sigma*eye(size(L)))\(epsilon*zx);
PsiJ_sigmaminus=-(epsilon*L+sigma*eye(size(L)))\(epsilon*zx);

Psi_sigmaplus=ones(size(I));
Psi_sigmaminus=ones(size(I));

Psi_sigmaplus(J)=PsiJ_sigmaplus;
Psi_sigmaminus(J)=PsiJ_sigmaminus;

MFHT=-epsilon*(Psi_sigmaplus-Psi_sigmaminus)/(2*sigma);

%% Compute the cumulant generating function F(x)=-\varepsilon log \psi (x)
F = -epsilon * log(Psi);

%% Compute the mean first hitting time E_x[\tau] 
F_sigmaplus=-epsilon*log(Psi_sigmaplus);
F_sigmaminus=-epsilon*log(Psi_sigmaminus);
MFHT_=(F_sigmaplus-F_sigmaminus)/(2*sigma);

%% Compute the optimal feedback control c^\ast(x)=-\sqrt{2} \nabla F(x)

%Define the vectors \nabla F and c inside the disk of radius 2
gradF=zeros(length(K),1);
c=zeros(size(gradF));
%Compute gradient of F using centered finite differences
for i=1:length(K)
    gradF(i)=(F((K(i)+1))-F(K(i)-1))/(x(K(i)+1)-x(K(i)-1));
    c(i)=-sqrt(2)*gradF(i);
end


%% Do plots

%Plot psi
h1=figure(1);
clf
plot(x,Psi)
xlabel('x')
ylabel('\psi(x)')
title('Plot of \psi(x) = E_x [ exp( - \tau / \epsilon ) ]')
%title_root_=strcat(title_root,'_mass');
%title_fig1='Figures/psiplot';
%print(h1,'-dpng',title_fig1)
%saveas(gcf,strcat(title_fig1,'.fig'));
%print(h1,'-depsc',title_fig1)


%Plot MFHT computed by taking centered finite differences
h2=figure(2);
clf
plot(x,MFHT,'-r',x,MFHT_,'--b')
h = legend('MHFT_\psi','MFHT_F',1);
xlabel('x')
ylabel(' E_x[\tau] ' )
title('Plot of mean first hitting time E_x [ \tau ] ')
%title_fig2='Figures/mfhtplot';
%print(h2,'-dpng',title_fig2)
%saveas(gcf,strcat(title_fig2,'.fig'));
%print(h2,'-depsc',title_fig2)

% KK is the indicator function for (-1,2)
% We choose this set because we want to know the control policy for points
% that will not directly lead to S
KK=find((x>-1)&(x<2));

KK_=find((x(K)>-1)&(x(K)<2));

%Plot cumulant generating function (CGF) and gradient of CGF for (-1,2)
h3=figure(3);
clf
plot(x(KK),F(KK),'-r',x(KK),c(KK_),'-b',x(KK),V(KK),'-k')
h = legend('F(x)','c^\ast(x)','V(x)',1);
xlabel('x')
title('Cmlt gen fn F(x), opt. control c^\ast(x)=-sqrt(2)F(x), pot. energy fn V(x)')
%title_fig3='Figures/cgf_control_poteng_plot';
%print(h3,'-dpng',title_fig3)
%saveas(gcf,strcat(title_fig3,'.fig'));
%print(h3,'-depsc',title_fig3)

%Plot original, optimally tilted pot. energies V(x) and G(x)=V(x)+2F(x)
h4=figure(4);
clf
plot(x(K),V(K),'-r',x(K),(V(K)+2*F(K)),'-b')
h = legend('V(x)','V(x)+2F(x)',1);
xlabel('x')
title('Potential energy landscapes')
%title_fig4='Figures/potengfuns_plot';
%print(h4,'-dpng',title_fig4)
%saveas(gcf,strcat(title_fig4,'.fig'));
%print(h4,'-depsc',title_fig4)

%% Save all the data
%save MFHT_PDEFiniteElement.mat
