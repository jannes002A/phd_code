%Skript für ein 2D Potential mit Homotopie
%clear all
%close all

%lambda=1:-.05:0.5;
lambda=[0.02,0.01];
exit = zeros(1,length(lambda));
for k=1:length(lambda)
    ex=0;
    k
    
    xplot = linspace(-2,2);
    yplot = linspace(-2,2);

    [Xplot,Yplot] = meshgrid(xplot,yplot);


beta=4;

    xplotHom = linspace(-2/lambda(k),2/lambda(k));
    yplotHom = linspace(-2/lambda(k),2/lambda(k));

    [XplotHom,YplotHom] = meshgrid(xplotHom,yplotHom);

V = @(x,y) 3.*exp(-x.^2-(y-1/3).^2)-3.*exp(-x.^2-(y-5/3).^2) - 5.*exp(-(x-1).^2-y.^2)-5.*exp(-(x+1).^2-y.^2) + ...
    1/5.*x.^4 +1/5*(y-1/3).^4;

dxV =@(x,y) -6.*x.*exp(-x.^2-(y-1/3).^2)+6.*x.*exp(-x.^2-(y-5/3).^2) + 10.*(x-1).*exp(-(x-1).^2-y.^2)+10.*(x+1).*exp(-(x+1).^2-y.^2) + ...
    4/5.*x.^3;

dyV=@(x,y) -6.*(y-1/3).*exp(-x.^2-(y-1/3).^2)+6.*(y-5/3).*exp(-x.^2-(y-5/3).^2) +10.*y.*exp(-(x-1).^2-y.^2)+10.*y.*exp(-(x+1).^2-y.^2) + ...
     +4/5*(y-1/3).^3;


xc=0;
yc=0;
sigma=0.01;
E=@(x,y) 1 / (2 * sqrt(2*pi))* exp(- ((x-xc).^2 + (y-yc).^2)./(2*sigma^2));

HV = conv2(V(Xplot,Yplot),E(Xplot,Yplot));


        HV=zeros(size(Xplot));

        for i=1:length(Xplot)
            for j=1:length(Yplot)
            HV(i,j)=HomotopieTest2D(XplotHom(i,j),YplotHom(i,j),0.025);
            end
        end


        figure(1)
        surf(Xplot,Yplot,HV)
        title('Homotopie')
        figure(2)
        surf(Xplot,Yplot,V(Xplot,Yplot))
        title('Orginal')
        figure(3)
        %imagesc(abs(V(Xplot,Yplot)-HV))
        surf(Xplot,Yplot,V(Xplot,Yplot)-HV)
        title('Differenzbild')
%%
N = 1000; %dt=1/N;
M = 1000; %Anzahl der Trajektorieren

dt=1/2^11;
%epsilon=sqrt(dt/beta);


%dW=sqrt(dt)*randn(1,N);
%Xem = zeros(1,N);
for i=1:M
   
%         Xzerox = -1;
%         Xzeroy = 0;
%         Xtempx = Xzerox;
%         Xtempy = Xzeroy;
%         Xemx = zeros(1,N);
%         Xemy = zeros(1,N);

        Xzerox1 = -1/lambda(k);
        Xzeroy1 = 0/lambda(k);
        Xtempx1 = Xzerox1;
        Xtempy1 = Xzeroy1;
        Xemx1 = zeros(1,N);
        Xemy1 = zeros(1,N);
% 
%         dWx = sqrt(dt)*randn(1,N);
%         dWy = sqrt(dt)*randn(1,N);

        dW1x = sqrt(dt)*randn(1,N);
        dW1y = sqrt(dt)*randn(1,N);
% 
%         Btx=zeros(1,N);
%         Bty=zeros(1,N);
        Btx1=zeros(1,N);
        Bty1=zeros(1,N);

 
    for j=1:N
        %Winc = sum(dW((j-1)+1:j));
        %Xtemp = Xtemp -(dpoth(Xtemp) + dbasisfun(Xtemp))*dt + sqrt(dt)*epsilon*randn; 
%         Btx(j) = sum(dWx((j-1)+1:j));
%         Bty(j) = sum(dWy((j-1)+1:j));
%         Xtempx = Xtempx -(dxV(Xtempx,Xtempy))*dt +Btx(j);% epsilon*randn;
%         Xtempy = Xtempy -(dyV(Xtempx,Xtempy))*dt +Bty(j);% epsilon*randn;
%         Xemx(j) = Xtempx;
%         Xemy(j) = Xtempy;
%      
     
        Btx1(j) = sum(dW1x((j-1)+1:j));
        Bty1(j) = sum(dW1y((j-1)+1:j));
        [~,dxf_l,dyf_l]=HomotopieTest2D(Xtempx1,Xtempy1,lambda(k));
        Xtempx1 = Xtempx1 -(dxf_l)*dt +Btx1(j);% epsilon*randn;
        Xtempy1 = Xtempy1 -(dyf_l)*dt +Bty1(j);% epsilon*randn;
        Xemx1(j) = Xtempx1;
        Xemy1(j) = Xtempy1;
     
     %if (isnan(Xemx(j)) || isnan(Xemy(j)) || isnan(Xemx1(j)) || isnan(Xemy1(j)))
%          isnan(Xemx(j)) 
%          isnan(Xemy(j)) 
%          isnan(Xemx1(j)) 
%          isnan(Xemy1(j))
     if (isnan(Xemx1(j)) || isnan(Xemy1(j)))

         j
         return;
     end
     if (Xemx1(j) > 0 && Xemy1(j) > 1)
        ex=ex+1;
        break;
     end
    end
 end
 
 exit(k)=ex/M;

figure(4)
surf(Xplot,Yplot,V(Xplot,Yplot)), hold on
plot3([Xzerox Xemx],[Xzeroy Xemy],V([Xzerox Xemx],[Xzeroy Xemy]),'r+'),%hold off
hold off
legend('Pot','Pot(X_{SDE}))')
title('Original Potential')
view(85,70)



HV=zeros(size(XplotHom));

for i=1:length(XplotHom)
    for j=1:length(YplotHom)
    HV(i,j)=HomotopieTest2D(XplotHom(i,j),YplotHom(i,j),lambda(k));
    end
end

HVP=zeros(1,N+1);
HVP(1)  =   HomotopieTest2D(Xzerox1,Xzeroy1,lambda(k));
for i=2:N+1
    HVP(i) = HomotopieTest2D(Xemx1(i-1),Xemy1(i-1),lambda(k));
end

figure(5)
surf(XplotHom,YplotHom,HV), hold on
plot3([Xzerox1 Xemx1],[Xzeroy1 Xemy1],HVP,'r+'),%hold off
hold off
legend('Pot','Pot(X_{SDE}))')
title('Homotopie des 2D Potentials')
view(85,70)

end