%Skript für ein 2D Potential mit HRX (homotopy replica exchange)
%clear all
%close all

%lambda=1:-.05:0.5;
lambda=0.7;
exit = zeros(1,length(lambda));
beta=10;

V = @(x,y) 3.*exp(-x.^2-(y-1/3).^2)-3.*exp(-x.^2-(y-5/3).^2) - 5.*exp(-(x-1).^2-y.^2)-5.*exp(-(x+1).^2-y.^2) + ...
    1/5.*x.^4 +1/5*(y-1/3).^4;

dxV =@(x,y) -6.*x.*exp(-x.^2-(y-1/3).^2)+6.*x.*exp(-x.^2-(y-5/3).^2) + 10.*(x-1).*exp(-(x-1).^2-y.^2)+10.*(x+1).*exp(-(x+1).^2-y.^2) + ...
    4/5.*x.^4;

dyV=@(x,y) -6.*(y-1/3).*exp(-x.^2-(y-1/3).^2)+6.*(y-5/3).*exp(-x.^2-(y-5/3).^2) +10.*y.*exp(-(x-1).^2-y.^2)+10.*y.*exp(-(x+1).^2-y.^2) + ...
     +4/5*(y-1/3).^3;

N = 10000; %dt=1/N;
M = 1; %Anzahl der Trajektorieren

dt=0.1;
count=zeros(1,N);
swap=0;

for i=1:M


        Xzerox1 = 0;
        Xzeroy1 = 0;
        Xtempx1 = Xzerox1;
        Xtempy1 = Xzeroy1;
        Xemx1 = zeros(1,N);
        Xemy1 = zeros(1,N);


        Xzerox2 = 0;
        Xzeroy2 = 0;
        Xtempx2 = Xzerox1;
        Xtempy2 = Xzeroy1;
        Xemx2 = zeros(1,N);
        Xemy2 = zeros(1,N);


 
    for j=1:N
    
        %time evolution of the SDE in a 2D potential
        Btx1 = sqrt(2*dt/beta)*randn;
        Bty1 = sqrt(2*dt/beta)*randn;
        V1=V(Xtempx1,Xtempy1);
        Xtempx1 = Xtempx1 -(dxV(Xtempx1,Xtempy1))*dt +Btx1;
        Xtempy1 = Xtempy1 -(dyV(Xtempx1,Xtempy1))*dt +Bty1;
        
        if Xtempx1 > 2 
             Xtempx1 = -2 + abs(2-Xhtempx1);
        end
        if Xtempx1 < -2
             Xtempx1 = 2 - abs(Xtempx1+2);
        end
        if Xtempy1 > 2 
             Xtempy1= -2 + abs(2-Xtempy1);
        end
        if Xtempy1 < -2
             Xtempy1 = 2 - abs(Xtempy1+2);
        end
     
        Btx2 = sqrt(2*dt/beta)*randn;
        Bty2 = sqrt(2*dt/beta)*randn;
        [V2,dxf_l,dyf_l]=ConBrownianMotion(Xtempx2,Xtempy2,lambda,V,dxV,dyV);
        Xtempx2 = Xtempx2 -(dxf_l)*dt +Btx1;
        Xtempy2 = Xtempy2 -(dyf_l)*dt +Bty1;
        
        if Xtempx2 > 2 
             Xtempx2= -2 + abs(2-Xtempx2);
        end
        if Xtempx2 < -2
             Xtempx2 = 2 - abs(Xtempx2+2);
        end
        if Xtempy2 > 2 
             Xtempy2= -2 + abs(2-Xtempy2);
        end
        if Xtempy2 < -2
             Xtempy2 = 2 - abs(Xtempy2+2);
        end
        
        
        if swap==0
            Xemx1(j) = Xtempx1;
            Xemy1(j) = Xtempy1;
            
            Xemx2(j) = Xtempx2;
            Xemy2(j) = Xtempy2;
        else
            Xemx1(j) = Xtempx2;
            Xemy1(j) = Xtempy2;
            
            Xemx2(j) = Xtempx1;
            Xemy2(j) = Xtempy1;
        end
        
        
        V21=ConBrownianMotion(Xtempx1,Xtempy1,lambda,V,dxV,dyV);
        V12=V(Xtempx2,Xtempy2);
        
        %Metroplis Hastings criterium to check for possible swapping
        alpha = min(1,((exp(-beta *V1)*exp(-beta*V2))/...
         (exp(-beta *V21)*exp(-beta*V12)))) ;
        
        %choose random number
         u=rand;
     
        count(j)= (u< alpha);
        
        % swapping the coordiantes if the energy allows it 
        if (u<alpha)
            Xhx=Xtempx1;
            Xhy=Xtempy1;
            Xtempx1=Xtempx2;
            Xtempy1=Xtempy2;
            Xtempx2=Xhx;
            Xtempy2=Xhy;
            swap=1;
        else
            swap=0;
        end
        
        
     if (isnan(Xemx1(j)) || isnan(Xemy1(j)))

         disp(j)
         return;
     end
     
    end
end

accrate= (sum(count)/N)*100;
fprintf('AccRate %2.4f \n',accrate)
 

xplot = linspace(-2,2);
yplot = linspace(-2,2);

[Xplot,Yplot] = meshgrid(xplot,yplot);

figure(1)
surf(Xplot,Yplot,V(Xplot,Yplot)), hold on
plot3([Xzerox1 Xemx1],[Xzeroy1 Xemy1],V([Xzerox1 Xemx1],[Xzeroy1 Xemy1]),'r+'),%hold off
hold off
view(85,70)



