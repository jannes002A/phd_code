%Skript für ein 2D Potential mit BoltzmannReweigting
%clear all
%close all

%lambda=1:-.05:0.5;
%lambda=[0.02,0.01];
lambda=0.1;
N = 10000; %dt=1/N;
M = 1; %Anzahl der Trajektorieren

dt=0.1;
beta=4;
 
    
    xplot = linspace(-2,2);
    yplot = linspace(-2,2);

    [Xplot,Yplot] = meshgrid(xplot,yplot);
    
    xhplot = linspace(-2/lambda,2/lambda);
    yhplot = linspace(-2/lambda,2/lambda);

    [Xhplot,Yhplot] = meshgrid(xhplot,yhplot);



V = @(x,y) 3.*exp(-x.^2-(y-1/3).^2)-3.*exp(-x.^2-(y-5/3).^2) - 5.*exp(-(x-1).^2-y.^2)-5.*exp(-(x+1).^2-y.^2) + ...
    1/5.*x.^4 +1/5*(y-1/3).^4;

dxV =@(x,y) -6.*x.*exp(-x.^2-(y-1/3).^2)+6.*x.*exp(-x.^2-(y-5/3).^2) + 10.*(x-1).*exp(-(x-1).^2-y.^2)+10.*(x+1).*exp(-(x+1).^2-y.^2) + ...
    4/5.*x.^3;

dyV=@(x,y) -6.*(y-1/3).*exp(-x.^2-(y-1/3).^2)+6.*(y-5/3).*exp(-x.^2-(y-5/3).^2) +10.*y.*exp(-(x-1).^2-y.^2)+10.*y.*exp(-(x+1).^2-y.^2) + ...
     +4/5*(y-1/3).^3;


 
 
        HV=zeros(size(Xplot));

        for i=1:length(Xplot)
            for j=1:length(Yplot)
            HV(i,j)=HomotopieTest2D(Xplot(i,j)/lambda,Yplot(i,j)/lambda,lambda);
            end
        end

% simulation of two trajectories based with the same Brownian motion
% (Xtemp,Ytemp) is the trajectory in the smoothed potential
% (Xhtemp,Yhtemp) is the trajectory in the orignial potential

for i=1:M
   


        Xzero0 = -1;
        Yzero0 = 0;
        Xtemp = zeros(1,N);
        Ytemp = zeros(1,N);
        Xtemp(1) = Xzero0;
        Ytemp(1) = Yzero0;
        
        Xhzero0 = -1;
        Yhzero0 = 0;
        Xhtemp = zeros(1,N);
        Yhtemp = zeros(1,N);
        Xhtemp(1) = Xzero0;
        Yhtemp(1) = Yzero0;
    
     for j=2:N

      
     
        Btx = sqrt((2*dt)/beta)*randn(1);
        Bty = sqrt((2*dt)/beta)*randn(1);
        [~,dxf_l,dyf_l]=HomotopieTest2D(Xhtemp(j-1)/lambda,Yhtemp(j-1)/lambda,lambda);        
        Xhtemp(j) = Xhtemp(j-1) -(dxf_l)*dt +Btx;
        Yhtemp(j) = Yhtemp(j-1) -(dyf_l)*dt +Bty;
        
        %periodic boundary conditions on the square [-2,2]x[-2,2]
        if Xhtemp(j) > 2 
             Xhtemp(j)= -2 + abs(2-Xhtemp(j));
        end
        if Xhtemp(j) < -2
             Xhtemp(j) = 2 - abs(Xhtemp(j)+2);
        end
        if Yhtemp(j) > 2 
             Yhtemp(j)= -2 + abs(2-Yhtemp(j));
        end
        if Yhtemp(j) < -2
             Yhtemp(j) = 2 - abs(Yhtemp(j)+2);
        end
              
        
        Xtemp(j) = Xtemp(j-1) -dxV(Xtemp(j-1),Ytemp(j-1))*dt +Btx;
        Ytemp(j) = Ytemp(j-1) -dyV(Xtemp(j-1),Ytemp(j-1))*dt +Bty;
        

        %periodic boundary conditions on the square [-2,2]x[-2,2]
        if Xtemp(j) > 2 
             Xtemp(j)= -2 + abs(2-Xtemp(j));
        end
        if Xtemp(j) < -2
             Xtemp(j) = 2 - abs(Xtemp(j)+2);
        end
        if Ytemp(j) > 2 
             Ytemp(j)= -2 + abs(2-Ytemp(j));
        end
        if Ytemp(j) < -2
             Ytemp(j) = 2 - abs(Ytemp(j)+2);
        end

    end
 


end


    
         figure(1)
         subplot(2,3,1)
         surf(Xplot,Yplot,V(Xplot,Yplot))
         title('Orginal Potential')
         subplot(2,3,2)
         surf(Xplot,Yplot,exp(-beta.*V(Xplot,Yplot))./sum(sum(exp(-beta*V(Xplot,Yplot)))))
         title('Orginal Boltzmann')
         subplot(2,3,4)
         surf(Xplot,Yplot,HV)
         title('Smoothed Potential')
         subplot(2,3,5)
         surf(Xplot,Yplot,exp(-beta.*HV)./sum(sum(exp(-beta.*HV))))
         title('Smoothed Boltzmann')
    
         subplot(2,3,3)
         XY=[Xtemp' Ytemp'];
         plot(XY(:,1),XY(:,2),'.')
         title('Original Dynamcis')
    
         subplot(2,3,6)
         XYh=[Xhtemp' Yhtemp'];
         plot(XYh(:,1),XYh(:,2),'.')
         title('Smoothed Dynamics')
         %%
         
         figure(2)
         surf(Xplot,Yplot,HV-V(Xplot,Yplot))
    
        num_cluster=100;      
         
        [idx,C]=kmeans(XY,num_cluster);
        [idxh,Ch]=kmeans(XYh,num_cluster);
    
        figure(3)
        for i=1:num_cluster
            plot(XY(idx==i,1),XY(idx==i,2),'.','MarkerSize',12)
            plot(C(:,1),C(:,2),'kx','MarkerSize',1,'LineWidth',3); hold on
        end
        hold off
        
        figure(4)
        for i=1:num_cluster
            plot(XYh(idxh==i,1),XYh(idxh==i,2),'.','MarkerSize',12)
            plot(Ch(:,1),Ch(:,2),'kx','MarkerSize',1,'LineWidth',3); hold on
        end
        hold off
        
        C_min= min(C);
        C_max= max(C);
        
        Ch_min= min(Ch);
        Ch_max= max(Ch);
        
        cplot_min=min(C_min,Ch_min);
        cplot_max=max(C_max,Ch_max);
        
        cxplot= linspace(cplot_min(1),cplot_max(1));
        cyplot= linspace(cplot_min(2),cplot_max(2));
        
        [blzxplot,blzyplot] = meshgrid(cxplot,cyplot);
        
        Vplot=V(blzxplot,blzyplot);
        
        oboltz= exp(-beta.*Vplot)./sum(sum(exp(-beta*Vplot)));
        
        figure(5)
        surf(blzxplot,blzyplot,oboltz)
        title('Orginal Boltzmann')
        HVS=zeros(size(blzhxplot));
        
        for i=1:length(blzhxplot)
            for j=1:length(blzhyplot)
            HVS(i,j)=HomotopieTest2D(blzxplot(i,j)/lambda,blzyplot(i,j)/lambda,lambda);
            end
        end
        
        Diff=V(blzxplot,blzyplot)-HVS;
        pboltz = exp(-beta.*HVS)./sum(sum(exp(-beta*HVS)));
        cboltz = exp(-beta.*HVS).*exp(beta.*Diff)./sum(sum(exp(-beta*HVS).*exp(beta*Diff)));

        figure(6)
        surf(blzxplot,blzyplot,oboltz)
        title('Perturbed Boltzmann')
        
        figure(7)
        surf(blzxplot,blzyplot,pboltz)
        title('Corrected Boltzmann')
        
%         figure(8)
%         surf(blzhxplot,blzhyplot,exp(-beta.*(Vplot-HVS))./sum(sum(exp(-beta.*(Vplot-HVS)))))
%         title('Difference o.Boltzmann and p.Boltzmann')
%         
%         figure(9)
%         surf(blzhxplot,blzhyplot,exp(-beta.*(Vplot-(HVS+Diff)))./sum(sum(exp(-beta.*(Vplot-(HVS+Diff))))))
%         title('Difference o.Boltzmann and c.Boltzmann')
        %%
        figure(10)
        imagesc(exp(-beta.*Vplot))
        
        figure(11)
        imagesc(exp(-beta.*HVS))
        
        figure(12)
        imagesc(exp(-beta.*(HVS+Diff)))
        