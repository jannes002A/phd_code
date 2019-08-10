%Skipt zur Simulation mit Homotopie
clear Xem
close all
% %Range
x=linspace(-0.5,1.5);
%e.g. 
poth =@(y) 8.*y.^4-44/3.*y.^3 + 2.*y.^2 +11/3.*y+1;
dpoth =@(y) 32.*y.^3-44.*y.^2+4.*y.^2+11/3;
pot = poth(x);
dpot= dpoth(x);
t=0.08;
hpoth =@(y) poth(y)+96.*t.^2+(4-88.*y+96.*y.^2).*t;
hpot=hpoth(x);

N=100000;
Xzero=-0.25;
dt=1/1000;

beta = 3;

 
    
    dhpoth =@(y) dpoth(y)+(-88 +192.*y).*t; 
    
    %homopoth =@(y) poth(y)+24.*t(i).^4+(2-44.*y+48.*y.^2).*t(i).^2;
    %dhomopoth =@(y) dpoth(y)+(-44 +48.*y).*t(i).^2; 
    
     Xtemp=zeros(1,N);
     Xtemp(1)=Xzero;
     
     
    for j=2:N
       Xtemp(j) = Xtemp(j-1) -(dhpoth(Xtemp(j-1))).*dt + sqrt(2*dt/beta)*randn;
    end
    
    %diff = hpoth-poth
    diff=@(y) 96.*t.^2+(4-88.*y+96.*y.^2).*t;
    sdiff= sum(diff(x));
%     
%     for j=1:N
%         Xpert (j) = (Xtemp(j)*exp(-beta*diff(Xtemp(j))))/sdiff;
%     end
    Xhelp=zeros(1,length(Xtemp));
    Xhelp(1)=Xzero;
    Xh=Xzero;
    j=1;
    %correct Normalization
    % Backwards jumps had to be excluded
    for i=1:length(Xtemp)
        if (Xh < Xtemp(i))
            Xhelp(i)= Xh;
            Xh=Xtemp(i);
            j=j+1;
        end
    end
    
    
    orgBoltz=exp(-beta*poth(Xtemp))./sum(exp(-beta*poth(Xtemp)));
    pertBoltz = exp(-beta*hpoth(Xtemp))./sum(exp(-beta*hpoth(Xtemp)));
    corrBoltz = (exp(-beta*(hpoth(Xtemp))).*exp(beta*diff(Xtemp)) )./sum(exp(-beta*(hpoth(Xtemp))).*exp(beta*diff(Xtemp)));
    
    subplot(2,3,1)
    plot(Xtemp,orgBoltz); 
    title('Orgiginal Boltzmann Sampled')
    %hold on
    %hold off
    
    subplot(2,3,2)
    plot(Xtemp,pertBoltz)
    title(['Perturbed Boltzmann Sampled t= ',num2str(t)])
    
    subplot(2,3,3)
    plot(Xtemp,corrBoltz);
    title('Corrected Boltzmann Sampled')
    % wenn man die richtige Konstante nimmt dann wird auch die richtige
    % Boltzmannverteilung gesampelt sum(exp(-beta*(hpoth(x))).*exp(beta*diff(x))
    
    subplot(2,3,4)
    plot(x,exp(-beta*poth(x))./sum(exp(-beta*poth(x))))
    title('Orgiginal Boltzmann')
    
    subplot(2,3,5)
    plot(x,exp(-beta*hpoth(x))./sum(exp(-beta*hpoth(x))))
    title('Perturbed Boltzmann')
    
%     subplot(2,3,6)
%     plot(x,exp(-beta*(poth(x)))./sum(exp(-beta*(poth(x))))-exp(-beta*(hpoth(x)))./sum(exp(-beta*(hpoth(x)))))
%     title('Difference Distribution')
    
    subplot(2,3,6)
    plot(Xtemp,orgBoltz-corrBoltz)
    title('Difference Sampled')
   
    