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
t=0.06;
hpoth =@(y) poth(y)+96.*t.^2+(4-88.*y+96.*y.^2).*t;
hpot=hpoth(x);
N=100000;
x=linspace(-0.5,1.5);

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
%     for i=1:length(Xtemp)
%         if (Xh < Xtemp(i))
%             Xhelp(i)= Xh;
%             Xh=Xtemp(i);
%             j=j+1;
%         end
%     end
    
    [N,C]=hist(Xtemp,45);
    bin=zeros(1,length(Xtemp));
    binh=zeros(1,length(Xtemp));
    
    h=C(1)-C(2);
    C= [C(1)+h,C,C(end)-h];
    N= [0,N,0];
    for i=1:length(Xtemp)
        ch= abs(C-Xtemp(i));
        [~,bin(i)]=min(ch);
    end
    binval= exp(-beta.*hpoth(Xtemp));
    rew = exp(beta.*(diff(Xtemp)))./sum(exp(beta.*diff(Xtemp)).*exp(-beta.*poth(Xtemp)));
    %
    histval=zeros(1,length(C));
    for i=1:length(C)
        histval(i) = sum(rew(bin==i));
        binh(bin==i) = histval(i);
    end
    figure(1)
    plot(C,histval)
    
    figure(2)
    plot(Xtemp,exp(-beta.*poth(Xtemp))./sum(exp(-beta.*poth(x))))
    
    rew2 = exp(beta.*(diff(C))).*exp(-beta.*hpoth(C));
    
%     figure(2)
%     hist(Xtemp,50)
    
    figure(3)
    stem(C,rew2); 
    figure(4)
    stem(C,exp(-beta.*poth(C)));
    figure(5)
    stem(C,exp(-beta.*hpoth(C)))
    
%     figure(4)
%     plot(Xtemp,exp(-beta.*poth(Xtemp))./sum(exp(-beta.*poth(Xtemp))))
%     
%     figure(5)
%     plot(x,exp(-beta.*poth(x))./sum(exp(-beta.*poth(x))))
    
%     orgBoltz = exp(-beta*poth(C))./N.*sum(exp(-beta*poth(C)));
%     pertBoltz = exp(-beta*hpoth(C))./N.*sum(exp(-beta*hpoth(C)));
%     corrBoltz = (exp(-beta*(hpoth(C))).*exp(beta*diff(C)) )./N.*sum(exp(-beta*(hpoth(C))).*exp(beta*diff(C)));
%     
%     figure(1)
%     plot(C,orgBoltz)
%     
%     figure(2)
%     plot(C,pertBoltz)
%     
%     figure(3)
%     plot(C,corrBoltz)
    
%     orgBoltz = exp(-beta*poth(Xtemp))./sum(exp(-beta*poth(Xtemp)));
%     pertBoltz = exp(-beta*hpoth(Xtemp))./sum(exp(-beta*hpoth(Xtemp)));
%     corrBoltz = (exp(-beta*(hpoth(Xtemp))).*exp(beta*diff(Xtemp)) )./sum(exp(-beta*(hpoth(Xtemp))).*exp(beta*diff(Xtemp)));
%     
%     subplot(2,3,1)
%     plot(Xtemp,orgBoltz); 
%     title('Orgiginal Boltzmann Sampled')
%     %hold on
%     %hold off
%     
%     subplot(2,3,2)
%     plot(Xtemp,pertBoltz)
%     title(['Perturbed Boltzmann Sampled t= ',num2str(t)])
%     
%     subplot(2,3,3)
%     plot(Xtemp,corrBoltz);
%     title('Corrected Boltzmann Sampled')
%     % wenn man die richtige Konstante nimmt dann wird auch die richtige
%     % Boltzmannverteilung gesampelt sum(exp(-beta*(hpoth(x))).*exp(beta*diff(x))
%     
%     subplot(2,3,4)
%     plot(x,exp(-beta*poth(x))./sum(exp(-beta*poth(x))))
%     title('Orgiginal Boltzmann')
%     
%     subplot(2,3,5)
%     plot(x,exp(-beta*hpoth(x))./sum(exp(-beta*hpoth(x))))
%     title('Perturbed Boltzmann')
%     
% %     subplot(2,3,6)
% %     plot(x,exp(-beta*(poth(x)))./sum(exp(-beta*(poth(x))))-exp(-beta*(hpoth(x)))./sum(exp(-beta*(hpoth(x)))))
% %     title('Difference Distribution')
%     
%     subplot(2,3,6)
%     plot(Xtemp,orgBoltz-corrBoltz)
%     title('Difference Sampled')
   
    