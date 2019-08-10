function [hit,t,k] = Homotopy_butan(start_pos,par_hom,nSamples)  

%Skipt for Beispiel von Leap Frog

catoms = [start_pos(:,1),start_pos(:,5),start_pos(:,8),start_pos(:,11)];
torsionAngButan(catoms);
%pos=start_pos;
hit=0;

%rng(12345)
% Initialisierung
delta = 2.5*10^-4;
L = 20;
%nSamples=100; %M
% varplot=0;
%V = @(pos) pot_butan(pos);

K =@(p,w) (p'*(p./w))/2;

x=zeros(42,nSamples);
x0= start_pos(:);
x(:,1) = x0;

%Gewichtsmatrix
m=ones(3,14);
m(:,1)=12;m(:,5)=12;m(:,8)=12;m(:,11)=12;
m=m(:);
t=1;
%beta=0.05;
beta=0.1;

while (t<nSamples)
    t=t+1;
    %random Momentum
    p0 = sqrt(m/beta).*randn(42,1);
    %First Momentum Step
    [~,dV] = fast_KraftButan_hom(x(:,t-1),par_hom);
    pStar = p0 - delta/2 * dV;
    %First Postion Step
    xStar = x(:,t-1) + delta*pStar./m;
    
    for i=1:L-1
        %Sampling
        [~,dVStar] = fast_KraftButan_hom(xStar,par_hom);
        pStar = pStar - delta * dVStar;
        xStar = xStar + delta*pStar./m;
    end
    % Last Momentum Step
    [~,dVStar] = fast_KraftButan_hom(xStar,par_hom);
    pStar = pStar -delta/2 *dVStar;
    
    x(:,t) = xStar;
    
    %Evalutate Energies
    U0 = fast_KraftButan_hom(x(:,t-1),par_hom);
    UStar = fast_KraftButan_hom(xStar,par_hom);
    
    K0 = K(p0,m);
    KStar=K(pStar,m);
    
    %Acceptance criterion
    alpha = min(1,exp(beta*((U0 + K0)-(UStar + KStar))));
    
    u=rand;
    
    
    if u < alpha
        x(:,t) = xStar;
    else
       x(:,t) = x(:,t-1);
    end
    
    %Winkel
    pos=reshape(x(:,t),3,14);
    cpos=[pos(:,1),pos(:,5),pos(:,8),pos(:,11)];
    angle = torsionAngButan(cpos);
    
    if abs(angle) <= 120
        hit=1;
        %fprintf('Number of Steps %4i \n',t)
        break
    end
    b(t) = (u < alpha);
%     if varplot==1
%     Plots
%         pos = reshape(x(:,t),3,14);
%         figure(1);clf
%         plot3(pos(1,:),pos(2,:),pos(3,:),'r*'); hold on
%         Verdindung des Backbone
%         line(pos(1,[1,5,8,11]),pos(2,[1,5,8,11]),pos(3,[1,5,8,11]))
%         Verbingungen zum 1. C Atom
%         line(pos(1,[1,2]),pos(2,[1,2]),pos(3,[1,2]),'Color','k');
%         line(pos(1,[1,3]),pos(2,[1,3]),pos(3,[1,3]),'Color','k');
%         line(pos(1,[1,4]),pos(2,[1,4]),pos(3,[1,4]),'Color','k');
%         Verbingungen zum 2. C Atom
%         line(pos(1,[5,6]),pos(2,[5,6]),pos(3,[5,6]),'Color','k');
%         line(pos(1,[5,7]),pos(2,[5,7]),pos(3,[5,7]),'Color','k');
%         Verbingungen zum 3. C Atom
%         line(pos(1,[8,9]),pos(2,[8,9]),pos(3,[8,9]),'Color','k');
%         line(pos(1,[8,10]),pos(2,[8,10]),pos(3,[8,10]),'Color','k');
%         Verbingungen zum 4. C Atom
%         line(pos(1,[11,12]),pos(2,[11,12]),pos(3,[11,12]),'Color','k');
%         line(pos(1,[11,13]),pos(2,[11,13]),pos(3,[11,13]),'Color','k');
%         line(pos(1,[11,14]),pos(2,[11,14]),pos(3,[11,14]),'Color','k');
%         axis equal
%         pause(0.1)
%         hold off
%     end
end
k=sum(b);
% ang=zeros(1,nSamples);
% sum(b(b==1))

% for i=1:nSamples
%     pos=reshape(x(:,i),3,14);
%     cpos=[pos(:,1),pos(:,5),pos(:,8),pos(:,11)];
%     ang(i)=torsionAngButan(cpos);
% end