%LFM für Butan

%Skipt for Beispiel von Leap Frog
load start_butan2.txt -ascii
catoms = [start_butan2(:,1),start_butan2(:,5),start_butan2(:,8),start_butan2(:,11)];
torsionAngButan(catoms)
pos=start_butan2;
hit=0;

rng(12345)
% Initialisierung
delta = 2.5*10^-4;
L = 20;
%nSamples=1000; %M
%varplot=1;
%V = @(pos) pot_butan(pos);
nS=1E6;
nT=1;


K =@(p,w) (p'*(p./w))/2;

%x=zeros(42,nSamples);

%Gewichtsmatrix
m=ones(3,14);
m(:,1)=12;m(:,5)=12;m(:,8)=12;m(:,11)=12;
m=m(:);
t=1;
beta=0.05;
angle=zeros(nS,nT);

%helpr=0.1*ones(42,1);

for nt=1:nT
    x0= start_butan2(:);
    x = x0;

    for ns=1:nS
        t=t+1;
        %random Momentum
        p0 = sqrt(m/beta).* randn(42,1);
        %First Momentum Step
        [~,dV] = pot_butan(x);
        pStar = p0 - delta/2 * dV;
        %First Postion Step
        xStar = x + delta*pStar./m;

        for i=1:L-1
            %Sampling
            [~,dVStar] = pot_butan(xStar);
            pStar = pStar - delta * dVStar;
            xStar = xStar + delta*pStar./m;
        end
        % Last Momentum Step
        [f,dVStar] = pot_butan(xStar);
        pStar = pStar -delta/2 *dVStar;

        %x(:,t) = xStar;

        %Evalutate Energies
        U0 = pot_butan(x);
        UStar = pot_butan(xStar);

        K0 = K(p0,m);
        KStar=K(pStar,m);

        %Acceptance criterion
        alpha = min(1,exp(beta*((U0 + K0)-(UStar + KStar))));

        u=rand;

        %b(t) = (u < alpha);

        if u < alpha
            x = xStar;
        %else
           %x(:,t) = x(:,t-1);
        end

        %Winkel
        pos=reshape(x,3,14);
        cpos=[pos(:,1),pos(:,5),pos(:,8),pos(:,11)];
        angle(ns,nt) = torsionAngButan(cpos);

    %     if abs(angle) <= 120
    %         hit=1;
    %         fprintf('Number of Steps %4i \n',t)
    %         break
    %     end

    %     if varplot==1
    %     %Plots
    %         pos = reshape(x(:,t),3,14);
    %         figure(1);clf
    %         plot3(pos(1,:),pos(2,:),pos(3,:),'r*'); hold on
    %         % Verdindung des Backbone
    %         line(pos(1,[1,5,8,11]),pos(2,[1,5,8,11]),pos(3,[1,5,8,11]))
    %         % Verbingungen zum 1. C Atom
    %         line(pos(1,[1,2]),pos(2,[1,2]),pos(3,[1,2]),'Color','k');
    %         line(pos(1,[1,3]),pos(2,[1,3]),pos(3,[1,3]),'Color','k');
    %         line(pos(1,[1,4]),pos(2,[1,4]),pos(3,[1,4]),'Color','k');
    %         % Verbingungen zum 2. C Atom
    %         line(pos(1,[5,6]),pos(2,[5,6]),pos(3,[5,6]),'Color','k');
    %         line(pos(1,[5,7]),pos(2,[5,7]),pos(3,[5,7]),'Color','k');
    %         % Verbingungen zum 3. C Atom
    %         line(pos(1,[8,9]),pos(2,[8,9]),pos(3,[8,9]),'Color','k');
    %         line(pos(1,[8,10]),pos(2,[8,10]),pos(3,[8,10]),'Color','k');
    %         % Verbingungen zum 4. C Atom
    %         line(pos(1,[11,12]),pos(2,[11,12]),pos(3,[11,12]),'Color','k');
    %         line(pos(1,[11,13]),pos(2,[11,13]),pos(3,[11,13]),'Color','k');
    %         line(pos(1,[11,14]),pos(2,[11,14]),pos(3,[11,14]),'Color','k');
    %         axis equal
    %         pause(0.1)
    %         hold off
    %     end
    end
    % x(:,nSamples)
    % ang=zeros(1,nSamples);
    % %sum(b(b==1))
    % 
    % for i=1:nSamples
    %     pos=reshape(x(:,i),3,14);
    %     cpos=[pos(:,1),pos(:,5),pos(:,8),pos(:,11)];
    %     ang(i)=torsionAngButan(cpos);
    % end
    % 
    % figure(1)
    % hist(ang)
end
 save('angle0rng','angle','-v7.3');