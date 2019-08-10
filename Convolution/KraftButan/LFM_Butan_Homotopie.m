load start_butan2.txt -ascii
%catoms = [start_pos(:,1),start_pos(:,5),start_pos(:,8),start_pos(:,11)];
%torsionAngButan(catoms);

    %rng(12345)
    % Initialisierung
    delta = 5*10^-5;
    L = 20;
    %nSamples=100; %M
    varplot=0;
    %V = @(pos) pot_butan(pos);

    K =@(p,w) (p'*(p./w))/2;
    
    m=ones(3,14);
    m(:,1)=12;m(:,5)=12;m(:,8)=12;m(:,11)=12;
    m=m(:);
    
    x0= start_butan2(:);
    beta=0.02;
    %lambda=[0.15,0.2,0.25];
    lambda=0.2;
    nS=1E6;
    nT=1;
    %times=zeros(length(lambda),nS);
    
    
for l=1:length(lambda)
    disp(l)
    
    angle=zeros(nS,nT);
     
  
    for nt=1:nT
        t=1;
        %angle = 180;
        x = x0; 
        disp(nt)
        for ns=1:nS
        %while (angle >= 120)
            %fprintf('Sampling %3i \n', t )
            t=t+1;
            %random Momentum
            p0 = sqrt(m/beta).*randn(42,1);
            %First Momentum Step
            [~,dV] = fast_KraftButan_hom(x,lambda(l));
            pStar = p0 - delta/2 * dV;
            %First Postion Step
            xStar = x + delta*pStar./m;

            for i=1:L-1
                %Sampling
                [~,dVStar] = fast_KraftButan_hom(xStar,lambda(l));
                pStar = pStar - delta * dVStar;
                xStar = xStar + delta*pStar./m;
            end
            % Last Momentum Step
            [~,dVStar] = fast_KraftButan_hom(xStar,lambda(l));
            pStar = pStar -delta/2 *dVStar;



            %Evalutate Energies
            U0 = fast_KraftButan_hom(x,lambda(l));
            UStar = fast_KraftButan_hom(xStar,lambda(l));

            K0 = K(p0,m);
            KStar=K(pStar,m);

            %Acceptance criterion
            alpha = min(1,exp(beta*((U0 + K0)-(UStar + KStar))));

            u=rand;


            if u < alpha
                x = xStar;
            %else x=x;
            end

            %Winkel
            pos=reshape(x,3,14);
            cpos=[pos(:,1),pos(:,5),pos(:,8),pos(:,11)];
            angle(ns,nt) = torsionAngButan(cpos);
            
%             if t>1E5
%               x=x0;
%               t=1;
%               angle=180;
%               fprintf('RESTART')
%             end

%             b(t) = (u < alpha);
%             if varplot==1
%             %Plots
%                 pos = reshape(x,3,14);
%                 figure(1);clf
%                 plot3(pos(1,:),pos(2,:),pos(3,:),'r*'); hold on
%                 %Verdindung des Backbone
%                 line(pos(1,[1,5,8,11]),pos(2,[1,5,8,11]),pos(3,[1,5,8,11]))
%                 %Verbingungen zum 1. C Atom
%                 line(pos(1,[1,2]),pos(2,[1,2]),pos(3,[1,2]),'Color','k');
%                 line(pos(1,[1,3]),pos(2,[1,3]),pos(3,[1,3]),'Color','k');
%                 line(pos(1,[1,4]),pos(2,[1,4]),pos(3,[1,4]),'Color','k');
%                 %Verbingungen zum 2. C Atom
%                 line(pos(1,[5,6]),pos(2,[5,6]),pos(3,[5,6]),'Color','k');
%                 line(pos(1,[5,7]),pos(2,[5,7]),pos(3,[5,7]),'Color','k');
%                 %Verbingungen zum 3. C Atom
%                 line(pos(1,[8,9]),pos(2,[8,9]),pos(3,[8,9]),'Color','k');
%                 line(pos(1,[8,10]),pos(2,[8,10]),pos(3,[8,10]),'Color','k');
%                 %Verbingungen zum 4. C Atom
%                 line(pos(1,[11,12]),pos(2,[11,12]),pos(3,[11,12]),'Color','k');
%                 line(pos(1,[11,13]),pos(2,[11,13]),pos(3,[11,13]),'Color','k');
%                 line(pos(1,[11,14]),pos(2,[11,14]),pos(3,[11,14]),'Color','k');
%                 axis equal
%                 pause(0.01)
%                 hold off
%             end
        end
        %name=['Angle',num2str(l)];
        name='angle1rng';
        save(name,'angle','-v7.3');
        %times(l,h)=t;
    end
end
LFM_Butan
%save('Times.mat','times')