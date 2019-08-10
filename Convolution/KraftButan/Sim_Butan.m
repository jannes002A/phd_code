% Skipt für die Simulation eines Butans

function [pos,f,ang] = Sim_Butan(pos,steps,flag)

[n,m]=size(pos);
N=10^7; dt=1/N;
epsilon = 5*10^2;
ang=zeros(1,steps);

for i=1:steps
        [f,gradV] = fast_pot_butan(pos);
        %[f,gradV,~] = pot_butan(pos);
        norm(gradV);
        pos(:) = pos(:) - gradV*dt + sqrt(dt)*epsilon*randn(42,1); 
        pos =  reshape(pos,n,m);
        
        if flag==1
        figure(1)
        plot3(pos(1,:),pos(2,:),pos(3,:),'ro'); hold on
        % Verdindung des Backbone
        line(pos(1,[1,5,8,11]),pos(2,[1,5,8,11]),pos(3,[1,5,8,11]),'Color','k')
        % Verbingungen zum 1. C Atom
        line(pos(1,[1,2]),pos(2,[1,2]),pos(3,[1,2]),'Color','b');
        line(pos(1,[1,3]),pos(2,[1,3]),pos(3,[1,3]),'Color','b');
        line(pos(1,[1,4]),pos(2,[1,4]),pos(3,[1,4]),'Color','b');
        % Verbingungen zum 2. C Atom
        line(pos(1,[5,6]),pos(2,[5,6]),pos(3,[5,6]),'Color','b');
        line(pos(1,[5,7]),pos(2,[5,7]),pos(3,[5,7]),'Color','b');
        % Verbingungen zum 3. C Atom
        line(pos(1,[8,9]),pos(2,[8,9]),pos(3,[8,9]),'Color','b');
        line(pos(1,[8,10]),pos(2,[8,10]),pos(3,[8,10]),'Color','b');
        % Verbingungen zum 4. C Atom
        line(pos(1,[11,12]),pos(2,[11,12]),pos(3,[11,12]),'Color','b');
        line(pos(1,[11,13]),pos(2,[11,13]),pos(3,[11,13]),'Color','b');
        line(pos(1,[11,14]),pos(2,[11,14]),pos(3,[11,14]),'Color','b');
        axis equal
        pause(0.1)
        hold off
        end
        
        % Winkel
        cpos=[pos(:,1),pos(:,5),pos(:,8),pos(:,11)];
        ang(i) = torsionAngButan(cpos);
end


