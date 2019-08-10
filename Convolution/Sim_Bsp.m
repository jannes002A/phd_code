

pos=[0,0; 0,0.5; 0,-.5; .5,0; -.5,0];
steps=10;
[n,m]=size(pos);
N=10^7; dt=1/N;
epsilon = 5*10^1;

for i=1:steps
        [f,gradV,~] = fast_pot_bsp(pos);
        norm(gradV);
        pos(:) = pos(:) - gradV*dt + sqrt(dt)*epsilon*randn(12,1); 
        pos =  reshape(pos,n,m);
        figure(1)
        plot(pos(:,1),pos(:,2),'o','MarkerSize',10,'MarkerEdgeColor','g','MarkerFaceColor','g'); 
        axis equal
        pause(0.1)
        hold off
end
