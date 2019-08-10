%Skript für ein 2D Potential mit Homotopie
%clear all
%close all


epsilon=0.6;
dt=1/1000;
lambda=[0.02,0.025];
Ns=100;
times=zeros(length(lambda),Ns);

for k=1:length(lambda)
    fprintf('Lambda %3i \n',k)
    
    for i=1:Ns
        fprintf('Sampling %3i \n', i )
        
        Xtempx = -1;
        Xtempy = 0;

        j=1;

        while (Xtempx <= 0 && Xtempy <= 1)

                [~,dxf_l,dyf_l]=HomotopieTest2D(Xtempx,Xtempy,lambda(k));
                Xtempx = Xtempx -(dxf_l)*dt +epsilon*sqrt(dt)*randn(1);
                Xtempy = Xtempy -(dyf_l)*dt +epsilon*sqrt(dt)*randn(1);

                j=j+1;

            if (isnan(Xtempx) || isnan(Xtempy))
                error('X=nan!')
            end

            % restart if j is too big
            if j>1E6
            Xtempx = -1;
            Xtempy = 0;

            j=1;
            fprintf('Restarting \n')
            end
        
        
        
        end
        
        times(k,i)=j*dt;
    end
    
    

end
%%
% Die Annahme, dass eine gößere Glättung zu einer kürzeren Stoppzeit führt
% ist manchmal schwer einzuhalten wenn die Varianz zu hoch ist

mtimes=mean(times,2);
vartimes=var(times,0,2);
rates=1./mtimes;


l=linspace(eps,lambda(end));


p=[lambda',[1 1 ]']\log(rates);

f2=@(x) exp(p(1).*x+p(2));
figure(5)
plot(l,f2(l));hold on
plot(lambda,rates,'rx');hold off
