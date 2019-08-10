%Skript für ein 2D Potential mit Homotopie
%clear all
%close all


epsilon=0.6;
dt=1/1000;
lambda=0;
Ns=100;
times=zeros(Ns);



    
    for i=1:Ns
        fprintf('Sampling %3i \n', i )
        
        Xtempx = -1;
        Xtempy = 0;

        j=1;

        while (Xtempx <= 0 && Xtempy <= 1)

                [~,dxf_l,dyf_l]=HomotopieTest2D(Xtempx,Xtempy,0);
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
        
        times(i)=j*dt;
    end
    
    



mtimes=mean(times);
vartimes=var(times);
rates=1./mtimes;

intminus = mtimes-1.96*sqrt(vartimes/Ns);
intplus = mtimes+1.96*sqrt(vartimes/Ns);