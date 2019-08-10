%Skipt für Girsanov aus dem Oksendal für Stopzeiten

N=1000;
dt=1/1000000;
R=100;



Ex=zeros(1,R);
Ey=zeros(1,R);
Eyg=zeros(1,R);
Eygg=zeros(1,R);
G=zeros(1,R);
ht=zeros(1,R);


a=1;
beta=3;
set=2.3;


for rep=1:R
    
    Gd=0;
    Gs=0;
    X=2; 
    Y=X;
    i=1;
    first=0;

    for i=1:N-1
    %while( (Y >set ) == 0)

        dBt=sqrt(dt)*randn;
        X = X + a*X*dt + sqrt(2/beta)*dBt;
        Y = Y + sqrt(2/beta)*dBt;

        Gs = Gs - sqrt(beta/2)*(a*Y)*dBt;
        Gd = Gd + beta/2*(a*Y)^2*dt;

         if (X > set && first == 0)
             ht(rep)=i;
             gs=Gs;
             gd=Gd;
             first=1; 
         end
        
%         i=i+1;

    end
    
    Ex(rep)=exp(-beta*ht(rep)*dt);
    Ey(rep)=exp(-beta*i*dt);
    Eyg(rep)=exp(-beta*ht(rep)*dt)*exp(gs+gd);
    Eygg(rep)=(exp(-beta*ht(rep)*dt)*exp(gs+gd))/exp(gs+gd);
    G(rep)=exp(gs+gd);   
    
    
    
%     % die Funktion f ist die Identität
%     Ex(rep)=X(N);
%     Ey(rep)=Y(N);
%     Eyg(rep)=Y(N)*exp(Gs(N)+Gd(N));
%     Eygg(rep)=(Y(N)*exp(Gs(N)+Gd(N)))/exp(Gs(N)+Gd(N));
%     G(rep)=exp(Gs(N)+Gd(N));
    
end

fprintf('Erwartungswert X %2.4f \n',mean(Ex))
fprintf('Erwartungswert Y %2.4f \n',mean(Ey))
fprintf('Erwartungswert YG %2.4f \n',mean(Eyg))
fprintf('Erwartungswert YG/G %2.4f \n',mean(Eyg)/mean(G))
fprintf('Erwartungswert G %2.4f \n',mean(G))