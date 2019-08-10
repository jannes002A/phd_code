%Skipt für Girsanov aus dem Oksendal

N=15000;
dt=1/1000;
R=100;



Ex=zeros(1,R);
Ey=zeros(1,R);
Eyg=zeros(1,R);
Eygg=zeros(1,R);
G=zeros(1,R);
Z=zeros(1,R);
hittime=zeros(1,R);


a=.1;
beta=3;
set=2.3;

%Wenn sqrt(beta/2)*(a*Y(i+1)) > 1 ist denn muss die Zeitdiskretisierung klein gewählt werden 

for rep=1:R
    
    Gs = zeros(1,N);
    Gd = zeros(1,N);
    Zt = zeros(1,N);
    X= zeros(1,N);
    Y= zeros(1,N);
    X(1)=2; Y(1)=X(1);
    i=1;
    

    for i=1:N-1
   

        dBt=sqrt(dt)*randn;
        X(i+1) = X(i) + a*X(i)*dt + sqrt(2/beta)*dBt;
        Y(i+1) = Y(i) + sqrt(2/beta)*dBt;

        Gs(i+1) = Gs(i) + sqrt(beta/2)*(a*Y(i+1))*dBt;
        Gd(i+1) = Gd(i) - beta/2*(a*Y(i+1))^2*dt;
        Zt(i+1) = Zt(i) + sqrt(beta/2)*(a*Y(i+1))*dBt;


    end
    
    Ex(rep)=X(N);
    Ey(rep)=Y(N);
    Eyg(rep)=Y(N)*exp(Gs(N)+Gd(N));
    Eygg(rep)=(Y(N)*exp(Gs(N)+Gd(N)))/exp(Gs(N)+Gd(N));
    G(rep)=exp(Gs(N)+Gd(N));   
    Z(rep)=1+Zt(N);
    
    
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
fprintf('Erwartungswert Z %2.4f \n',mean(Z))