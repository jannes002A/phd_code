%Driver Butan


% Ausführen der Lep Frog Simulation für das Butan
% easy
LFM_Butan

%% 
% Simulation Butan mit overdamped Langevin
% Sim_Butan(Position der Atome, Anzahl der Durchläufe, Zeichnen(0=nein,1=ja))

load start_butan2.txt 

x_opt=start_butan2;

[F,grad,ang]=Sim_Butan(x_opt,100,1);
