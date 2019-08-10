

% Skipt für die Siumlation von Butan
% Distanzmatrix
%   C1 H1 H2 H3 C2 H4 H5 C3 H6 H7 C4 H8 H9 H10

D=[ 0 1.06 1.06 1.06 1.54 4.59 4.59 6.33 0 0 0 0 0 0;
    1.06 0 2.99 2.99 6.33 0 0 0 0 0 0 0 0 0;
    1.06 2.99 0 2.99 6.33 0 0 0 0 0 0 0 0 0;
    1.06 2.99 2.99 0 6.33 0 0 0 0 0 0 0 0 0;
    1.54 6.33 6.33 6.33 0 1.06 1.06 1.54 4.59 4.59 6.33 0 0 0;
    4.59 0 0 0 1.54 00 2.99 4.59 0 0 0 0 0 0;
    2.99 0 0 0 1.06 2.99 0 4.59 0 0 0 0 0 0;  
    6.33 0 0 0 1.54 4.59 4.59 0 1.06 1.06 1.54 4.59 4.59 4.59;
    0 0 0 0 4.59 0 0 1.06 0 2.99 2.99 0 0 0;
    0 0 0 0 4.59 0 0 1.06 2.99 0 4.59 0 0 0;
    0 0 0 0 6.33 0 0 1.54 4.59 4.59 0 1.06 1.06 1.06;
    0 0 0 0 0 0 0 4.59 0 0 1.06 0 2.99 2.99;
    0 0 0 0 0 0 0 4.59 0 0 1.06 2.99 00 2.99;
    0 0 0 0 0 0 0 4.59 0 0 1.06 2.99 2.99 00];  

D=D.^2;

atom = distance_geometry(D,3);
atom=atom';
figure(1)
plot3(atom(1,:),atom(2,:),atom(3,:),'r*'); hold on
line(atom(1,[1,5,8,11]),atom(2,[1,5,8,11]),atom(3,[1,5,8,11])); hold off

%%

atom=randn(3,14);

Sim_Butan(atom,1)

%%

[d60, d180]=torsion_distance(1.06, 1.54, 1.54); 

%%
for i=1:1000
    [x_opt,f]= Sim_Butan(x_opt,1000);
end
%%

load start_butan2.txt 

x_opt=start_butan2;

[F,grad,ang]=Sim_Butan(x_opt,100,1);
%[F_hom,grad_hom]=Sim_Butan(x_opt,100);

%[F,grad]=pot_butan(x_opt);

% WW = [00 12 13 21 22 23 31 32 33 41 42 43];
% for i=1:14
% [f,fx,fy,fz] = KraftButan_hom(5,5,x_opt(:,1),x_opt(:,2),WW(i),1)
% end
% tic
% [F_hom,grad_hom]=pot_butan(x_opt);
% toc
%[F,grad]=fast_pot_butan(x_opt);

%tic
%[F_hom,grad_hom]=fast_KraftButan_hom(x_opt,0);
%toc

%%

x=[0,0,0];
y=[1,1,1];

K=5;
T=5;
KB=0;



  
  for i=1:12
      [f,fx,fy,fz]= KraftButan(K,T,x,y,WW(1,i));
      KB=KB+f;
  end
  
 KB;
 
 %%
 load start_butan2.txt 

x_opt=start_butan2;

[f,fg]=pot_butan(x_opt);
%%
tic
LFM_Butan
toc