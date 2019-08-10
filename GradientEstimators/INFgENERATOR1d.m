%%%
%%% Aufruf: INFgENERATOR1d
%%%
%%% discretization and spectral analysis for the
%%% high friction Langevin (backwards) equation,
%%% namely the infinitessimal genrator, in the
%%% case of a 1d potential
 
%%% parameters %%%

clear all

 script_runs = 0; %Me: why do this?

if ~script_runs
    
 clear all
 clear functions
 script_runs = 0;
 
%%% choose a method: 'LMP' or 'hfLMP' or 'Lioville'
%%% and a potential: 'sym1wellPot' or 'sym2wellPot' or 'asym2wellPot' or 
%%%                  'asym3wellPot' or 'paper3wellPot'

 METHOD    = 'hfLMP';  % Me: high friction Langevin Markov Process
 POT       = 'asym2wellPot'; 
 NO_BOXES  = 1500;            % choose NO_BOXES-2 in the interval [lBoundary,rBoundary]    

 %Me: this routine just switches between the different potentials being
 %used
 switch POT
    case 'sym1wellPot'; 
       	GAMMA     = 2.0;      % friction constant for the process 
       	SIGMA     = 1.0;      % standard deviation for the Brownian motion
	    BETA      = 2*GAMMA/SIGMA^2;  % inverse temperature
       	R_BOUND   = 3.0;      % and the infinte boxes (-infty,lBoundary)and (R_BOUND,infty)
       	L_BOUND   = -R_BOUND;            
        V         = 's1';
    case 'sym2wellPot'; 
       	GAMMA     = 2.0;      
       	SIGMA     = 1.0; % 1.5
	    BETA      = 2*GAMMA/SIGMA^2;  % inverse temperature
       	R_BOUND   = 2.0; % 2.6
       	L_BOUND   = -R_BOUND; 
        V         = 's2';
    case 'asym2wellPot';      
       	GAMMA     = 1.0;      
       	SIGMA     = 1.0/sqrt(2);      
	    BETA      = 2*GAMMA/SIGMA^2;  % inverse temperature
       	R_BOUND   = 3.0;      
       	L_BOUND   = -R_BOUND;            
        V         = 'a2';
    case 'Asym2wellPot';      
        SIGMA     = 1.0;  % default 1.5;      
        GAMMA     = 1.3; % default 2.25;       
        BETA      = 2*GAMMA/SIGMA^2; 
       	R_BOUND   = 2.0;      
       	L_BOUND   = -R_BOUND;            
        V         = 'A2';
    case 'asym3wellPot'; 
        epsi      = 0.5;
        BETA      = 1/epsi; % default 2.0 
        SIGMA     = sqrt(2*epsi);  % default 1.0;      
        GAMMA     = 1; 
       	R_BOUND   = 5.0;
       	L_BOUND   = -R_BOUND-0.5;            
        V         = 'a3';
    case 'asym5wellPot'; 
        epsi      = 0.4;
        BETA      = 1/epsi; % default 2.0 
        SIGMA     = sqrt(2*epsi);  % default 1.0;      
        GAMMA     = 1; 
       	R_BOUND   = 5.0;
       	L_BOUND   = -R_BOUND-0.5;            
        V         = 'a5';
    case 'rugged5wellPot'; 
        epsi      = 0.1;
        BETA      = 1/epsi; % default 2.0 
        SIGMA     = sqrt(2*epsi);  % default 1.0;      
        GAMMA     = 1; 
       	R_BOUND   = 5.0;
       	L_BOUND   = -R_BOUND-0.5;            
        V         = 'ra5';
    case 'paper3wellPot'; 
       	GAMMA     = 1.0;        
       	SIGMA     = 0.75;
	    BETA      = 2*GAMMA/SIGMA^2;  % inverse temperature
        R_BOUND   = 1.8;        
       	L_BOUND   = -R_BOUND;            
	    V         = 'p3';   
    case 'diffusive'; 
       	GAMMA     = 1.0;        
       	SIGMA     = 0.8;
	    BETA      = 2*GAMMA/SIGMA^2;  % inverse temperature
        R_BOUND   = 8.62;        
       	L_BOUND   = -3;            
        V         = 'diff';
    otherwise;
	error('\n\n Unkown potential in HighFrictionLangevin.m ! \n\n'); 
 end;
 
 %Me: why add this here? Seems unnecessary
 BETA      = 2*GAMMA/SIGMA^2;  % inverse temperature

end; % of if ~script_runs 
 
%%% information about the chosen parameters %%%

%Me: the following stuff is just for printing information to the screen
 if strcmp(METHOD,'LMP') 
     error('\n\n METHOD not yet implemented! \n\n');
     methodName = 'Langevin Markov process';
 elseif strcmp(METHOD,'hfLMP') 
     methodName = 'HIGH FRICTION Langevin Markov process';
 else error('\n\n Unknown METHOD! \n\n');
 end;    

 fprintf('   \n');
 fprintf('   METHOD      =  %s    \n',methodName);
 fprintf('   POTENTIAL   =  %s    \n',POT);
 fprintf('   GAMMA       =  %4.2f \n',GAMMA);
 fprintf('   SIGMA       =  %4.2f \n',SIGMA);
 fprintf('   BETA        =  %2.2f \n',BETA);
 fprintf('   interval    =  [%2.2f,%2.2f] \n',L_BOUND,R_BOUND);
 fprintf('   gridpoints  =  %d \n',NO_BOXES);
 fprintf('   \n');
  
%%%------------  discretize infinitisimal generator -----------%%%
%Me: important stuff follows. This is possibly the most important part of
%this entire code.

 dx = (R_BOUND-L_BOUND)/(NO_BOXES-1); %Me: spatial discretization
 x  = (L_BOUND:dx:R_BOUND)'; %Me: the spatial grid (num nodes=NO_BOXES)

 %Me: the following is just for printing messages about how long it takes
 %to compute thte inf. gen. nothing particularly important
 fprintf('   compute inf. generator               = ');
 seconds = cputime;
 e = ones(NO_BOXES,1);
 
 %Me: the following is for printing error messages. 
 if strcmp(METHOD,'LMP') 
      error('\n\n METHOD not jet implemented in InfGenerator.m ! \n\n');
      
      %Construction of discretized inf. gen ( a matrix 'A')
 elseif strcmp(METHOD,'hfLMP') 
    A = spdiags(SIGMA^2/(2*GAMMA^2)*e*[1 -2 1]/(dx^2),[-1 0 1],NO_BOXES,NO_BOXES) + ...
     	  spdiags(-1/GAMMA*dPot(x,POT)*[-1 0 1]/(2*dx),[1 0 -1],NO_BOXES,NO_BOXES)';
      %Impose boundary conditions.
      A(1,1) = -A(1,2); 
      A(NO_BOXES,NO_BOXES) = -A(NO_BOXES,NO_BOXES-1);
 elseif strcmp(METHOD,'Liouville') 
      error('\n\n METHOD not jet implemented in InfGenerator.m ! \n\n');
 else error('\n\n Unknown METHOD in InfGenerator.m ! \n\n');
 end;
 
%_____ force A to have row sum equal to zero! 
%_____ rowSum = full(sum(A'));
% figure(4); s=full(sum(A')); plot(s(1:NO_BOXES)); title([METHOD]);

%Me: computation time printouts
 time = cputime-seconds;
 fprintf('%4.2f \n',time);
 
%%%------- solve eigenvalue problem for LEFT eigenfunctions -----------%%%%

 options.disp  = 0;     % display no during iteration
 options.issym = 1;     % matrix is symmetric 
 no_lEV = 10;           % Me: number of left eigenvectors that we want 
 
 fprintf('   solve for LEFT eigenfunctions        = ');
 tic;
 %Me: below we want the ten left eigenvectors corresponding to the ten
 %smallest eigenvalues 
 [lEV,EW] = eigs(A',no_lEV,'sm',options);  % 'sm' = ev smallest in magnitude
 time = toc; minutes = floor(time/60); seconds = round(time-60*minutes);
 fprintf('%d min %2d sec\n', minutes, seconds);
 
 %Me: extract vector of eigenvalues
 EW = diag(EW);
 %Me: the invariant density is the left eigenvector corresponding to the
 %smallest eigenvalue
 invDensity = lEV(:,1); 
 %Me: Z is the partition function / normalization constant 
 Z = sum(invDensity)*dx;
 %Me: normalize the left eigenvector
 invDensity = 1/Z * invDensity;

 %Me: compute the canonical density. Pot is a function in the directory 
 %    which computes the potential energy, transforms it to the Boltzmann.
 canDensity = exp(-BETA*Pot(x,POT)) / ( sum( exp(-BETA*Pot(x,POT)) )*dx); 
 
 %Me: L1_dev is the L1 norm of the difference between invariant and
 %canonical densities
 L1_dev = sum(abs( invDensity - canDensity )*dx);
 fprintf('   L1 norm error between inv. densities = %6.4e \n',L1_dev);

 % normalizing w.r.t. the invariant density
 
 %Me: L1_norm of what? Guess: weighted L1_norm of the eigenvectors.
 L1_norm = sum(abs(lEV).*(invDensity*ones(1,no_lEV))*dx);
 
 %Me: modify the left eigenvectors so they are norm 1 and all same sign
 lEV = lEV * diag(1./L1_norm) * diag(sign(sum(lEV))); % last term for sign
 
 %Me: find the L1 error between the eigenvector lEV'*A and diag(EW)*lEV'
 L1_err = sum(abs(lEV'*A-diag(EW)*lEV')'*dx);
 
%%%%%%%%%%%%%%%%% graphics of potential and can. density  %%%%%%%%%%%%%%%%%%%%%%%%%%

 dz = 0.01;
 z  = (L_BOUND:dz:R_BOUND)';

 figure(1); clf % 'clf' clears the current figure window
 
 subplot(2,1,1); 
 plot(z,Pot(z,POT)); 
 hold on; 
 plot(x,0,'r:');
 if strcmp(METHOD,'hfLMP') % the Schroedinger potential
  U_pot = 1/(2*SIGMA^2)*dPot(z,POT).^2 - 1/(2*GAMMA)*ddPot(z,POT);
  plot(z,U_pot-min(U_pot),'r--'); 
 end;
 hold off; 
 title(methodName);
 ylabel('potential');
 axis([L_BOUND R_BOUND 0 min([max([Pot(L_BOUND,POT) Pot(R_BOUND,POT)]) 8])]);
 
 subplot(2,1,2);
 canDensity = exp(-BETA*Pot(z,POT)); 
 canDensity = 1/(sum(canDensity)*dz)*canDensity;
 plot(z,canDensity);
 hold on;
 plot(x,invDensity,'r--');
 if min(invDensity)<0
     plot(x(find(invDensity<0)),0,'r*');
 end;
 hold off;
 titletext = sprintf('BETA = %2.2f, gamma = %4.2f,  #boxes = %d.  dx = %2.4f',...
        BETA,GAMMA,NO_BOXES,dx);
 title(titletext);
 ylabel('can. density');
 axis([L_BOUND R_BOUND 0 1.2*max(canDensity)]);
 
%%%%%%%%%%%%%%%%% graphics of LEFT eigenfunctions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
 maxValue = 1.1*max(max(lEV)); minValue = 1.1*min(min(lEV));
 l=4;
 figure(2); clf;
 for k=[1:l]
    subplot(l,1,k);
    if k>1 
	plot(x,0,'r:',x,lEV(:,k),'b'); 
    %Me: plot a red line at 0 and a blue line for the kth eigenfunction
    else
	warning off
	plot(x,lEV(:,1),'b'); 
    %Me: plot the left eigenfunction corresponding to the first eigenvalue.
    %But why do we need the warning? Why need if k>1?  
	warning on 
    end;
    titletext = sprintf('left eigenfunction of A for \\lambda_%d = %5.4f',k,EW(k));
    title(titletext);
    axis([L_BOUND R_BOUND minValue(1) maxValue(1)]);
 end;

% testing against exact eigenvectors in case of the parabolic potential
 if strcmp(POT,'sym1wellPot')
    Pot_s1_exact_EF;
 end;
  %Me: I do not understand what this does. Why do we need this?

%%%%%% direct solver !!!!

 directSolver = 0;
 %Me: why do we even need the following 'if' loop? It doesn't even run
 %since we set directSolver=0. If you set directSolver=1, then the routine
 %below in the if loop prints out figures which correspond to the figures
 %produced by the preceding loop (flipped about the x-axis)
 if directSolver

   [LEV,LEW] = eig(full(A')); 
   %Me: 'full' converts sparse A' to full A'. Just a storage thing. Why?
   
   LEW = diag(LEW); [ew,index] = sort(-LEW); 
   LEW = LEW(index(1:no_lEV));  LEV = LEV(:,index(1:no_lEV)); 
   %Me:sort eigenvalues, choose lEV smallest eigenvalues and associated
   %eigenfunctions
   
   L1_norm = sum(abs(LEV).*(invDensity*ones(1,no_lEV))*dx); 
   %Me: L1_norm is a row of the invDensity-norms of the eigenfunctions
   LEV = LEV * diag(1./L1_norm) * diag(sign(sum(lEV))); % last term for sign
   %Me: normalize the eigenfunctions so they have L1-norm 1
   %Me: if sum(lEV(:,i))=1, then the sum of the ith eigenvector is
   %positive; if 0, then zero; if -1, then negative. So the command above
   %modifies the entries of the eigenfunctions so that their sum is always 
   %positive and equal to 1.
   
   maxValue = 1.1*max(max(LEV)); minValue = 1.1*min(min(LEV)); 
   %Me: min(min(LEV)) is negative
   l=4;
   figure(5); clf;
   for k=[1:l]
 	subplot(l,1,k);
 	if k>1 
 	      plot(x,0,'r:',x,LEV(:,k),'b'); 
 	else
 	      warning off
 	      plot(x,LEV(:,k),'b');
 	      warning on 
 	end;
 	titletext = sprintf('Left eigenfunction of A for \\lambda_%d = %5.4f',k,LEW(k));
 	title(titletext);
 	axis([L_BOUND R_BOUND minValue maxValue]);
   end;

 end;
 
%%%------- solve eigenvalue problem for RIGHT eigenfunctions -----------%%%%

if ~script_runs

 options.disp  = 0;     % display only 5 ev during iteration
 options.issym = 1;     % matrix is symmetric 
 no_rEV = 10;
 
 fprintf('   solve for RIGHT eigenfunctions       = ');
 tic;
 [rEV,ew] = eigs(A,no_rEV,'sm',options);  % 'sm' = ev smallest in magnitude

 time = toc; minutes = floor(time/60); seconds = round(time-60*minutes);
 fprintf('%d min %2d sec\n', minutes, seconds);
 
 EW = diag(ew); 
 Max_norm = max(abs(rEV));
 rEV = rEV * diag(1./Max_norm) * diag(sign(sum(rEV))); % last term for sign
 
 index = [1:no_rEV];
 max_norm_err = max(abs(A*rEV-rEV*diag(EW)));

 fprintf('\n     inf. generator      L1 norm error(lEV)     max norm error(rEV) ');
 fprintf('\n   -----------------    --------------------   --------------------');
 fprintf('\n     %+10.5f             %10.2e             %10.2e',[EW';L1_err;max_norm_err] );
 fprintf('\n          ...    ');
 fprintf('\n'   );
 fprintf('\n   max row sum of A    = %+4.3e     ',full(max(abs(sum(A')))));
 fprintf('\n   min diag entry of A = %+4.3e \n\n',full(min(diag(A))));

else

 fprintf('\n     inf. generator      L1 norm error(lEV) ');
 fprintf('\n   -----------------    --------------------');
 fprintf('\n     %+10.4f             %10.2e',[EW';L1_err] );
 fprintf('\n          ...    ');
 fprintf('\n'   );
 fprintf('\n   max row sum of A    = %+4.3e     ',full(max(abs(sum(A')))));
 fprintf('\n   min diag entry of A = %+4.3e \n\n',full(min(diag(A))));

end; % if ~script_runs    

% addpath('/home/ritz/schuette/huisinga/eMD/matlab/markov/');
%  rEV = SkalierungDerEV(rEV,'PosNegMaxNorm',invDensity,0);
%  kEV = kEV*diag(sign(sum(kEV)));   % to put the sign right
 
% [IndVek,C] = FastInvMengen8(P,rEV,k,invDichte);
% plot1dVek(IndVek(:,1:k),diag(C),2);     % k gefundene Strukturen

%%%%%%%%%%%%%%%%% graphics of RIGHT eigenfunctions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~script_runs

 figure(3); clf;
 r = 4;
 for k=[1:r]
     
     subplot(r,1,k);
     titletext = sprintf('right eigenfunction of A for \\lambda_%d = %5.4f',k,EW(k));
     if k>1 
 	 plot(x,0,'r:',x,rEV(:,k),'b'); 
     else
 	 warning off;
 	 plot(x,rEV(:,k),'b');
 	 warning on;
     end;
     title(titletext);
     axis([L_BOUND R_BOUND -1.1 1.1]);
     
 end;

end; % if ~script_runs 
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% saving results  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


PRINTING = 0;
if PRINTING

    No = 1;
    pathname = '/home/ritz/biocomp/huisinga/eMD/matlab/toolbox/figures';
    filename = sprintf('%s/POT%s_%s_A_b%2.2f_s%3.2f',pathname,V,METHOD,BETA,SIGMA);

% save varios figures

    aS = 3;   %Achsen-st�rke
    lS = 3;   %Linien-st�rke
    fS = 18;  %fontsize   
    figNo = 10;     

    fid = fopen(sprintf('%s_data',filename),'w');
    
    fprintf(fid,'   \n');
    fprintf(fid,'   METHOD        =  %s \n',methodName);
    fprintf(fid,'   POTENTIAL     =  %s    \n',POT);
    fprintf(fid,'   GAMMA         =  %4.2f \n',GAMMA);
    fprintf(fid,'   SIGMA         =  %4.2f \n',SIGMA);
    fprintf(fid,'   BETA          =  %2.2f \n',BETA);
    fprintf(fid,'   interval      =  [%2.2f,%2.2f] \n',L_BOUND,R_BOUND);
    fprintf(fid,'   gridpoints    =  %d \n',NO_BOXES);
    fprintf(fid,'   \n');
    fprintf(fid,'   L1 dev. between inv. densities = %4.2e \n',L1_dev);
    fprintf(fid,'\n     inf. generator      L1 norm error(lEV)     max norm error(rEV) ');
    fprintf(fid,'\n   -----------------    --------------------   --------------------');
    fprintf(fid,'\n     %+10.4f             %10.2e             %10.2e',[EW';L1_err;max_norm_err] );
    fprintf(fid,'\n          ...    ');
    fprintf(fid,'\n'   );
    fprintf(fid,'\n   max row sum of A    = %+4.3e     ',full(max(abs(sum(A')))));
    fprintf(fid,'\n   min diag entry of A = %+4.3e \n\n',full(min(diag(A))));
    
    status = fclose(fid);
    if status==-1
	error('\n\n Error while writing %s.data in InfGenerator! \n\n',filename);
    end;
    
% Potential  
    figure(figNo); clf
    plot(z,Pot(z,POT),'LineWidth',lS);
    axis([L_BOUND R_BOUND 0 min([max([Pot(L_BOUND,POT) Pot(R_BOUND,POT)]) 5])]);
    xlabel('x'); ylabel('V(x)');
    set(gca,'FontSize',fS); set(gca,'LineWidth',aS); set(get(gca,'title'),'Fontsize',fS);  
    set(get(gca,'xlabel'),'Fontsize',fS); set(get(gca,'ylabel'),'Fontsize',fS); 

    eval(sprintf('print -depsc -f%d %s_Pot_%d.eps',figNo,filename,No));
    figNo = figNo+1;
  
% inv&canDensity   
    figure(figNo); clf
    plot(z,canDensity,'LineWidth',lS);
    hold on;
    plot(x,invDensity,'r--','LineWidth',lS);
    hold off;
    axis([L_BOUND R_BOUND 0 1.2*max(invDensity)]);
    xlabel('x'); ylabel('inv. density');
    set(gca,'FontSize',fS); set(gca,'LineWidth',aS); set(get(gca,'title'),'Fontsize',fS);  
    set(get(gca,'xlabel'),'Fontsize',fS); set(get(gca,'ylabel'),'Fontsize',fS); 

    eval(sprintf('print -depsc -f%d %s_InvDen_%d.eps',figNo,filename,No));
    figNo = figNo+1;

% k LEFT eigenfunction
    for k=[1:4]
       figure(figNo); clf
       plot(x,lEV(:,k),'LineWidth',lS);
       hold on; plot(x,zeros(1,length(x)),'r'); hold off;
       titletext = sprintf('\\lambda_%d = %5.4f',k,EW(k));
       axis([L_BOUND R_BOUND minValue maxValue]);
       xlabel('x'); ylabeltext = sprintf('u_%d(x)',k); ylabel(ylabeltext);
       set(gca,'FontSize',fS); set(gca,'LineWidth',aS);   
       set(get(gca,'xlabel'),'Fontsize',fS); set(get(gca,'ylabel'),'Fontsize',fS); 
     
       eval(sprintf('print -depsc -f%d %s_lEF%d_%d.eps',figNo,filename,k,No));
       figNo = figNo+1;
    end;
 
% k RIGHT eigenfunction
    for k=[1:4]
       figure(figNo); clf
       plot(x,rEV(:,k),'LineWidth',lS);
       hold on; plot(x,zeros(1,length(x)),'r'); hold off;
       titletext = sprintf('\\lambda_%d = %5.4f',k,EW(k));
       axis([L_BOUND R_BOUND -1.1 1.1]);
       xlabel('x'); ylabeltext = sprintf('h_%d(x)',k); ylabel(ylabeltext);
       set(gca,'FontSize',fS); set(gca,'LineWidth',aS);   
       set(get(gca,'xlabel'),'Fontsize',fS); set(get(gca,'ylabel'),'Fontsize',fS); 
     
       eval(sprintf('print -depsc -f%d %s_rEF%d_%d.eps',figNo,filename,k,No));
       figNo = figNo+1;
    end;
end;

%Me: the command below saves all the data to 'Lr5Well-0.1.mat'
%save asym2wellPot

% preparing for calling the ModInfGenerator1d  routines
 
%  h = rEV(:,2); Lambda = EW(2);