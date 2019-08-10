clear 
clc
close all

%% Parameter sets that worked (no nonpositive values of MGF)
%% domain parameters fixed at L_BOUND=[-3.0 -1.5]; R_BOUND=[2.0 1.5]; 
%% discretisation parameters fixed at dx1=0.05; dx2=0.05
%% landscape parameters fixed at:
%%   c=2.0, a=0.8, w_peak_support=0.6; mult_const_inner=2.0; outer=2.0*innr
%%   htzw=0.1*pi, hwog=0.05; hwig=hwog; multconst1=0.65*epsi; g=0.2*epsi;
%%   multconst2=0.03*epsi; w_add_const=0.1

%% KEY PARAMETERS: (largest admissible value of sigma tried, given epsi)
%% sigma=0.0048, epsi=0.0005
%% sigma=0.009, epsi=0.001
%% sigma=0.015, epsi=0.002
%% sigma=0.037, epsi=0.004
% Warning: some pathological MFPT surfaces may result (e.g. sigma=0.0048)

dx1 = 0.05; dx2 = 0.05;
 
% boundaries of the square [L_BOUND(1) R_BOUND(1)]x[L_BOUND(2) R_BOUND(2)]
R_BOUND   = [2.0 1.5];        
L_BOUND   = [-3.0 -1.5];            

NO_GRIDPOINTS1  = ceil((R_BOUND(1)-L_BOUND(1))/dx1)+1;  
NO_GRIDPOINTS2  = ceil((R_BOUND(2)-L_BOUND(2))/dx2)+1;  
NO_GRIDPOINTS = NO_GRIDPOINTS1*NO_GRIDPOINTS2;
dx1 = (R_BOUND(1)-L_BOUND(1))/(NO_GRIDPOINTS1-1); % small correction
dx2 = (R_BOUND(2)-L_BOUND(2))/(NO_GRIDPOINTS2-1); 
dx = dx1*dx2;

%% DEFINE POTENTIAL

epsi=0.001;
c=2.0;
a=0.8;

%smaller values of w_peak_support -> more delta-like
w_peak_support=0.6; 
Ind_Set1=@(x)(1.0*gt(x,a-w_peak_support)).*...
    (1.0*lt(x,a+w_peak_support));
mult_const_inner=2.8;
%% Note: in the original, the 
%     ( cos(blah) + 1.0 )
% part was taken to the first power. We take it to the second power so that
% the potential is twice differentiable. If the potential is only once
% differentiable, then the second derivative is discontinuous, and we get a
% Gibbs phenomenon in solving the pde by finite differences later. This
% leads to negative values of the moment generating function, which is
% problematic.
wfun_inner=@(x)mult_const_inner*Ind_Set1(x).*...
    (cos((pi/w_peak_support)*(x-a))+1.0).^2;
% The price of additional regularity is more pain when it comes to defining
% the function in Matlab. 
gradwfun_inner=@(x)2.0*mult_const_inner*Ind_Set1(x).*...
    (-pi/w_peak_support).*sin((pi/w_peak_support)*(x-a)).*...
    (cos((pi/w_peak_support)*(x-a))+1.0);
lapwfun_inner=@(x)2.0*mult_const_inner*Ind_Set1(x).*...
    (-pi/w_peak_support)^2.*(...
    (sin((pi/w_peak_support)*(x-a))).^2-...
    (cos((pi/w_peak_support)*(x-a))).^2-...
    cos((pi/w_peak_support)*(x-a)));

Ind_Set2=@(x)(1.0*gt(x,-a-w_peak_support)).*...
    (1.0*lt(x,-a+w_peak_support));

mult_const_outer=1.0*mult_const_inner;

wfun_outer=@(x)mult_const_outer*Ind_Set2(x).*...
    (cos((pi/w_peak_support*(x+a)))+1.0).^(2);
gradwfun_outer=@(x)2.0*mult_const_outer*Ind_Set2(x).*...
    (-pi/w_peak_support).*sin((pi/w_peak_support)*(x+a)).*...
    (cos((pi/w_peak_support)*(x+a))+1.0);
lapwfun_outer=@(x)2.0*mult_const_outer*Ind_Set2(x).*...
    (-pi/w_peak_support)^2.*(...
    (sin((pi/w_peak_support)*(x+a))).^2-...
    (cos((pi/w_peak_support)*(x+a))).^2-...
    cos((pi/w_peak_support)*(x+a)));

x=L_BOUND(1):dx1:R_BOUND(1);

figure;
subplot(2,3,1)
plot(x,wfun_inner(x))
xlabel('x')
ylabel('wfun_i(x)')
title('Inner membrane function')
xlim([x(1) x(length(x))])

fd=(wfun_inner(x(3:length(x)))-wfun_inner(x(1:(length(x)-2))))/(2.0*dx1);
subplot(2,3,2)
plot(x,gradwfun_inner(x))
hold on
plot(x(2:(length(x)-1)),fd,'r--')
xlabel('x')
ylabel('grad wfun_i(x)')
title('Gradient, i.m.f.')
xlim([x(1) x(length(x))])

fd=(gradwfun_inner(x(3:length(x)))-...
    gradwfun_inner(x(1:(length(x)-2))))/(2.0*dx1);
subplot(2,3,3)
plot(x,lapwfun_inner(x))
hold on
plot(x(2:(length(x)-1)),fd,'r--')
xlabel('x')
ylabel('lap wfun_i(x)')
title('Laplacian, i.m.f')
xlim([x(1) x(length(x))])

subplot(2,3,4)
plot(x,wfun_outer(x))
xlabel('x')
ylabel('wfun_o(x)')
title('Outer membrane function')
xlim([x(1) x(length(x))])

fd=(wfun_outer(x(3:length(x)))-wfun_outer(x(1:(length(x)-2))))/(2.0*dx1);
subplot(2,3,5)
plot(x,gradwfun_outer(x))
hold on
plot(x(2:(length(x)-1)),fd,'r--')
xlabel('x')
ylabel('grad wfun_o(x)')
title('Gradient, o.m.f.')
xlim([x(1) x(length(x))])

fd=(gradwfun_outer(x(3:length(x)))-...
    gradwfun_outer(x(1:(length(x)-2))))/(2.0*dx1);
subplot(2,3,6)
plot(x,lapwfun_outer(x))
hold on
plot(x(2:(length(x)-1)),fd,'r--')
xlabel('x')
ylabel('lap wfun_o(x)')
title('Laplacian, o.m.f.')
xlim([x(1) x(length(x))])

%% boundaries
htzw=0.05*pi; % half transition zone width
twopi_by_tzw=pi/htzw;

% Note that outer gap is centred at y=0
ogloc=0.7;   % outer gap location (centre)

% outer gap left transition function. Note that we have also performed the
% modification below of taking the important bit to the second power.
og_ltf=@(y)(-(sin(twopi_by_tzw*((y-ogloc)+htzw))./twopi_by_tzw+...
    (y-ogloc))./(2.0*htzw).*...
    lt(y,ogloc).*gt((y-ogloc),-2.0*htzw));
% gradient of outer gap left transition function
grad_og_ltf=@(y)(-(cos(twopi_by_tzw*((y-ogloc)+htzw))+1.0)./(2.0*htzw).*...
    lt(y,ogloc).*gt((y-ogloc),-2.0*htzw));

% outer gap right transition function
og_rtf=@(y)((sin(twopi_by_tzw*((y-ogloc)-htzw))./twopi_by_tzw+...
    (y-ogloc))./(2.0*htzw).*...
    lt((y-ogloc),2.0*htzw).*gt(y,ogloc));
% gradient of outer gap right transition function
grad_og_rtf=@(y)((cos(twopi_by_tzw*((y-ogloc)-htzw))+1.0)./(2.0*htzw).*...
    lt((y-ogloc),2.0*htzw).*gt(y,ogloc));

% outer gap non transition zones
og_ntz=@(y)(1.0*le((y-ogloc),-2.0*htzw)+1.0*ge((y-ogloc),2.0*htzw));
% outer gap function
og=@(y)(og_ltf(y)+og_rtf(y)+og_ntz(y));    
% gradient of outer gap function
grad_og=@(y)(grad_og_ltf(y)+grad_og_rtf(y));

igloc=-0.9;   % inner gap location (centre) 
% inner gap left transition function
ig_ltf=@(y)(-(sin(twopi_by_tzw*((y-igloc)+htzw))./twopi_by_tzw+...
    (y-igloc))./(2.0*htzw).*...
    lt(y,igloc).*gt((y-igloc),-2.0*htzw));
% gradient of inner gap left transition function
grad_ig_ltf=@(y)(-(cos(twopi_by_tzw*((y-igloc)+htzw))+1.0)./(2.0*htzw).*...
    lt(y,igloc).*gt((y-igloc),-2.0*htzw));

% inner gap right transition function
ig_rtf=@(y)((sin(twopi_by_tzw*((y-igloc)-htzw))./twopi_by_tzw+...
    (y-igloc))./(2.0*htzw).*...
    lt((y-igloc),+2.0*htzw).*gt(y,igloc));
% gradient of inner gap right transition function
grad_ig_rtf=@(y)((cos(twopi_by_tzw*((y-igloc)-htzw))+1.0)./(2.0*htzw).*...
    lt((y-igloc),2.0*htzw).*gt(y,igloc));

% inner gap non transition zones
ig_ntz=@(y)(1.0*le((y-igloc),-2.0*htzw)+1.0*ge((y-igloc),2.0*htzw));
% inner gap function
ig=@(y)(ig_ltf(y)+ig_rtf(y)+ig_ntz(y));  
% gradient of inner gap function
grad_ig=@(y)(grad_ig_ltf(y)+grad_ig_rtf(y));

% figure(73);
% subplot(2,2,1)
% plot(x,ig(x))
% xlabel('y')
% ylabel('ig(y)')
% title('Inner gap function')
% xlim([x(1) x(length(x))])
% 
% subplot(2,2,2)
% plot(x,grad_ig(x))
% xlabel('y')
% ylabel('grad ig(y)')
% title('Gradient, i.g.f.')
% xlim([x(1) x(length(x))])
% 
% subplot(2,2,3)
% plot(x,og(x))
% xlabel('y')
% ylabel('og(y)')
% title('Outer gap function')
% xlim([x(1) x(length(x))])
% 
% subplot(2,2,4)
% plot(x,grad_og(x))
% xlabel('y')
% ylabel('grad og(y)')
% title('Gradient, o.g.f.')
% xlim([x(1) x(length(x))])

%% POTENTIAL
[X,Y]=meshgrid(x,x);

% larger values of multconst1 make the energetic barrier between 
% [outside membrane] and [inside membrane] larger
multconst1=0.65*epsi; 
% larger values of g will drive the system to the right
g=0.2*epsi; 
% larger values of multconst2 just make the height of the barriers that
% make up the external and internal boundaries larger. The value of
% multconst2 should be large enough to ensure that virtually all 
% transitions from [outside membrane] to [inside membrane] occur via the 
% channel. However, excessively high values of multconst2 lead to extreme
% Gibbs phenomena of the moment gen. funct. at the external membrane, which
% may lead to negative values. 
multconst2=0.04*epsi;
% quad_y_mult_const scales the quadratic-in-y part of the potential
w_add_const=0.1;

%% Parameter sets where all values of the moment gen. funct. strictly +ve:
% multconst1=0.3*epsi, g=0.2*epsi; multconst2=0.01*epsi; w_add_const=0.1;
% same as above, except multconst2=0.02*epsi up to 0.055*epsi; 
% for multconst2=0.06*epsi, there are 16 points for which MGF<=0.

% multconst1=0.5*epsi, g=0.2*epsi; multconst2=0.055*epsi; w_add_const=0.1
% yielded 34 bad points at which MGF <=0
% same parameter set as above but multconst1=0.4*epsi yielded 0 bad points

% multconst1=0.6*epsi, g=0.2*epsi; multconst2=0.05epsi; w_add_const=0.1
% yielded 34 bad points at which MGF <=0
% same parameter set as above but multconst2=0.045*epsi yielded 38 bad pts
% same parameter set as above but multconst2=0.04*epsi yielded 38 bad pts
% same parameter set as above but multconst2=0.03*epsi yielded 0 bad pts

% multconst1=0.7*epsi, g=0.2*epsi; multconst2=0.05epsi; w_add_const=0.1
% yielded 34 bad points at which MGF <=0

% Define the potential and its gradient
Pot=@(x,y)(multconst1.*(x.^2-c^2).^2-g.*x+...
    multconst2.*(w_add_const.*(y.^2)+...
    (wfun_inner(x).^2).*ig(y)+(wfun_outer(x).^2).*og(y)));
gradPot=@(x,y)[(4.0*multconst1).*(x.^2-c^2).*x-g+...
    (2.0*multconst2)*(wfun_inner(x).*gradwfun_inner(x).*ig(y)+...
    wfun_outer(x).*gradwfun_outer(x).*og(y));...
    multconst2.*(2.0.*w_add_const.*y+(wfun_inner(x).^2).*grad_ig(y)+...
    (wfun_outer(x).^2).*grad_og(y))];

dPot_dx=@(x,y)((4.0*multconst1).*(x.^2-c^2).*x-g.*ones(size(x))+...
    (2.0*multconst2).*(wfun_inner(x).*gradwfun_inner(x).*ig(y)+...
    wfun_outer(x).*gradwfun_outer(x).*og(y)));
dPot_dy=@(x,y)(multconst2*((2.0*w_add_const).*y+...
    ((wfun_inner(x)).^2).*grad_ig(y)+...
    ((wfun_outer(x)).^2).*grad_og(y)));

minV_y_eq_0=min(Pot(x,0.0));
% figure(1);clf
% surf(X,Y,Pot(X,Y)/epsi,'Linestyle','none')
% title('Potential')
% xlabel('x')
% ylabel('y')
% zlabel('V(x,y)/\epsilon')

% figure(2);clf
% surf(X,Y,dPot_dx(X,Y),'Linestyle','none')
% title('x-derivative of Potential')
% xlabel('x')
% ylabel('y')
% zlabel('dV(x,y)/dx')
% 
% figure(3);clf
% surf(X,Y,dPot_dy(X,Y),'Linestyle','none')
% title('y-derivative of Potential')
% xlabel('x')
% ylabel('y')
% zlabel('dV(x,y)/dy')

figure(71);
subplot(2,2,1)
surf(X,Y,Pot(X,Y),'Linestyle','none'),view(-30,45)
title('Potential energy surface')
xlabel('x')
ylabel('y')
zlabel('V(x,y)')
xlim([x(1) x(length(x))])
ylim([x(1) x(length(x))])

subplot(2,2,2)
surf(X,Y,Pot(X,Y),'Linestyle','none'),view(0,90)
title('Potential energy surface')
xlabel('x')
ylabel('y')
zlabel('V(x,y)')
xlim([x(1) x(length(x))])
ylim([x(1) x(length(x))])

subplot(2,2,3)
plot(x,(Pot(x,0.0)-minV_y_eq_0)/epsi,'b')
xlabel('x')
ylabel('\epsilon^{-1}V(x,y)')
title('\epsilon^{-1}V(x,y=0.0) as a function of x')
xlim([x(1) x(length(x))])

xlb=-2.0;
xrb=2.0;
subplot(2,2,4)
plot(x,(Pot(xlb,x)-minV_y_eq_0)/epsi,'k-',...
    x,(Pot(0.0,x)-minV_y_eq_0)/epsi,'k--',...
    x,(Pot(xrb,x)-minV_y_eq_0)/epsi,'k-.',...
    x,(Pot(-a,x)-minV_y_eq_0)/epsi,'r-',...
    x,(Pot(a,x)-minV_y_eq_0)/epsi,'r--',...
    'Linewidth',2.0)
legend(['a=',num2str(xlb)],'a=0.0',['a=',num2str(xrb)],...
    'a=ext','a=int','Location','Best')
xlabel('y')
ylabel('\epsilon^{-1}V(x,y)')
xlim([x(1) x(length(x))])
title('\epsilon^{-1}V(x=a,y) as a function of y')

clear X Y x

%%%%%%%%%%%%%%%%%%%%%%%%%%%  discretize infinitisimal generator  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[x1 x2]  = meshgrid((L_BOUND(1):dx1:R_BOUND(1)),(L_BOUND(2):dx2:R_BOUND(2)));
dx  = dx1*dx2;

% reshape the matrices x1 and x2 into row vectors of length NO_GRIDPOINTS,
% stack them together to form a 2 by NO_GRIDPOINTS matrix, then take the
% transpose to get a NO_GRIDPOINTS by 2 matrix
x   = [reshape(x1,1,NO_GRIDPOINTS);reshape(x2,1,NO_GRIDPOINTS)]';

tic;
fprintf('   compute inf. generator               = ');

%L is a NO_GRIDPOINTS by NO_GRIDPOINTS matrix with 5*NO_GRIDPOINTS nonzero
%entries
L = sparse([],[],0,NO_GRIDPOINTS,NO_GRIDPOINTS,5*NO_GRIDPOINTS); 

% all this scaling stuff is for the Langevin equation, I presume
weightL1 = epsi; 
weightL2 = epsi;

 % will be used to make 2d laplacian
 laplace1 = [1 -2 1]/dx1^2; laplace2 = [1 -2 1]/dx2^2; 
 % will be used to make 2d gradient
 grad1 = [-1 0 1]/(2*dx1); grad2 = [-1 0 1]/(2*dx2); 
 
 for index1=[2:NO_GRIDPOINTS1-1]
    for index2=[2:NO_GRIDPOINTS2-1]
        % 2-vector
        x12 = [x1(index2,index1); x2(index2,index1)]; 
        % evaluate the gradient of potential
        dPot = gradPot(x12(1),x12(2)); 
	
        % row index (1-vector)
        rows = NO_GRIDPOINTS2*(index1-1)+index2; 
        % column indices (3-vector)
        columns = NO_GRIDPOINTS2*([index1-1 index1 index1+1]-1)+index2; 
        % populate the matrix L
        L(rows,columns) = L(rows,columns) + [weightL1*laplace1-dPot(1)*grad1] ;
        
        % the command below does not yield the same result as the command
        % some lines above. 
        columns = NO_GRIDPOINTS2*(index1-1)+[index2-1 index2 index2+1];
        L(rows,columns) = L(rows,columns) + [weightL2*laplace2-dPot(2)*grad2] ;       
    end;
 end;
 % The two L(rows,columns)= blah lines just fill up the infinitesimal
 % generator. The two columns commands matter because we have to encode a
 % partial differential operator on a two-dimensional space into a matrix

 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% start: boundary conditions %%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
 % These appear to be Neumann boundary conditions with right
 % hand side equal to zero, i.e. no flux conditions. They are certainly not
 % Dirichlet boundary conditions.
 
 % bottom left corner of square domain
 for index1=[1]
     laplace1 = [-1 1]/dx1^2; grad1 = [-1 1]/(2*dx1);
     for index2=[1]
	laplace2 = [-1 1]/dx2^2; grad2 = [-1 1]/(2*dx2); 
	x12 = [x1(index2,index1); x2(index2,index1)];
	dPot = gradPot(x12(1),x12(2)); 
	
	rows = NO_GRIDPOINTS2*(index1-1)+index2;
        columns = NO_GRIDPOINTS2*([index1 index1+1]-1)+index2;
	L(rows,columns) = L(rows,columns) + [weightL1*laplace1-dPot(1)*grad1] ;
        
	rows = NO_GRIDPOINTS2*(index1-1)+index2;
	columns = NO_GRIDPOINTS2*(index1-1)+[index2 index2+1];
	L(rows,columns) = L(rows,columns) + [weightL2*laplace2-dPot(2)*grad2] ;
        
    end;
 end;

 % bottom edge minus corners
 for index1=[1]
     laplace1 = [-1 1]/dx1^2; grad1 = [-1 1]/(2*dx1);
    for index2=[2:NO_GRIDPOINTS2-1]
        laplace2 = [1 -2 1]/dx2^2; grad2 = [-1 0 1]/(2*dx2); 
	x12 = [x1(index2,index1); x2(index2,index1)];
	dPot = gradPot(x12(1),x12(2)); 
	
	rows = NO_GRIDPOINTS2*(index1-1)+index2;
        columns = NO_GRIDPOINTS2*([index1 index1+1]-1)+index2;
	L(rows,columns) = L(rows,columns) + [weightL1*laplace1-dPot(1)*grad1] ;
        
	rows = NO_GRIDPOINTS2*(index1-1)+index2;
	columns = NO_GRIDPOINTS2*(index1-1)+[index2-1 index2 index2+1];
	L(rows,columns) = L(rows,columns) + [weightL2*laplace2-dPot(2)*grad2] ;
        
    end;
 end;

 % bottom right corner
 for index1=[1]
     laplace1 = [-1 1]/dx1^2; grad1 = [-1 1]/(2*dx1);
    for index2=[NO_GRIDPOINTS2]
        laplace2 = [1 -1]/dx2^2; grad2 = [-1 1]/(2*dx2); 
	x12 = [x1(index2,index1); x2(index2,index1)];
	dPot = gradPot(x12(1),x12(2)); 
	
	rows = NO_GRIDPOINTS2*(index1-1)+index2;
        columns = NO_GRIDPOINTS2*([index1 index1+1]-1)+index2;
	L(rows,columns) = L(rows,columns) + [weightL1*laplace1-dPot(1)*grad1] ;
        
	rows = NO_GRIDPOINTS2*(index1-1)+index2;
	columns = NO_GRIDPOINTS2*(index1-1)+[index2-1 index2];
	L(rows,columns) = L(rows,columns) + [weightL2*laplace2-dPot(2)*grad2] ;
        
    end;
 end;

 % right edge minus corners
 for index1=[2:NO_GRIDPOINTS1-1]
    laplace1 = [1 -2 1]/dx1^2; grad1 = [-1 0 1]/(2*dx1); 
    for index2=[NO_GRIDPOINTS2]
        laplace2 = [1 -1]/dx2^2; grad2 = [-1 1]/(2*dx2); 
	x12 = [x1(index2,index1); x2(index2,index1)];
	dPot = gradPot(x12(1),x12(2)); 
	
	rows = NO_GRIDPOINTS2*(index1-1)+index2;
        columns = NO_GRIDPOINTS2*([index1-1 index1 index1+1]-1)+index2;
	L(rows,columns) = L(rows,columns) + [weightL1*laplace1-dPot(1)*grad1] ;
        
	rows = NO_GRIDPOINTS2*(index1-1)+index2;
	columns = NO_GRIDPOINTS2*(index1-1)+[index2-1 index2];
	L(rows,columns) = L(rows,columns) + [weightL2*laplace2-dPot(2)*grad2] ;
        
    end;
 end;

 % top right corner
 for index1=[NO_GRIDPOINTS1]
    laplace1 = [1 -1]/dx1^2; grad1 = [-1 1]/(2*dx1); 
    for index2=[NO_GRIDPOINTS2]
        laplace2 = [1 -1]/dx2^2; grad2 = [-1 1]/(2*dx2); 
	x12 = [x1(index2,index1); x2(index2,index1)];
	dPot = gradPot(x12(1),x12(2)); 
	
	rows = NO_GRIDPOINTS2*(index1-1)+index2;
        columns = NO_GRIDPOINTS2*([index1-1 index1]-1)+index2;
	L(rows,columns) = L(rows,columns) + [weightL1*laplace1-dPot(1)*grad1] ;
        
	rows = NO_GRIDPOINTS2*(index1-1)+index2;
	columns = NO_GRIDPOINTS2*(index1-1)+[index2-1 index2];
	L(rows,columns) = L(rows,columns) + [weightL2*laplace2-dPot(2)*grad2] ;
        
    end;
 end;

 % top edge minus corners
 for index1=[NO_GRIDPOINTS1]
    laplace1 = [1 -1]/dx1^2; grad1 = [-1 1]/(2*dx1); 
    for index2=[2:NO_GRIDPOINTS2-1]
        laplace2 = [1 -2 1]/dx2^2; grad2 = [-1 0 1]/(2*dx2); 
	x12 = [x1(index2,index1); x2(index2,index1)];
	dPot = gradPot(x12(1),x12(2)); 
	
	rows = NO_GRIDPOINTS2*(index1-1)+index2;
        columns = NO_GRIDPOINTS2*([index1-1 index1]-1)+index2;
	L(rows,columns) = L(rows,columns) + [weightL1*laplace1-dPot(1)*grad1] ;
        
	rows = NO_GRIDPOINTS2*(index1-1)+index2;
	columns = NO_GRIDPOINTS2*(index1-1)+[index2-1 index2 index2+1];
	L(rows,columns) = L(rows,columns) + [weightL2*laplace2-dPot(2)*grad2] ;
        
    end;
 end;

 % top left corner
 for index1=[NO_GRIDPOINTS1]
    laplace1 = [1 -1]/dx1^2; grad1 = [-1 1]/(2*dx1); 
    for index2=[1]
        laplace2 = [-1 1]/dx2^2; grad2 = [-1 1]/(2*dx2); 
	x12 = [x1(index2,index1); x2(index2,index1)];
	dPot = gradPot(x12(1),x12(2)); 
	
	rows = NO_GRIDPOINTS2*(index1-1)+index2;
        columns = NO_GRIDPOINTS2*([index1-1 index1]-1)+index2;
	L(rows,columns) = L(rows,columns) + [weightL1*laplace1-dPot(1)*grad1] ;
        
	rows = NO_GRIDPOINTS2*(index1-1)+index2;
	columns = NO_GRIDPOINTS2*(index1-1)+[index2 index2+1];
	L(rows,columns) = L(rows,columns) + [weightL2*laplace2-dPot(2)*grad2] ;
        
    end;
 end;

 % left edge minus corners
 for index1=[2:NO_GRIDPOINTS1-1]
    laplace1 = [1 -2 1]/dx1^2; grad1 = [-1 0 1]/(2*dx1); 
    for index2=[1]
        laplace2 = [-1 1]/dx2^2; grad2 = [-1 1]/(2*dx2); 
	x12 = [x1(index2,index1); x2(index2,index1)];
	dPot = gradPot(x12(1),x12(2)); 
	
	rows = NO_GRIDPOINTS2*(index1-1)+index2;
        columns = NO_GRIDPOINTS2*([index1-1 index1 index1+1]-1)+index2;
	L(rows,columns) = L(rows,columns) + [weightL1*laplace1-dPot(1)*grad1] ;
        
	rows = NO_GRIDPOINTS2*(index1-1)+index2;
	columns = NO_GRIDPOINTS2*(index1-1)+[index2 index2+1];
	L(rows,columns) = L(rows,columns) + [weightL2*laplace2-dPot(2)*grad2] ;
        
    end;
 end;

 %%%%%%%%%%%%%%%%%%%%%%%%%% end:  boundary conditions  %%%%%%%%%%%%%%%%%%%%%

 %_____ force A to have row sum equal to zero! 
 %_____ rowSum = full(sum(A'));
 % figure(4); s=full(sum(A')); plot(s(1:NO_GRIDPOINTS)); title([DYNAMIC]);

 time = toc; minutes = floor(time/60); seconds = round(time-60*minutes);
 fprintf('%d min %2d sec\n', minutes, seconds);
 
%  %% Invariant density
%  
%  options.disp  = 0;  % display no during iteration
%  options.issym = 1;  % matrix symmetric?
%  no_lEV = 3;
%  
%  tic;
%  fprintf('   solve for LEFT eigenfunctions        = ');
%  lEV = zeros(NO_GRIDPOINTS,no_lEV); 
%  [lEV,ew] = eigs(L',no_lEV,'sm',options);  % 'sm' = ev smallest in magnitude
%  EW = diag(ew); 
%  
%  time = toc; minutes = floor(time/60); seconds = round(time-60*minutes);
%  fprintf('%d min %2d sec\n', minutes, seconds);
% 
%  invDensity = lEV(:,1);
%  invDensity = invDensity/(sum(invDensity)*dx);
% 
%  canDensity = exp(-Pot(x(:,1),x(:,2))/epsi);
%  Z=sum(canDensity)*dx;
%  canDensity=canDensity/Z;
%  
%  maxValue = max(max(invDensity),max(canDensity));
%  
%  L1_dev = sum(abs( invDensity - canDensity )*dx);
%  fprintf('   L1 norm error between inv. densities = %.2e \n',L1_dev);
%  
% figure;
% subplot(2,2,1)
% surf(x1,x2,reshape(Pot(x(:,1),x(:,2)),NO_GRIDPOINTS2,NO_GRIDPOINTS1),'Linestyle','none');
% title('Potential')
% 
% subplot(2,2,2)
% surf(x1,x2,reshape(canDensity,NO_GRIDPOINTS2,NO_GRIDPOINTS1),'Linestyle','none');                 
% title('Canonical density')

%% First mean passage time
deadset_idcs=find(x(:,1)>=2.0);
deadset=x(deadset_idcs,:);

aliveset_idcs=find(x(:,1)<2.0);
aliveset=x(aliveset_idcs,:);

% deadset_idcs=vertcat(find(x(:,1)==R_BOUND(1) & x(:,2)<R_BOUND(2)),...
%     find(x(:,1)>L_BOUND(1)& x(:,2)==R_BOUND(2)),...
%     find(x(:,1)==L_BOUND(1) & x(:,2)>L_BOUND(2)),...
%     find(x(:,1)<R_BOUND(1)& x(:,2)==L_BOUND(2)));
% deadset=x(deadset_idcs,:);
% 
% aliveset_idcs=find((x(:,1)<R_BOUND(1)&x(:,1)>L_BOUND(1))&...
%     (x(:,2)<R_BOUND(2) & x(:,2)>L_BOUND(2)));
% aliveset=x(aliveset_idcs,:);

%the running cost function is f(x)=1.0 (the lt(z1,3.0) is included in order
%to enforce dependence on the state in some way. It could be replaced with
%any other suitable function of z1 or z2
rnngcostfunk=@(z1,z2)(1.0*lt(z1,3.0)); 
f=rnngcostfunk(aliveset(:,1),aliveset(:,2));

% the terminal cost function is g(x)=0.0. The multiplication with z1 is
% again merely to enforce dependence on the state, and can be replaced with
% any other suitable function of z1 or z2
tmnlcostfunk=@(z1,z2)(0.0*z1);
g=tmnlcostfunk(deadset(:,1),deadset(:,2));

% Create psi(x)=E^x[exp(-sigma W)]
psi1=ones(NO_GRIDPOINTS,1);

%% Warning: If sigma is too large, then the result will not be meaningful
% sigma=0.01 worked
sigma1=0.005;
% Calculate values of psi over killing set
psi1(deadset_idcs)=exp(-sigma1*g);

%% Solve PDE over aliveset, since we do not know the values of psi here
L_aliveset=L(aliveset_idcs,aliveset_idcs);
psi1(aliveset_idcs)=(L_aliveset-sigma1*diag(f))\...
    (-L(aliveset_idcs,deadset_idcs)*psi1(deadset_idcs));

figure(1);clf
surf(x1,x2,reshape(psi1,NO_GRIDPOINTS2,NO_GRIDPOINTS1),'Linestyle','none')
title(['MGF, |{(x_i,y_i) s.t. MGF(x_i,y_i)<=0}|=',num2str(length(find(psi1<=0)))])
xlabel('x')
ylabel('y')
zlabel('MGF(x,y)')

figure(2);clf
surf(x1,x2,reshape(Pot(x1,x2)/epsi,NO_GRIDPOINTS2,NO_GRIDPOINTS1),'Linestyle','none')
title('V(x,y)/\epsilon')
xlabel('x')
ylabel('y')
zlabel('V(x,y)/\epsilon')

if(isequal(length(find(psi1<=0)),0))
    F1=-log(psi1)/sigma1;
    F1_reshaped=reshape(F1,NO_GRIDPOINTS2,NO_GRIDPOINTS1);
    figure(3);clf
    plot(x1(1,:),F1_reshaped(1,:),'k-',x1(1,:),F1_reshaped(15,:),'k--',...
        x1(1,:),F1_reshaped(30,:),'k-.','Linewidth',2)
    xlabel('x')
    xlim([x1(1,1) x1(1,NO_GRIDPOINTS1)])
    legend(['a=',num2str(x2(1,1))],['a=',num2str(x2(15,1))],...
        ['a=',num2str(x2(30,1))])
    ylabel('F(x,y=a)')
    title('F(x,y=a)')

    figure(4);clf
    surf(x1,x2,reshape(F1,NO_GRIDPOINTS2,NO_GRIDPOINTS1),'Linestyle','none')
    title('F(x,y)')
    xlabel('x')
    ylabel('y')
    zlabel('F(x,y)')

    psi2=zeros(size(psi1));
    sigma2=(1.0-sigma1)*sigma1;
    psi2(deadset_idcs)=exp(-sigma2*g);
    psi2(aliveset_idcs)=(L_aliveset-sigma2*diag(f))\...
        (-L(aliveset_idcs,deadset_idcs)*psi2(deadset_idcs));
    F2=-log(psi2)/sigma2;
    MFPT=-(F1-F2)./sigma1^2;
    
    figure(5);clf
    surf(x1,x2,reshape(MFPT,NO_GRIDPOINTS2,NO_GRIDPOINTS1),'Linestyle','none')
    title('MFPT(x,y)')
    xlabel('x')
    ylabel('y')
    zlabel('MFPT(x,y)')
    
    
    dF1_dx=(F1_reshaped(2:(NO_GRIDPOINTS2-1),3:NO_GRIDPOINTS1)-...
        F1_reshaped(2:(NO_GRIDPOINTS2-1),1:(NO_GRIDPOINTS1-2)))/(2.0*dx1);
    dF1_dy=(F1_reshaped(3:NO_GRIDPOINTS2,2:(NO_GRIDPOINTS1-1))-...
        F1_reshaped(1:(NO_GRIDPOINTS2-2),2:(NO_GRIDPOINTS1-1)))/(2.0*dx2);
    
    figure(6);clf
    surf(x1(2:(NO_GRIDPOINTS2-1),(2:NO_GRIDPOINTS1-1)),...
        x2(2:(NO_GRIDPOINTS2-1),(2:NO_GRIDPOINTS1-1)),...
        -2.0*epsi*sigma1*dF1_dx,'Linestyle','none')
    title('Opt. ctrl, x-component')
    xlabel('x')
    ylabel('y')
    zlabel('copt(x,y), x-component')
    
    figure(7);clf
    surf(x1(2:(NO_GRIDPOINTS2-1),(2:NO_GRIDPOINTS1-1)),...
        x2(2:(NO_GRIDPOINTS2-1),(2:NO_GRIDPOINTS1-1)),...
        -2.0*epsi*sigma1*dF1_dy,'Linestyle','none')
    title('Opt. ctrl, y-component')
    xlabel('x')
    ylabel('y')
    zlabel('copt(x,y), y-component')
    
    int_x=reshape(x1(2:(NO_GRIDPOINTS2-1),(2:NO_GRIDPOINTS1-1)),...
        (NO_GRIDPOINTS2-2)*(NO_GRIDPOINTS1-2),1);
    int_y=reshape(x2(2:(NO_GRIDPOINTS2-1),(2:NO_GRIDPOINTS1-1)),...
        (NO_GRIDPOINTS2-2)*(NO_GRIDPOINTS1-2),1);
    dPot_dx_int=dPot_dx(int_x,int_y);
    dPot_dy_int=dPot_dy(int_x,int_y);
    dPot_dx_int_resh=reshape(dPot_dx_int,NO_GRIDPOINTS2-2,NO_GRIDPOINTS1-2);
    dPot_dy_int_resh=reshape(dPot_dy_int,NO_GRIDPOINTS2-2,NO_GRIDPOINTS1-2);
    
    figure(8);clf
    surf(x1(2:(NO_GRIDPOINTS2-1),(2:NO_GRIDPOINTS1-1)),...
        x2(2:(NO_GRIDPOINTS2-1),(2:NO_GRIDPOINTS1-1)),...
        -2.0*epsi*sigma1*dF1_dx-dPot_dx_int_resh,'Linestyle','none')
    title('gradient of opt. tilted potential, x-component')
    xlabel('x')
    ylabel('y')
    zlabel('grad. opt. tilt. pot(x,y), x-component')    
    
    figure(9);clf
    surf(x1(2:(NO_GRIDPOINTS2-1),(2:NO_GRIDPOINTS1-1)),...
        x2(2:(NO_GRIDPOINTS2-1),(2:NO_GRIDPOINTS1-1)),...
        -2.0*epsi*sigma1*dF1_dy-dPot_dy_int_resh,'Linestyle','none')
    title('gradient of opt. tilted potential, y-component')
    xlabel('x')
    ylabel('y')
    zlabel('grad. opt. tilt. pot(x,y), y-component')
    
    dopt_tilt_pot_dx=-2.0*epsi*sigma1*dF1_dx-dPot_dx_int_resh;
        
    figure(10);clf
    plot(x1(1,2:(NO_GRIDPOINTS1-1)),dopt_tilt_pot_dx(1,:),'k-',...
        x1(1,2:(NO_GRIDPOINTS1-1),1),...
        dopt_tilt_pot_dx(ceil((NO_GRIDPOINTS2-2)/2.0),:),'k--')
    xlabel('x')
    xlim([x1(1,2) x1(1,NO_GRIDPOINTS1-1)])
    legend(['a=',num2str(x2(2,1))],...
        ['a=',num2str(x2(ceil((NO_GRIDPOINTS2-2)/2.0),1))])
    ylabel('d/dx opt tilt pot(x,y=a)')
    title('d/dx opt tilt pot(x,y=a)')
    
    
end



