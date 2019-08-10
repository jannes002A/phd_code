function laplace = DV(q,POTENTIAL)
%Me: this function computes the laplacian of a specific potential given by
%POTENTIAL. The other input 'q' refers to the point (configuration) at
%which the Laplacian of the potential is computed.
 switch POTENTIAL
    case 'sym1wellPot'; % harm. Oszillator
        laplace = 1;
    case 'sym2wellPot'; % symmetric double well potential
        laplace = 12*q.^2 - 4;
    case 'asym2wellPot';  % asymmetric double well potential    
        laplace =  12*q.^2 - 4;
    case 'asym3wellPot'; % asymmetric three well potential
        laplace = ( 15*q.^4 - 180*q.^2 + 238) / 200;
    case 'asym5wellPot'; % asymmetric three well potential
        laplace = ( 15*q.^4 - 180*q.^2 + 238) / 200;
    case 'rugged5wellPot'; % asymmetric three well potential
        laplace = ( 15*q.^4 - 180*q.^2 + 238) / 200;
    case 'diffusive'; % diffusive two well potential
        a = 5.62;
        x=q;
        laplace = 4*(x.^2-1) + 8*x.^2;
        I=find(x>a);
        laplace(I) = 4*((x(I)-a).^2-1) + 8*(x(I)-a).^2;
        m=4; epsi=0.2;
        b = a/(2*pi*m);
        J=find((x>=0)&(x<=a));
        laplace(J) = -epsi*epsi * cos(x(J)/b)/b/b;
 
    case 'paper3wellPot'; % paper three well potential
        c = 0.1;
        laplace = 2*(4*q.*(q.^2 - 1) + c ).^2 + 2*((q.^2-1).^2 -1 + c*q ).*(12*q.^2 - 4);
    otherwise;
	fprintf('\n\n Unkown potential! \n\n'); return
 end;

