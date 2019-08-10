function grad = DV(q,POTENTIAL)
%Me: this function computes the gradient of a specific potential given by
%POTENTIAL. The other input 'q' refers to the point (configuration) at
%which the Laplacian of the potential is computed.
 switch POTENTIAL
    case 'sym1wellPot'; % harm. Oszillator
        grad = q;
    case 'sym2wellPot'; % symmetric double well potential
        grad = 4*q.*(q.*q-1);
    case 'asym2wellPot';  % asymmetric double well potential    
        grad =  4*q.*(q.*q-1) - 0.2;
    case 'asym3wellPot'; % asymmetric three well potential
        grad = ( 3*q.^5 - 60*q.^3 + 238*q + 28) / 200; 
    case 'asym5wellPot'; % asymmetric three well potential
        grad = ( 3*q.^5 - 60*q.^3 + 238*q + 28) / 200 + 0.6*(q+2).*exp(-0.5*(q+2).^2/(0.2)^2)/(0.2)^2 + 0.7*(q-1.8).*exp(-0.5*(q-1.8).^2/(0.2)^2)/(0.2)^2; 
    case 'rugged5wellPot'; % asymmetric three well potential
        grad = ( 3*q.^5 - 60*q.^3 + 238*q + 28) / 200 + 0.6*(q+2).*exp(-0.5*(q+2).^2/(0.2)^2)/(0.2)^2 + 0.7*(q-1.8).*exp(-0.5*(q-1.8).^2/(0.2)^2)/(0.2)^2+ 0.05*20*cos(20*q) + 0.025*15*cos(15*q);      
     case 'diffusive'; % diffusive two well potential
        a = 5.62;
        x=q;
        grad = 4*(x.^2-1).*x;
        I=find(x>a);
        grad(I) = 4*((x(I)-a).^2-1).*(x(I)-a);
        m=4; epsi=0.2;
        b = a/(2*pi*m);
        J=find((x>=0)&(x<=a));
        grad(J) = -epsi * sin(x(J)/b)/b;
        
    case 'paper3wellPot'; % paper three well potential
        c = 0.1;
        grad = 2*( (q.^2-1).^2 -1 + c*q ).*( 4*q.*(q.^2-1) + c);
    otherwise;
	fprintf('\n\n Unkown potential! \n\n'); return
 end;

 