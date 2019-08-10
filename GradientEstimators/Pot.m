function pot_energie=Pot(q,POTENTIAL)
%Me : original was pot_energie=V(q,POTENTIAL)
%%% Note: q has to be either a 1x? or 
%%% for the 'eric' potential a 2x? matrix.

 switch POTENTIAL
    case 'sym1wellPot'; % harm. Oszillator
        pot_energie = 0.5*q.*q;
    case 'sym2wellPot'; % symmetric double well potential
        pot_energie = (q.*q-1).^2;
    case 'asym2wellPot';  % asymmetric double well potential    
        pot_energie = (q.*q-1).^2 - 0.2*q + 0.3;
    case 'asym3wellPot'; % asymmetric three well potential
        pot_energie = ( 0.5*q.^6 - 15*q.^4 + 119*q.^2 + 28*q + 50) / 200;
    case 'asym5wellPot'; % asymmetric three well potential
        pot_energie = ( 0.5*q.^6 - 15*q.^4 + 119*q.^2 + 28*q + 50) / 200 - 0.6*exp(-0.5*(q+2).^2/(0.2)^2)- 0.7*exp(-0.5*(q-1.8).^2/(0.2)^2);
    case 'rugged5wellPot'; % asymmetric three well potential
        pot_energie = ( 0.5*q.^6 - 15*q.^4 + 119*q.^2 + 28*q + 50) / 200 - 0.6*exp(-0.5*(q+2).^2/(0.2)^2)- 0.7*exp(-0.5*(q-1.8).^2/(0.2)^2) + 0.05*sin(20*q) + 0.025*sin(15*q);
    case 'paper3wellPot'; % paper three well potential
        c = 0.1;
        pot_energie = ( (q.^2-1).^2 -1 + c*q ).^2;
    case 'diffusive'; % two well diffusive potential
        c = 0.1;
        a = 5.62;
        x=q;
        pot_energie = (x.^2-1).^2;
        I=find(x>a);
        pot_energie(I) = ((x(I)-a).^2-1).^2;
        m=4; epsi=0.2;
        b = a/(2*pi*m);
        J=find((x>=0)&(x<=a));
        pot_energie(J) = epsi * cos(x(J)/b)+1-epsi;

    case 'eric4wellPot'; % Eric Barth's 2d example
        pot_energie = ( (q(2,:)+.4).^4 - 20*q(2,:).^2 + ...  
		      (q(1,:)-.1).^4 - 20*q(1,:).^2 + ...
		       q(2,:).*q(1,:) + 290.4 + ...
		      10*cos(5*q(2,:)).*sin(5*q(1,:)) )/100;
    otherwise;
	error('\n\n Unkown potential! \n\n');
 end;
