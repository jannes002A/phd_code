%Potential und Kraftfeld Butan
%# codegen

function [f,fx,fy,fz] = KraftButan(K,T,x,y,bound)

% Unterschidelichen Verbindungen 
%      H H H H
%      | | | |
%    H-C-C-C-C-H
%      | | | |
%      H H H H

%      4 7 9  12
%      | | |  |
%    2-1-5-8- 11-13
%      | | |  |
%      3 6 10 14

% Tabellarische Zusammenfassung der unterscheidlichen Potentialle

%     C1 H1 H2 H3 C2 H4 H5 C3 H6 H7 C4 H8 H9 H10
%  C1 -  BW BW BW BW B! B! B! TL TL TL LJ LJ LJ
%  H1 12  - B! B! B! TL TL TL LJ LJ LJ LJ LJ LJ
%  H2 12 21  - B! B! TL TL TL LJ LJ LJ LJ LJ LJ
%  H3 12 21 21  - B! TL TL TL LJ LJ LJ LJ LJ LJ
%  C2 13 23 23 23  - BW BW BW B! B! B! TL TL TL
%  H4 22 31 31 31 12  - B! B! TL TL TL LJ LJ LJ
%  H5 22 31 31 31 12 21  - B! TL TL TL LJ LJ LJ   
%  C3 23 42 42 42 13 22 22  - BW BW BW B! B! B!
%  H6 32 41 41 41 22 31 31 12  - B! B! TL TL TL
%  H7 32 41 41 41 22 31 31 12 21  - B! TL TL TL
%  C4 33 42 42 42 23 42 42 13 22 22  - BW BW BW
%  H8 42 41 41 41 32 43 43 22 31 31 12  - B! B!
%  H9 42 41 41 41 32 43 43 22 31 31 12 21  - B!
% H10 42 41 41 41 32 43 43 22 31 31 12 21 21 -   

% In der unteren Matrix sind alle Wirkungen auf alle anderen Atome
% zusammengestellt. Die erste Zahl gibt an welcher Bindungstyp vorhanden
% ist (1-> BW , 2-> B! (1,3 Bindung), 3 -> TL (1,4 Bildung), 4-> LJ) und
% die zweite Zahl gibt an von welcher Art die Wechselwirkung ist (1-> HH, 
% 2-> CH, 3 -> CC)


%     C1 H1 H2 H3 C2 H4 H5 C3 H6 H7 C4 H8 H9 H10
%  C1 00 12 12 12 13 22 22 23 32 32 33 42 42 42
%  H1 12 00 21 21 23 31 31 42 41 41 42 41 41 41
%  H2 12 21 00 21 23 31 31 42 41 41 42 41 41 41
%  H3 12 21 21 00 23 31 31 42 41 41 42 41 41 41
%  C2 13 23 23 23 00 12 12 13 22 22 23 32 32 32
%  H4 22 31 31 31 12 00 21 22 31 31 42 43 43 43
%  H5 22 31 31 31 12 21 00 22 31 31 42 43 43 43   
%  C3 23 42 42 42 13 22 22 00 12 12 13 22 22 22
%  H6 32 41 41 41 22 31 31 12 00 21 22 31 31 31
%  H7 32 41 41 41 22 31 31 12 21 00 22 31 31 31
%  C4 33 42 42 42 23 42 42 13 22 22 00 12 12 12
%  H8 42 41 41 41 32 43 43 22 31 31 12 00 21 21
%  H9 42 41 41 41 32 43 43 22 31 31 12 21 00 21
% H10 42 41 41 41 32 43 43 22 31 31 12 21 21 00     

% Potential
      % Bindung Winkel CC
%     f = +K*(r)^2 -(1.54)^2)^2 

    % Bindung Winkel CH
%     f = +K*(r)^2 -(1.06)^2)^2 

%     Bindung Winkel 1-3 C C
%     dCC = 1.54^2+1.54^2-2*1.54^2*cos(109.5*(2*pi/360))
%     +K*(r)^2 -(1.54^2+1.54^2-2*1.54^2*cos(109.5*(2*pi/360))))^2
%     Bindung Winkel 1-3 C H
%     dCH = 1.54^2+1.06^2-2*1.54*1.06*cos(109.5*(2*pi/360))
%     +K*(r)^2 -(1.54^2+1.06^2-2*1.06^2*cos(109.5*(2*pi/360))))^2
%     Bindung Winkel 1-3 H H
%     dHH = 1.06^2+1.06^2-2*1.06*1.06*cos(109.5*(2*pi/360))
%     +K*(r)^2 -(1.06^2+1.06^2-2*1.06*1.06*cos(109.5*(2*pi/360))))^2

%     % Lennard Jones CC
%     +0.01*(285800/(r)^2)^6 -372.5/(r)^2)^3)  
%     % Lennard Jones CH
%     +0.01*(37430/(r)^2)^6 - 127.4/(r)^2)^3)   
%     % Lennard Jones HH
%     +0.01*(7220/(r)^2)^6 - 76/(r)^2)^3)   

%     % Torsions 
%     d1CC = d60
%     d1CC=(1.54^2+2*1.54^2*cos(109.5*(2*pi/360))*cos(109.5*(2*pi/360))-1/2*sin(109.5*(2*pi/360))*sin(109.5*(2*pi/360))+1.54^2-2*1.54^2*cos(109.5*(2*pi/360)));
%     d2 CC = d180
%     d2CC=1.54^2+2*1.54^2+1.54^2-2*1.54^2*cos(109.5*(2*pi/360));
%     d1HH = d60
%     d1HH=(1.06^2+2*(1.54*1.06)*cos(109.5*(2*pi/360))*cos(109.5*(2*pi/360))-1/2*sin(109.5*(2*pi/360))*sin(109.5*(2*pi/360))+1.06^2-2*(1.54*1.06)*cos(109.5*(2*pi/360)));
%     d2HH = d180
%     d2HH=1.06^2+2*(1.54*1.06) + 1.06^2 - 2*(1.54+1.06)*cos(109.5*(2*pi/360));
%     d1CH = d60
%     d1CH = (1.06^2+ 2*1.54*1.06*cos(109.5*(2*pi/360))*cos(109.5*(2*pi/360))-1/2*sin(109.5*(2*pi/360))*sin(109.5*(2*pi/360))+1.54^2-2*(1.54*1.06)*cos(109.5*(2*pi/360)));
%     d2CH = d180
%     d2CH = 1.06^2+2*(1.54*1.06)+1.54^2-2*1.54*1.06*cos(109.5*(2*pi/360));
%     + T*((r)^2)-dCH60)^2*((r)^2)-dCH180)^2

%  Übergangeparamter: 1. Abstand der beiden Atome wird als r^2 Übergeben um
%  sich die Wurzel zu sparen
%                     2. Bindungsart
%                     3. Bindungstyp


r = (x(1)-y(1))^2+(x(2)-y(2))^2+(x(3)-y(3))^2;
lj=0.01;

if bound == 12
    % para(2) == 1 wäre BW HH das nicht existiert
    % da es keine direkte Bindung zwischen zwei Wasserstoffatomen gibt
    % BW CH
    f = K*(r - 1.1^2)^2;
elseif bound == 13
    % BW CC
    f = K*(r - 1.54^2)^2;
elseif bound == 21
    % B! (1,3 Wechselwirkung) HH
    f = K*(r - 1.78^2)^2;
elseif bound == 22
    % B! (1,3 Wechselwirkung) CH
    f = K*(r - 2.19^2)^2;
elseif bound == 23
    % B! (1,3 Wechselwirkung) CC
    f = K*(r - 2.6^2)^2;
elseif bound == 31
    %TJ (1,4 Wechslewirkung) HH
    %dHH60 =  1.06^2 + 2*1.54*1.06*((cos(109.5*(2*pi/360))*cos(109.5*(2*pi/360))) - 1/2 * (sin(109.5*(2*pi/360))*sin(109.5*(2*pi/360)))) + 1.06^2 - 2*1.54*1.06*cos(109.5*(2*pi/360));
    %dHH180 = 1.06^2 + 2*1.54*1.06 + 1.06^2 - 2*1.54*1.06*cos(109.5*(2*pi/360));
    % Abstände bestimmt mit der Funktion von Marcus
    % HCCH
    % [d60, d180]=torsion_distance(1.06, 1.54, 1.06) 
    dHH60 = 2.5;
    dHH180 = 3.1;
    f = lj*(7220/r^6 - 76/r^3) +  T*(r-dHH60^2)^2 *(r-dHH180^2)^2;
elseif bound == 32
    % TJ (1,4 Wechslewirkung) CH 
    % dCH60 =  1.1^2 + 2*1.54*1.1*((cos(109.5*(2*pi/360))*cos(109.5*(2*pi/360))) - 1/2*(sin(109.5*(2*pi/360))*sin(109.5*(2*pi/360)))) + 1.54^2 - 2* 1.54 * 1.1*cos(109.5*(2*pi/360));
    % dCH180 = 1.1^2 + 2*1.54*1.1 + 1.54^2 -2*1.54*1.1*cos(109.5*(2*pi/360));
    % Abstände bestimmt mit der Funktion von Marcus
    % CCCH
    % [d60, d180]=torsion_distance(1.54, 1.54, 1.06) 
    dCH60 = 3.58;
    dCH180 = 8.1;
    f = lj*(37430/r^6 - 127.4/r^3) + T*(r-dCH60^2)^2*(r-dCH180^2)^2;
elseif bound == 33
    % TJ (1,4 Wechselwirkungen) CC
    %dCC60 =  1.54^2+2*1.54^2*((cos(109.5*(2*pi/360))*cos(109.5*(2*pi/360))) - 1/2*(sin(109.5*(2*pi/360))*sin(109.5*(2*pi/360))) ) + 1.54^2 - 2*1.54^2*cos(109.5*(2*pi/360));
    %dCC180 = 1.54^2 + 2*1.54^2 + 1.54^2 - 2*1.54^2 * cos(109.5*(2*pi/360));
    % CCCC
    % [d60, d180]=torsion_distance(1.54, 1.54, 1.54)
    dCC60 = 4.74;
    dCC180 = 11.1;
    f = lj*(285800/r^6 -372.5/r^3) + T*(r-dCC60^2)^2*(r-dCC180^2)^2; 
elseif bound == 41
    % LJ HH
    f = (7220/r^6 - 76/r^3);
elseif bound == 42
    % LJ CH
    f = (37430/r^6 - 127.4/r^3);
elseif bound == 43
    % LJ CC
    f=(285800/r^6 - 372.5/r^3); 
else
    f = 0;
end
% Ableitung

if bound == 12
    % para(2) == 1 wäre BW HH das nicht existiert
    % da es keine direkte Bindung zwischen zwei Wasserstoffatomen gibt
    % BW CH
    fx = 4*(x(1)-y(1))*K*(r - 1.1^2);
    fy = 4*(x(2)-y(2))*K*(r - 1.1^2);
    fz = 4*(x(3)-y(3))*K*(r - 1.1^2);
elseif bound == 13
    % BW CC
    fx = 4*(x(1)-y(1))*K*(r - 1.54^2);
    fy = 4*(x(2)-y(2))*K*(r - 1.54^2);
    fz = 4*(x(3)-y(3))*K*(r - 1.54^2);
elseif bound == 21
    % B! (1,3 Wechselwirkung) HH
    fx = 4*(x(1)-y(1))*K*(r - 1.78^2);
    fy = 4*(x(2)-y(2))*K*(r - 1.78^2);
    fz = 4*(x(3)-y(3))*K*(r - 1.78^2);
elseif bound == 22
    % B! (1,3 Wechselwirkung) CH
    fx = 4*(x(1)-y(1))*K*(r -2.19^2);
    fy = 4*(x(2)-y(2))*K*(r -2.19^2);
    fz = 4*(x(3)-y(3))*K*(r -2.19^2);
elseif bound == 23
    % B! (1,3 Wechselwirkung) CC
    fx = 4*(x(1)-y(1))*K*(r - 2.6^2);
    fy = 4*(x(2)-y(2))*K*(r - 2.6^2);
    fz = 4*(x(3)-y(3))*K*(r - 2.6^2);
elseif bound == 31
    %TJ (1,4 Wechslewirkung) HH
    dHH60 = 2.5;
    dHH180 = 3.1;
    fx = lj*( (6*76*(x(1)-y(1)))/r^4  -  (12*7220*(x(1)-y(1)))/r^7 ) +  4*T*(x(1)-y(1))*((r-dHH60^2) * (r-dHH180^2)^2 + ((r-dHH60^2)^2 *(r-dHH180^2)));
    fy = lj*( (6*76*(x(2)-y(2)))/r^4  -  (12*7220*(x(2)-y(2)))/r^7 ) +  4*T*(x(2)-y(2))*((r-dHH60^2) * (r-dHH180^2)^2 + ((r-dHH60^2)^2 *(r-dHH180^2)));
    fz = lj*( (6*76*(x(3)-y(3)))/r^4  -  (12*7220*(x(3)-y(3)))/r^7 ) +  4*T*(x(3)-y(3))*((r-dHH60^2) * (r-dHH180^2)^2 + ((r-dHH60^2)^2 *(r-dHH180^2)));
elseif bound == 32
    %TJ (1,4 Wechslewirkung) CH 
    dCH60 = 3.58;
    dCH180 = 8.1;
    fx = lj*((6*127.4*(x(1)-y(1)))/r^4 - (12*37430*(x(1)-y(1)))/r^7 ) + 4*T*(x(1)-y(1))*((r-dCH60^2) * (r-dCH180^2)^2 + ((r-dCH60^2)^2 *(r-dCH180^2)));
    fy = lj*((6*127.4*(x(2)-y(2)))/r^4 - (12*37430*(x(2)-y(2)))/r^7 ) + 4*T*(x(2)-y(2))*((r-dCH60^2) * (r-dCH180^2)^2 + ((r-dCH60^2)^2 *(r-dCH180^2)));
    fz = lj*((6*127.4*(x(3)-y(3)))/r^4 - (12*37430*(x(3)-y(3)))/r^7 ) + 4*T*(x(3)-y(3))*((r-dCH60^2) * (r-dCH180^2)^2 + ((r-dCH60^2)^2 *(r-dCH180^2)));
elseif bound == 33
   % TJ (1,4 Wechselwirkungen) CC
   dCC60 = 4.74;
   dCC180 = 11.1;
   fx = lj*((6*372.5*(x(1)-y(1)))/r^4-(12*285800*(x(1)-y(1)))/r^7) + 4*T*(x(1)-y(1))*((r-dCC60^2) * (r-dCC180^2)^2 + ((r-dCC60^2)^2 *(r-dCC180^2)));
   fy = lj*((6*372.5*(x(2)-y(2)))/r^4-(12*285800*(x(2)-y(2)))/r^7) + 4*T*(x(2)-y(2))*((r-dCC60^2) * (r-dCC180^2)^2 + ((r-dCC60^2)^2 *(r-dCC180^2)));
   fz = lj*((6*372.5*(x(3)-y(3)))/r^4-(12*285800*(x(3)-y(3)))/r^7) + 4*T*(x(3)-y(3))*((r-dCC60^2) * (r-dCC180^2)^2 + ((r-dCC60^2)^2 *(r-dCC180^2)));
elseif bound == 41
   % LJ HH
   fx = (6*76*(x(1)-y(1)))/r^4  -  (12*7220*(x(1)-y(1)))/r^7;
   fy = (6*76*(x(2)-y(2)))/r^4  -  (12*7220*(x(2)-y(2)))/r^7;
   fz = (6*76*(x(3)-y(3)))/r^4  -  (12*7220*(x(3)-y(3)))/r^7;
elseif bound == 42
   % LJ CH
   fx = (6*127.4*(x(1)-y(1)))/r^4 - (12*37430*(x(1)-y(1)))/r^7;
   fy = (6*127.4*(x(2)-y(2)))/r^4 - (12*37430*(x(2)-y(2)))/r^7;
   fz = (6*127.4*(x(3)-y(3)))/r^4 - (12*37430*(x(3)-y(3)))/r^7;
elseif bound == 43
   % LJ CC
   fx = (6*372.5*(x(1)-y(1)))/r^4-(12*285800*(x(1)-y(1)))/r^7;
   fy = (6*372.5*(x(2)-y(2)))/r^4-(12*285800*(x(2)-y(2)))/r^7;
   fz = (6*372.5*(x(3)-y(3)))/r^4-(12*285800*(x(3)-y(3)))/r^7;
else
    fx = 0;
    fy = 0;
    fz = 0;
end