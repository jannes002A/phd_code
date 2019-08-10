function [ges_F, grad] = fast_pot_butan_hom(atom,lambda)

% Die Koodianten der Atome müssen als Zeilenvektor übergeben werden

% Wechelwirkungen von Butan
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

WW = [00 12 12 12 13 22 22 23 32 32 33 42 42 42;
      12 00 21 21 23 31 31 42 41 41 42 41 41 41;
      12 21 00 21 23 31 31 42 41 41 42 41 41 41;
      12 21 21 00 23 31 31 42 41 41 42 41 41 41;
      13 23 23 23 00 12 12 13 22 22 23 32 32 32;
      22 31 31 31 12 00 21 22 31 31 42 43 43 43;
      22 31 31 31 12 21 00 22 31 31 42 43 43 43;   
      23 42 42 42 13 22 22 00 12 12 13 22 22 22;
      32 41 41 41 22 31 31 12 00 21 22 31 31 31;
      32 41 41 41 22 31 31 12 21 00 22 31 31 31;
      33 42 42 42 23 42 42 13 22 22 00 12 12 12;
      42 41 41 41 32 43 43 22 31 31 12 00 21 21;
      42 41 41 41 32 43 43 22 31 31 12 21 00 21;
      42 41 41 41 32 43 43 22 31 31 12 21 21 00];

grad = zeros(3,14); 

[~,m]=size(atom);

ges_F = 0; 
K=5;
T=5;


rx = (x(1)-y(1))^2+(x(2)-y(2))^2+(x(3)-y(3))^2;
lj=0.01;

%Vektor mit dem Stützstellen


for i=1:m
    for j = i+1:m
        %[f,fx,fy,fz] = KraftButan_hom(K,T,atom(:,i),atom(:,j),WW(i,j),0);
        %[f,fx,fy,fz] = KraftButan(K,T,atom(:,i),atom(:,j),WW(i,j));
        
        x=atom(:,i);
        y=atom(:,j);
        
        bound=WW(i,j);
        
       
        con = 1/(sqrt(pi)*rx);
        lh = sqrt(2)*lambda;
        s= [-lh*2.02018287,-lh*0.95857246,0,lh*0.95857246,lh*2.02018287];
        r = rx+s;
        w = [0.01995324,0.39361932,0.94530872,0.39361932,0.01995324];

        abl_r1=4*(x(1)-y(1));
        abl_r2=4*(x(2)-y(2));
        abl_r3=4*(x(3)-y(3));

        switch bound
            case 12
            % para(2) == 1 wäre BW HH das nicht existiert
            % da es keine direkte Bindung zwischen zwei Wasserstoffatomen gibt
            % BW CH
            f = K.*(r - 1.1^2).^2;
            f_hom = con*(w.*r*f');
            %Ableitung 


            H1 = 4.*(x(1)-y(1)).*K.*(r - 1.1^2);
            H2 = 4.*(x(2)-y(2)).*K.*(r - 1.1^2);
            H3 = 4.*(x(3)-y(3)).*K.*(r - 1.1^2);

            fx_hom = -abl_r1/(sqrt(pi).*rx^2)* w.*r*f' + 1/(sqrt(pi)*rx)*(w.*(abl_r1)*f'+ w.*r*H1');
            fy_hom = -abl_r2/(sqrt(pi).*rx^2)* w.*r*f' + 1/(sqrt(pi)*rx)*(w.*(abl_r2)*f'+ w.*r*H2');
            fz_hom = -abl_r3/(sqrt(pi).*rx^2)* w.*r*f' + 1/(sqrt(pi)*rx)*(w.*(abl_r3)*f'+ w.*r*H3');
            case  13
            % BW CC
            f = K.*(r - 1.54^2).^2;
            f_hom = con*(w.*r*f');

            %Ableitung


            H1 = 4*(x(1)-y(1))*K*(r - 1.54^2);
            H2 = 4*(x(2)-y(2))*K*(r - 1.54^2);
            H3 = 4*(x(3)-y(3))*K*(r - 1.54^2);

            fx_hom = -abl_r1/(sqrt(pi).*rx^2)* w.*r*f' + 1/(sqrt(pi)*rx)*(w.*(abl_r1)*f'+ w.*r*H1');
            fy_hom = -abl_r2/(sqrt(pi).*rx^2)* w.*r*f' + 1/(sqrt(pi)*rx)*(w.*(abl_r2)*f'+ w.*r*H2');
            fz_hom = -abl_r3/(sqrt(pi).*rx^2)* w.*r*f' + 1/(sqrt(pi)*rx)*(w.*(abl_r3)*f'+ w.*r*H3');
            case 21
            % B! (1,3 Wechselwirkung) HH
            f = K.*(r - 1.78^2).^2;
            f_hom = con*(w.*r*f');

            %Ableitung


            H1 = 4*(x(1)-y(1))*K*(r - 1.78^2);
            H2 = 4*(x(2)-y(2))*K*(r - 1.78^2);
            H3 = 4*(x(3)-y(3))*K*(r - 1.78^2);

            fx_hom = -abl_r1/(sqrt(pi).*rx^2)* w.*r*f' + 1/(sqrt(pi)*rx)*(w.*(abl_r1)*f'+ w.*r*H1');
            fy_hom = -abl_r2/(sqrt(pi).*rx^2)* w.*r*f' + 1/(sqrt(pi)*rx)*(w.*(abl_r2)*f'+ w.*r*H2');
            fz_hom = -abl_r3/(sqrt(pi).*rx^2)* w.*r*f' + 1/(sqrt(pi)*rx)*(w.*(abl_r3)*f'+ w.*r*H3');
            case 22
            % B! (1,3 Wechselwirkung) CH
            f = K.*(r - 2.19^2).^2;
            f_hom = con*(w.*r*f');

            %Ableitung
            H1 = 4*(x(1)-y(1))*K*(r -2.19^2);
            H2 = 4*(x(2)-y(2))*K*(r -2.19^2);
            H3 = 4*(x(3)-y(3))*K*(r -2.19^2);

            fx_hom = -abl_r1/(sqrt(pi).*rx^2)* w.*r*f' + 1/(sqrt(pi)*rx)*(w.*(abl_r1)*f'+ w.*r*H1');
            fy_hom = -abl_r2/(sqrt(pi).*rx^2)* w.*r*f' + 1/(sqrt(pi)*rx)*(w.*(abl_r2)*f'+ w.*r*H2');
            fz_hom = -abl_r3/(sqrt(pi).*rx^2)* w.*r*f' + 1/(sqrt(pi)*rx)*(w.*(abl_r3)*f'+ w.*r*H3');
            case 23
            % B! (1,3 Wechselwirkung) CC
            f = K.*(r - 2.6^2).^2;
            f_hom = con*(w.*r*f'); 

            %Abelitung
            H1 = 4*(x(1)-y(1))*K*(r - 2.6^2);
            H2 = 4*(x(2)-y(2))*K*(r - 2.6^2);
            H3 = 4*(x(3)-y(3))*K*(r - 2.6^2);

            fx_hom = -abl_r1/(sqrt(pi).*rx^2)* w.*r*f' + 1/(sqrt(pi)*rx)*(w.*(abl_r1)*f'+ w.*r*H1');
            fy_hom = -abl_r2/(sqrt(pi).*rx^2)* w.*r*f' + 1/(sqrt(pi)*rx)*(w.*(abl_r2)*f'+ w.*r*H2');
            fz_hom = -abl_r3/(sqrt(pi).*rx^2)* w.*r*f' + 1/(sqrt(pi)*rx)*(w.*(abl_r3)*f'+ w.*r*H3');
            case 31
                     %TJ (1,4 Wechslewirkung) HH
            %dHH60 =  1.06^2 + 2*1.54*1.06*((cos(109.5*(2*pi/360))*cos(109.5*(2*pi/360))) - 1/2 * (sin(109.5*(2*pi/360))*sin(109.5*(2*pi/360)))) + 1.06^2 - 2*1.54*1.06*cos(109.5*(2*pi/360));
            %dHH180 = 1.06^2 + 2*1.54*1.06 + 1.06^2 - 2*1.54*1.06*cos(109.5*(2*pi/360));
            % Abstände bestimmt mit der Funktion von Marcus
            % HCCH
            % [d60, d180]=torsion_distance(1.06, 1.54, 1.06) 
            dHH60 = 2.5;
            dHH180 = 3.1;
            f = lj.*(7220./r.^6 - 76./r.^3) +  T.*(r-dHH60^2).^2 .*(r-dHH180^2).^2;
            f_hom = con*(w.*r*f');

            H1 = lj*( (6*76*(x(1)-y(1)))./r.^4  -  (12*7220*(x(1)-y(1)))./r.^7 ) +  4*T*(x(1)-y(1)).*((r-dHH60^2) .* (r-dHH180.^2).^2 + ((r-dHH60^2).^2 .*(r-dHH180^2)));
            H2 = lj*( (6*76*(x(2)-y(2)))./r.^4  -  (12*7220*(x(2)-y(2)))./r.^7 ) +  4*T*(x(2)-y(2)).*((r-dHH60^2) .* (r-dHH180.^2).^2 + ((r-dHH60^2).^2 .*(r-dHH180^2)));
            H3 = lj*( (6*76*(x(3)-y(3)))./r.^4  -  (12*7220*(x(3)-y(3)))./r.^7 ) +  4*T*(x(3)-y(3)).*((r-dHH60^2) .* (r-dHH180.^2).^2 + ((r-dHH60^2).^2 .*(r-dHH180^2)));

            fx_hom = -abl_r1/(sqrt(pi).*rx^2)* w.*r*f' + 1/(sqrt(pi)*rx)*(w.*(abl_r1)*f'+ w.*r*H1');
            fy_hom = -abl_r2/(sqrt(pi).*rx^2)* w.*r*f' + 1/(sqrt(pi)*rx)*(w.*(abl_r2)*f'+ w.*r*H2');
            fz_hom = -abl_r3/(sqrt(pi).*rx^2)* w.*r*f' + 1/(sqrt(pi)*rx)*(w.*(abl_r3)*f'+ w.*r*H3');
            case 32
            % TJ (1,4 Wechslewirkung) CH 
            % dCH60 =  1.1^2 + 2*1.54*1.1*((cos(109.5*(2*pi/360))*cos(109.5*(2*pi/360))) - 1/2*(sin(109.5*(2*pi/360))*sin(109.5*(2*pi/360)))) + 1.54^2 - 2* 1.54 * 1.1*cos(109.5*(2*pi/360));
            % dCH180 = 1.1^2 + 2*1.54*1.1 + 1.54^2 -2*1.54*1.1*cos(109.5*(2*pi/360));
            % Abstände bestimmt mit der Funktion von Marcus
            % CCCH
            % [d60, d180]=torsion_distance(1.54, 1.54, 1.06) 
            dCH60 = 3.58;
            dCH180 = 8.1;
            f = lj.*(37430./r.^6 - 127.4./r.^3) + T.*(r-dCH60^2).^2.*(r-dCH180^2).^2;
            f_hom = con*(w.*r*f');

            %Ableitung
            H1 = lj*((6*127.4*(x(1)-y(1)))./r.^4 - (12*37430*(x(1)-y(1)))./r.^7 ) + 4*T*(x(1)-y(1)).*((r-dCH60^2) .* (r-dCH180^2).^2 + ((r-dCH60^2).^2 .*(r-dCH180^2)));
            H2 = lj*((6*127.4*(x(2)-y(2)))./r.^4 - (12*37430*(x(2)-y(2)))./r.^7 ) + 4*T*(x(2)-y(2)).*((r-dCH60^2) .* (r-dCH180^2).^2 + ((r-dCH60^2).^2 .*(r-dCH180^2)));
            H3 = lj*((6*127.4*(x(3)-y(3)))./r.^4 - (12*37430*(x(3)-y(3)))./r.^7 ) + 4*T*(x(3)-y(3)).*((r-dCH60^2) .* (r-dCH180^2).^2 + ((r-dCH60^2).^2 .*(r-dCH180^2)));

            fx_hom = -abl_r1/(sqrt(pi).*rx^2)* w.*r*f' + 1/(sqrt(pi)*rx)*(w.*(abl_r1)*f'+ w.*r*H1');
            fy_hom = -abl_r2/(sqrt(pi).*rx^2)* w.*r*f' + 1/(sqrt(pi)*rx)*(w.*(abl_r2)*f'+ w.*r*H2');
            fz_hom = -abl_r3/(sqrt(pi).*rx^2)* w.*r*f' + 1/(sqrt(pi)*rx)*(w.*(abl_r3)*f'+ w.*r*H3');

            case 33
                       % TJ (1,4 Wechselwirkungen) CC
            %dCC60 =  1.54^2+2*1.54^2*((cos(109.5*(2*pi/360))*cos(109.5*(2*pi/360))) - 1/2*(sin(109.5*(2*pi/360))*sin(109.5*(2*pi/360))) ) + 1.54^2 - 2*1.54^2*cos(109.5*(2*pi/360));
            %dCC180 = 1.54^2 + 2*1.54^2 + 1.54^2 - 2*1.54^2 * cos(109.5*(2*pi/360));
            % CCCC
            % [d60, d180]=torsion_distance(1.54, 1.54, 1.54)
            dCC60 = 4.74;
            dCC180 = 11.1;
            f = lj.*(285800./r.^6 -372.5./r.^3) + T.*(r-dCC60^2).^2.*(r-dCC180^2).^2; 
            f_hom = con*(w.*r*f');

            H1 = lj*((6*372.5*(x(1)-y(1)))./r.^4-(12*285800*(x(1)-y(1)))./r.^7) + 4*T*(x(1)-y(1)).*((r-dCC60^2) .* (r-dCC180^2).^2 + ((r-dCC60^2).^2 .*(r-dCC180^2)));
            H2 = lj*((6*372.5*(x(2)-y(2)))./r.^4-(12*285800*(x(2)-y(2)))./r.^7) + 4*T*(x(2)-y(2)).*((r-dCC60^2) .* (r-dCC180^2).^2 + ((r-dCC60^2).^2 .*(r-dCC180^2)));
            H3 = lj*((6*372.5*(x(3)-y(3)))./r.^4-(12*285800*(x(3)-y(3)))./r.^7) + 4*T*(x(3)-y(3)).*((r-dCC60^2) .* (r-dCC180^2).^2 + ((r-dCC60^2).^2 .*(r-dCC180^2)));

            fx_hom = -abl_r1/(sqrt(pi).*rx^2)* w.*r*f' + 1/(sqrt(pi)*rx)*(w.*(abl_r1)*f'+ w.*r*H1');
            fy_hom = -abl_r2/(sqrt(pi).*rx^2)* w.*r*f' + 1/(sqrt(pi)*rx)*(w.*(abl_r2)*f'+ w.*r*H2');
            fz_hom = -abl_r3/(sqrt(pi).*rx^2)* w.*r*f' + 1/(sqrt(pi)*rx)*(w.*(abl_r3)*f'+ w.*r*H3');
            case 41
            % LJ HH
            f = (7220./r.^6 - 76./r.^3);
            f_hom = con*(w.*r*f');

            H1 = (6*76*(x(1)-y(1)))./r.^4  -  (12*7220*(x(1)-y(1)))./r.^7;
            H2 = (6*76*(x(2)-y(2)))./r.^4  -  (12*7220*(x(2)-y(2)))./r.^7;
            H3 = (6*76*(x(3)-y(3)))./r.^4  -  (12*7220*(x(3)-y(3)))./r.^7;

            fx_hom = -abl_r1/(sqrt(pi).*rx^2)* w.*r*f' + 1/(sqrt(pi)*rx)*(w.*(abl_r1)*f'+ w.*r*H1');
            fy_hom = -abl_r2/(sqrt(pi).*rx^2)* w.*r*f' + 1/(sqrt(pi)*rx)*(w.*(abl_r2)*f'+ w.*r*H2');
            fz_hom = -abl_r3/(sqrt(pi).*rx^2)* w.*r*f' + 1/(sqrt(pi)*rx)*(w.*(abl_r3)*f'+ w.*r*H3');
            case 42
            % LJ CH
            f = (37430./r.^6 - 127.4./r.^3);
            f_hom = con*(w.*r*f');

            H1 = (6*127.4*(x(1)-y(1)))./r.^4 - (12*37430*(x(1)-y(1)))./r.^7;
            H2 = (6*127.4*(x(2)-y(2)))./r.^4 - (12*37430*(x(2)-y(2)))./r.^7;
            H3 = (6*127.4*(x(3)-y(3)))./r.^4 - (12*37430*(x(3)-y(3)))./r.^7;

            fx_hom = -abl_r1/(sqrt(pi).*rx^2)* w.*r*f' + 1/(sqrt(pi)*rx)*(w.*(abl_r1)*f'+ w.*r*H1');
            fy_hom = -abl_r2/(sqrt(pi).*rx^2)* w.*r*f' + 1/(sqrt(pi)*rx)*(w.*(abl_r2)*f'+ w.*r*H2');
            fz_hom = -abl_r3/(sqrt(pi).*rx^2)* w.*r*f' + 1/(sqrt(pi)*rx)*(w.*(abl_r3)*f'+ w.*r*H3');

            case 43
            f=(285800./r.^6 - 372.5./r.^3); 
            f_hom = con*(w.*r*f'); 

            H1 = (6*372.5*(x(1)-y(1)))./r.^4-(12*285800*(x(1)-y(1)))./r.^7;
            H2 = (6*372.5*(x(2)-y(2)))./r.^4-(12*285800*(x(2)-y(2)))./r.^7;
            H3 = (6*372.5*(x(3)-y(3)))./r.^4-(12*285800*(x(3)-y(3)))./r.^7;

            fx_hom = -abl_r1/(sqrt(pi).*rx^2)* w.*r*f' + 1/(sqrt(pi)*rx)*(w.*(abl_r1)*f'+ w.*r*H1');
            fy_hom = -abl_r2/(sqrt(pi).*rx^2)* w.*r*f' + 1/(sqrt(pi)*rx)*(w.*(abl_r2)*f'+ w.*r*H2');
            fz_hom = -abl_r3/(sqrt(pi).*rx^2)* w.*r*f' + 1/(sqrt(pi)*rx)*(w.*(abl_r3)*f'+ w.*r*H3');
            otherwise
            f = 0;

            fx = 0;
            fy = 0;
            fz = 0;
        end
        ges_F = ges_F + f_hom;
        
        grad(1,i) = grad(1,i) + fx_hom;
        grad(2,i) = grad(2,i) + fy_hom;
        grad(3,i) = grad(3,i) + fz_hom;
        
        grad(1,j) = grad(1,j) - fx_hom;
        grad(2,j) = grad(2,j) - fy_hom;
        grad(3,j) = grad(3,j) - fz_hom;
    end
end

grad =  grad(:);