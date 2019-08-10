% function zur Bestimmung des Torsionswinkels

function [ang]=torsionAngButan(points)
% die nehem an, dass die Position der Atome als Matix eingegeben wird der
% Größe (4 x 3). Die Atome sind nach ihrem Auftreten nach geordnet

% Vektor aus C(2)-C(1)
x1 = points(:,2)-points(:,1);
% Vektor aus C(2)-C(3) = x2 => C(3)-C(2)=-x2
x2 = points(:,2)-points(:,3);
x3 = points(:,3)-points(:,4);

cros1 = cross(x1,x2);
% nehem hier -x2 da man eigentlich 
cros2 = cross(-x2,x3);

%Skalarprodukt zwischen den beiden Kreuzprodukten
ang = cros1'*cros2/(norm(cros1)*norm(cros2));

ang=acos(ang)*(180/(pi));