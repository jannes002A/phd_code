% dVbais 2nd derivative

function [dvbias]=dVbias2nd(x,sc)

dvbias=zeros(1,length(x));

for i=1:length(x)
    if (x(i)<=0)
        dvbias(i) = sc*(6.*x(i).^2-2 +2);
    else
        dvbias(i)=0;
    end
end