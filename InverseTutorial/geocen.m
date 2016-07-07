function [ gcolat ] = geocen( colat )
%GEOCEN Corrects geographic latitude to geocentric

fac=0.993305621334896;
gcolat=0.5*pi-atan(fac*cos(colat)/max(1.e-200,sin(colat)));
return

end

