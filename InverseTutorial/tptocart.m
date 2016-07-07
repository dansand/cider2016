function [ cart ] = tptocart( theta,phi )
%TPTOCART Converts theta phi coordinates on a unit sphere to Cartesian
%   Simple spherical to Cartesian conversion
    s=sin(theta);
    cart(1)=s*cos(phi);
    cart(2)=s*sin(phi);
    cart(3)=cos(theta);
    return
end

