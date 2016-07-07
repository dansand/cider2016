function [ indx ] = fblk( t,p,nlat,bsize,mlat,hsize )
%FBLK Returns the index of block containing point (t,p) using grid info
%   output from blks2d (nlat,bsize,mlat,hsize)
%   t and p are permitted to be vectors

    it=floor(t/bsize)+1;
    [I,J]=find(it<1|it>nlat);
    it(I,J)=1;
    j1=mlat(it)+1;
    indx=floor(p./hsize(it))+j1;
    indx(I,J)=-1;
    %if(it>=1&&it<=nlat)
    %    j1=mlat(it)+1;
    %    indx=floor(p/hsize(it))+j1;
    %else
    %    indx=-1;
    %end
    return
end

