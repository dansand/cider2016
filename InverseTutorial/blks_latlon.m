function [ lat,lon ] = blks_latlon(nblk,bsize,nlat,mlat,hsize)
%BLKS_LATLON Returns lat lon vectors for a model parameterization
%  input values nblk,bsize,nlat,mlat,hsize are return values from
%  blks2d

lat=zeros(nblk,1);
lon=zeros(nblk,1);

ib=0;
for i=1:nlat
lati=90.0-((i-0.5)*bsize);
    for j=1:(mlat(i+1)-mlat(i))
        ib=ib+1;
        lonj=(j-0.5)*hsize(i);
        lat(ib)=lati;
        lon(ib)=lonj;
    end
end

return
end

