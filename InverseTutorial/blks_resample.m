function [ lat,lon,mz ] = blks_resample(nblk,bsize,nlat,mlat,hsize,mest,sampling)
%BLKS_RESAMPLE Returns lat lon vectors for a model parameterization
%  input values nblk,bsize,nlat,mlat,hsize are return values from
%  blks2d.  output is resampled to have even bsize spacing in both latitude
%  and longitude
clear indx nlon colat lat lon

nlat_sample=round(180/sampling);
sampsize=180/nlat_sample;

nlon=2*nlat_sample;
colat=repmat((0.5*sampsize:sampsize:180)',1,nlon);
lat=90.0-colat;
lon=repmat((0.5*sampsize:sampsize:360),nlat_sample,1);
indx=round(fblk(colat,lon,nlat,bsize,mlat,hsize));
mz=mest(indx);

return
end

