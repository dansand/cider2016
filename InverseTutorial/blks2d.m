function [ nblk,bsize,nlat,mlat,hsize ] = blks2d(blocksz)
%BLKS2D Establishes block parameters
%  Divides a spherical surface into approximately blocksz x blocksz blocks
%  where blocksz is degrees measured at the equator
%  input
%    blocksz - the desired block size (in degrees)
%  outputs
%    nblk - total number of blocks in the model
%    bsize - the nominal output block size in degrees (may differ to be a
%            factor of 180)
%    nlat - the number of latitudinal samples
%    mlat - the block index limits for each latitudinal band of blocks
%    hsize - the dimension of each block in longitudinal degrees for each
%            latitudinal band of blocks

rconv=180./pi;

nlat=round(180./blocksz);
bsize=180./nlat;

ib=0;
tmax=0.0;
mlat=zeros(1,nlat+1);
hsize=zeros(1,nlat);
mlat(1)=0.0;
for i=1:nlat
    tmin=tmax;
    tmax=tmin+bsize;
    th=0.5*(tmin+tmax)/rconv;
    s1=sin(th);
    mlon=max(round((360./bsize)*s1),1);
    hsize(i)=360./mlon;
    ib=ib+mlon;
    mlat(i+1)=ib;
end
nblk=ib;
return
end

