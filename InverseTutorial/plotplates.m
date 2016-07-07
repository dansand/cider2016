function [ error ] = plotplates( lonmax )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
error=-1;
lat=[]; lon=[];

filename='nuvel_1_plates.txt'; %should probably not have this hard-coded
tolerance=0.9;

fid=fopen(filename,'r');
tline=fgetl(fid); %read first line
npts=0;
while ischar(tline)
    if(tline(1)==':') %segment boundary
        if(isempty(lat))
            tline=fgetl(fid);
            npts=0;
            continue
        else
            %plotm(lat,lon,'--k'); %removed mapping toolbox
            %wrap lons greater than lonmax
            lon(lon>lonmax)=lon(lon>lonmax)-360.;
            m_line(lon,lat,'color','black','linestyle','--');
            tline=fgetl(fid);
            lat=[]; lon=[];
            npts=0;
            continue
        end
    end 
    npts=npts+1;
    lonlat=sscanf(tline,'%f %f');
    lat(npts)=lonlat(2);
    lon(npts)=lonlat(1);
    if(npts>1)
        %segment if longitude wraps around
        if((lon(npts)>lonmax&&lon(npts-1)<=lonmax) ...
                ||(lon(npts)<lonmax-360&&lon(npts-1)>=lonmax-360) ...
                ||(abs(lon(npts)-lon(npts-1))>tolerance*360.))
            lon(lon>lonmax)=lon(lon>lonmax)-360.;
            m_line(lon(1:npts-1),lat(1:npts-1),'color','black','linestyle','--');
            tline=fgetl(fid);
            lat=[]; lon=[];
            npts=1;
            lat(npts)=lonlat(2);
            lon(npts)=lonlat(1);
            continue
        end
    end
    tline=fgetl(fid);
end

