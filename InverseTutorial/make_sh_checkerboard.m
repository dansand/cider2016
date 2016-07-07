function [ ckmodel ] = make_sh_checkerboard(l,nblk,bsize,nlat,mlat,hsize)
%MAKE_SH_CHECKERBOARD Makes a spherical harmonic checkerboard of order l
%   m is set to 1/2 l to give an approximate spherical harmonic derived
%   checkerboard
convert=pi/180.;

m=round(0.5*l);

ckmodel=zeros(nblk,1);
for i=1:nblk
    [th,ph]=iblk(i,bsize,nlat,mlat,hsize);
    c=cos(th*convert);
    s=sin(th*convert);
    x=shfcn(l,c,s);
    argm=m*ph*convert;
    %frp=x(m+1)*cos(argm);
    fip=x(m+1)*sin(argm);
    if(fip>0)
        ckmodel(i)=10;
    else
        ckmodel(i)=-10;
    end
end
return

end

