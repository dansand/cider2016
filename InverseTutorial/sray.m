function [ row,delt ] = sray( sth,sph,rth,rph,norb,nblk,nlat,bsize,mlat,hsize )
%SRAY Determine a row of G matrix according to surface wave ray
%   theory

    row=zeros(1,nblk);
    rad=pi/180.0;
    pi2=0.5*pi;
    
    % Determine Euler angles to rotate source and receiver to equator with
    % source at x axis
    if(norb<10)
        [alpha,beta,gamma,delta]=euler(sth,sph,rth,rph);
    else
        beta=sth;
        alpha=sph;
        gamma=0.;
    end
    calf=cos(alpha);
    salf=sin(alpha);
    cbet=cos(beta);
    sbet=sin(beta);
    cgam=cos(gamma);
    sgam=sin(gamma);
    if(norb>10)
        delta=2.0*pi;
    else
        if(mod(norb,2)==1)
            delta=delta+(norb-1)*pi;
        else
            delta=delta-norb*pi;
        end
    end
    delt=abs(delta)/rad;
%*** note that delta is negative for even norb
    dp=0.01*rad;
    if(delta<0.d0) 
        dp=-dp;
    end
    delp=abs(dp);
    %ps=0.;
%*** move around equator with small increments and add incremental path
%lengths to appropriate block index for row in G matrix
%Vectorize for matlab speed
    ps=dp:dp:delta;
    x=cos(ps);
    y=sin(ps);
    t1=cgam*x-sgam*y;
    t2=sgam*x+cgam*y;
    u1=cbet*t1;
    z=-sbet*t1;
    x=calf*u1-salf*t2;
    y=salf*u1+calf*t2;
    [az,el]=cart2sph(x,y,z);
    az(az<0)=az(az<0)+2.0*pi;
    t0=(pi2-el)/rad;
    p0=az/rad;
    ib=fblk(t0,p0,nlat,bsize,mlat,hsize);
    for i=1:size(ps,2)
        row(ib(i))=row(ib(i))+delp/rad;
    end
    %Do last step to delta if necessary
    if(abs(ps(end))<abs(delta))
        pslast=delta;
        delp=abs(pslast-ps(end));
        xlast=cos(pslast);
        ylast=sin(pslast);
        t1last=cgam*xlast-sgam*ylast;
        t2last=sgam*xlast+cgam*ylast;
        u1last=cbet*t1last;
        zlast=-sbet*t1last; 
        xlast=calf*u1last-salf*t2last;
        ylast=salf*u1last+calf*t2last;
        [azlast,ellast]=cart2sph(xlast,ylast,zlast);
        if(azlast<0)
            azlast=azlast+2.0*pi;
        end 
        t0last=(pi2-ellast)/rad;
        p0last=azlast/rad;
        iblast=fblk(t0last,p0last,nlat,bsize,mlat,hsize);
        row(iblast)=row(iblast)+delp/rad;
    end
    %while(ps~=delta)
    %    ps0=ps;
    %    ps=ps+dp;
    %    if(dp<0.&&ps<delta) 
    %        ps=delta;
    %    end
    %    if(dp>0.&&ps>delta) 
    %        ps=delta;
    %    end
    %    delp=abs(ps-ps0);
%*** convert to cartesian
    %    x=cos(ps);
    %    y=sin(ps);
%*** use euler angles to rotate back to geographic coordinates
    %    t1=cgam*x-sgam*y;
    %    t2=sgam*x+cgam*y;
    %    u1=cbet*t1;
    %    z=-sbet*t1;
    %    x=calf*u1-salf*t2;
    %    y=salf*u1+calf*t2;
    %    [az,el]=cart2sph(x,y,z);
    %    t0=(pi2-el)/rad;
    %    p0=az/rad;
    %    ib=fblk(t0,p0,nlat,bsize,mlat,hsize);
    %    row(ib)=row(ib)+delp/rad;      
    %end

end

