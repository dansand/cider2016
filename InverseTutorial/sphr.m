function [ sphcoord ] = sphr( cartcoord )
%SPHR compute spherical coordinates (r,t,p) 
%  from cartesian coordinates (x,y,z). t is colatitude and 
%  p is longitude in radians. p is returned in the interval (0,2*pi).
   
    sphcoord(1)=norm(cartcoord);
    if (sphcoord(1)<=0)
        sphcoord(2)=0.;
        sphcoord(3)=0.;
        return
    end
    t=sqrt(cartcoord(1)^2+cartcoord(2)^2);
    sphcoord(2)=atan2(t,cartcoord(3));
    if (cartcoord(1)~=0)
        sphcoord(3)=atan2(cartcoord(2),cartcoord(1));
        if (sphcoord(3)<0)
            sphcoord(3)=sphcoord(3)+2.0*pi;
        end
        return
    end
    if (cartcoord(2)<0)
        sphcoord(3)=1.5*pi;
    elseif (cartcoord(2)==0)
        sphcoord(3)=0.;
    else
        sphcoord(3)=0.5*pi;
    end
    return
end

