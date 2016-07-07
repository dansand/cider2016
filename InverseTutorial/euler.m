function [ alpha,beta,gamma,del ] = euler( sth,sph,rth,rph )
%EULER Finds Euler angles for a coordinate axis rotation
%   Finds the Euler angles alpha,beta,gamma that rotate the coordinate axes 
%   so that the z-axis is at the pole of the source-receiver great circle 
%   (s x r), and the x-axis is at the source. See Edmonds' Angular Momentum 
%   in Quantum Mechanics, page 7 for the angle conventions.
%     input: sth,sph = source coordinates in radians
%            rth,rph = receiver coordinates in radians
%     output: alpha,beta,gamma = euler angles in radians which rotate the 
%             original coordinate system to the one with the source-receiver 
%             great circle on the equator, the source at (PI/2,0). The minor 
%             arc to the receiver is in the positive phi direction.
%            del = source-receiver separation in radians.
%            pth,pph = source-receiver great circle pole location.

    % Get cartesian coordinates for source and receiver
    scart=tptocart(sth,sph);
    rcart=tptocart(rth,rph);
    del=acos(dot(scart,rcart));
    pcart=cross(scart,rcart);
    pth=atan2(sqrt(pcart(1)^2+pcart(2)^2),pcart(3));
    if(pcart(1)==0 && pcart(2)==0) 
        % special case of pole at z or -z
        pph=0.0;
    else
        pph=atan2(pcart(2),pcart(1));
    end
    alpha=pph;
    beta=pth;
    %  the x'' axis (call it t) is at pth+pi/2,pph
    ttheta=pth + pi/2.0;
    tcart=tptocart(ttheta,pph);
    %  the third Euler angle, gamma, rotates x'' to the source s.
    gamma=acos(dot(scart,tcart));
    %  form q = x'' x s to check the sign of gamma (q/|q| = +-p/|p|)
    qcart=cross(tcart,scart);
    sgn=dot(pcart,qcart);
    if(sgn<0.0) 
        gamma=-gamma;
    end
    return

end

