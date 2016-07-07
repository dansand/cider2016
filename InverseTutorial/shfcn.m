function [ x,dxdth ] = shfcn( l,c,s )
%SHFCN COmpute ordinary spherical harmonics
%  computes ordinary spherical harmonics (Edmonds eqn 2.5.29)
%    x(l,m,theta)*exp(imphi) are a set of orthonormal spherical 
%    harmonics where theta is colatitude and phi is longitude.
%  input:
%    l=harmonic degree
%    c=cos(theta)
%    s=sin(theta) 
%  output:
%    x(1) contains m=0, x(2) contains m=1, ... x(l+1) contains m=l
%    where m=azimuthal order 0.le.m.le.l .  
%    dxdth (theta derivatives) stored in the same way as x
%  calls no other routines
%       TGM
%  translated quickly to matlab from Guy Master's fortran
      lp1=l+1;
      fl2p1=l+lp1;
      con=sqrt(fl2p1/(4*pi));
      x=zeros(lp1,1);
      dxdth=zeros(lp1,1);
      
%*** handle special case of l=0
      if(l==0) 
        x(1)=con;
        dxdth(1)=0.d0;
        return
      end
      if(abs(s)<1.e-20)
%*** handle very small arguments
        x(1)=con;
        dxdth(2)=-0.5d0*con*sqrt(l*lp1);
        return
      end
%*** try m recursions ; first compute xll
      f=1.0;
      for i=1:l
        f=f*(i+i-1.d0)/(i+i);
      end
      if(log(s)*l>-679d0)
%*** use m recurrence starting from m=l
        cot=c/s;                      
        x(lp1)=con*sqrt(f)*(-s)^l;
        dxdxth(lp1)=x(lp1)*l*cot;
        for i=1:l
          m=lp1-i;
          mp1=m+1;   
          f=sqrt(i*(fl2p1-i));
          x(m)=-(dxdxth(mp1)+m*cot*x(mp1))/f;
          dxdxth(m)=(m-1)*x(m)*cot+x(mp1)*f;
        end
      else
%*** use the Libbrecht algorithm
        c2=c+c;
        x(lp1)=con*sqrt(f)*(-1.d0)^l;
        dxdxth(lp1)=0.d0;
        for i=1:l
          m=lp1-i;
          mp1=m+1;    
          f=sqrt(i*(fl2p1-i));
          x(m)=-(s*dxdxth(mp1)+m*c2*x(mp1))/f;
          dxdxth(m)=s*x(mp1)*f;
        end
%*** now convert back to ordinary spherical harmonics
        fac=1.d0;
        for i=2:lp1
          dxdxth(i)=(s*dxdxth(i)+(i-1)*c*x(i))*fac;
          x(i)=x(i)*s*fac;
          fac=fac*s;
        end
      end
      return

end

