function [ th,ph ] = iblk( blockindx,bsize,nlat,mlat,hsize )
%IBLK Returns theta, phi of center of block with index blockindx

for i=1:nlat
    j1=mlat(i)+1;
    j2=mlat(i+1);
    if(blockindx>=j1&&blockindx<=j2)
        th=(i-0.5)*bsize;
        ph=(blockindx-j1+0.5)*hsize(i);
        return
    end
end
th=-999; ph=-999;
return
end

