function [ theta,phi,found ] = getstninfo( stname )
%GETSTNTP Reads file stations and returns coordinates 
    persistent dsta dtype dlat dlon nsta
    
    rconv=180/pi;
    found=0; theta=0.0; phi=0.0;

    % check to see if file has not yet been read
    if isempty(dsta)
        fid=fopen('stations');
        C=textscan(fid,'%5c %2c %f %f %*[^\n]');
        dsta=C{1};
        dtype=C{2};
        dlat=C{3};
        dlon=C{4};
        nsta=size(dsta,1);
        fclose(fid);
    end
    
    for i=1:nsta
        if (dsta(i,1:5)==stname(1:5))
            found=1;
            theta=(90.0-dlat(i))/rconv;
            phi=dlon(i)/rconv;
            return
        end
    end
end

