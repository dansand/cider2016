function n = linecount(fid)
%LINECOUNT Counts the number of lines in a file
%  Note: rewinds file before and after calculation
    frewind(fid);
    n=0;
    tline=fgetl(fid);
    while ischar(tline)
        tline=fgetl(fid);
        n=n+1;
    end
    frewind(fid);
    return
end

