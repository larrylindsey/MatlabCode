function recPathToXML(path, name, color, fid)
colorstr = sprintf('%03.f %03.f %03.f', color(1), color(2), color(3));

if nargin > 3
    

    
    fprintf(fid, ['<ZContour name="%s" closed="false" border="%s"' ...
        ' fill="%s" mode="11"\n points="'], name, colorstr, colorstr);
    for i_p = 1:size(path,1)
        fprintf(fid, '%f %f %d,\n\t', path(i_p,1), path(i_p,2), path(i_p,3));
    end
    fprintf(fid, '"/>\n');

else

    fprintf(['<ZContour name="%s" closed="false" border="%s"' ...
        ' fill="%s" mode="11"\n points="'], name, colorstr, colorstr);
    for i_p = 1:size(path,1)
        fprintf('%f %f %d,\n\t', path(i_p,1), path(i_p,2), path(i_p,3));
    end
    fprintf('"/>\n');
end



end