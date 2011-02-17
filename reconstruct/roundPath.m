function [out sel] = roundPath(inpath)

if min(size(inpath)) ~= 1
    error('input path must be a vector');
end

rpath = round(inpath);
drpath = inpath - rpath;

sel = false(size(drpath));%logical(drpath == 0);

dsig = drpath(2:end) ./ (drpath(1:(end - 1)) + eps);

sel(find(dsig < 0) + 1) = true;

out = inpath(sel);