function [stat, e] = computestat(sname, cname, str, type, minArea)
% sname - stat name
% cname - class name
% str - the struct
% type - 'a' or 'm' for annotation or mito

cal = .002 * .002;

% population
pop = str.(cname).(type) * cal;
pop(pop < minArea * cal) = [];

totarea = str.totpix * cal;

if nargout > 1
    e = stderr(pop);
end

switch sname
    case 'n'
        stat = numel(pop);
    case 'fraction'
        % non capillary pixels
        noncappix = totarea - sum(str.capillary.a * cal);
        stat = sum(pop) / noncappix;
    case 'avgsize'
        stat = mean(pop);
    case 'nnorm'
        % non capillary pixels
        noncappix = totarea - sum(str.capillary.a * cal);
        stat = numel(pop) / noncappix;
    otherwise
        error('don''t know how to compute stat %s', sname);
end

end