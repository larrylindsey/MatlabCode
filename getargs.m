function [argstr argcell] = getargs(argname, argcell)

istomp = [];

for ii = 1:numel(argname)
    argstr.(argname{ii}) = [];
end

for ii = 1:numel(argcell)
    key = argcell{ii};
    if ischar(key) &&  ii < numel(argcell)
        for jj = 1:numel(argname)
            if strcmpi(key, argname{jj})
                argstr.(argname{jj}) = argcell{ii + 1};
                istomp = [istomp ii ii + 1]; %#ok<AGROW>
            end
        end
    end
end

argcell(istomp) = [];


end