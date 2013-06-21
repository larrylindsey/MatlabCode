function vfix = graph_cc(v, vis)
% vfix = graph_cc(v)
%  Find connected components, fixup the vertices so that connected components
%  all map to the smallest vertex index.
%
% v should be in [n 2]

if nargin < 2
    vis = false;
end

if isempty(v)
    vfix = v;
    return;
end

v = validate(v);

% vertex stack
vs = unique(v(:,1));

vfix = [];
if vis
    hh = waitbar(0, 'Processing Edges');
end
os = numel(vs);

while ~isempty(vs)
    % current vertex index
    [icv vs] = pop(vs);
    
    % next vertices stack    
    [vfix,~,vs] = bfs(vfix, icv, icv, vs, v);
    
    if vis
        waitbar(1 - (numel(vs) / os), hh);
    end
end

vfix = sortEdges(vfix);

if vis
    close(hh);
end

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [vfix inew vs] = bfs(vfix, icv, inew, vs, v)

if isempty(vfix)
    prev_sel = false;
else
    prev_sel = vfix(:,2) == icv;
end

if any(prev_sel)
   it = vfix(prev_sel,1);
   if it < inew
       cat(1, vfix, [it inew]);
       rep_sel = vfix(:,1) == inew;
       vfix(rep_sel,1) = it;
       inew = it;
   else
       cat(1, vfix, [inew it]);
       rep_sel = vfix(:,1) == it;
       vfix(rep_sel,1) = inew;
   end
end


vs(vs == icv) = [];
if inew ~= icv && notAdded(icv, vfix)
    vfix = cat(1, vfix, [inew icv]);
end
sel = v(:,1) == icv;
nvs = v(sel,2);

for ii = 1:numel(nvs)
    [vfix inew vs] = bfs(vfix, nvs(ii), inew, vs, v);
end

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function v = validate(v)
% If v is of size [n 2], ensure that:
% * v(i, 1) <= v(i, 2) for 1 <= i <= n, and
% * v(i, 1) <= v(j, 1) for 1 <= i <= j <= n

sel = v(:,1) > v(:,2);

vtemp = v(sel,:);
v(sel,1) = vtemp(sel,2);
v(sel,2) = vtemp(sel,1);

v = sortEdges(v);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function v = sortEdges(v)
n = max(v(:));
vm = v(:,1) + v(:,2) / n;
[~,isort] =  sort(vm);
v = v(isort,:);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ii, v] = pop(v)
ii = v(1);
v = v(2:end);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function yn = notAdded(icv, vfix)
if isempty(vfix)
    yn = true;
else
    yn = all(vfix(:,2) ~= icv);
end
end

