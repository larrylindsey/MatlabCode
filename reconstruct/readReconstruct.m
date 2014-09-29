function [serdoc secdoc] = readReconstruct(ser, sec)


%Read XML files into a matlab struct
if nargin < 2
    %expect ser to be a string.
    dot = find(ser == '.');
    prefix = ser(1:(dot - 1));
    
    serdoc = xmlToMatstruct(readXML(ser));
    
    secdoc = readSections(prefix);
elseif nargin == 2
    %expect ser, sec to be structs
    serdoc = ser;
    secdoc = sec;
end

%Convert Reconstruct data.
for i_sec = 1:numel(secdoc)
    fprintf('Converting section %d\n', secdoc(i_sec).index);
    clear newtrans;
    if isfield(secdoc(i_sec).section, 'Transform')
        for i_trans = 1:numel(secdoc(i_sec).section.Transform)
            newtrans(i_trans) = ...
                processTransform(secdoc(i_sec).section.Transform(i_trans));%#ok
            
            %Record the transform containing image information.
            if ~isempty(newtrans(i_trans).Image)
                secdoc(i_sec).section.transImageIndex = i_trans;
                
            end
        end
    else
        secdoc(i_sec).section.transImageIndex = 1;
        newtrans = repmat(struct('dim', 0,...
            'xcoef', 0, ...
            'ycoef', 0,...
            'Image', 0,...
            'Contour', 0,...
            'type', 0,...
            'T', 0,...
            'doTrans', 0,...
            'createTrans', 0,...
            'matrixFun', 0,...
            'param', 0,...
            'order', 0,...
            'monomialOrder', 0,...
            'ndim', 0,...
            'data', 0,...
            'inv', 0), 0);        
    end
    secdoc(i_sec).section.Transform = newtrans;
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function trans = convertTransform(trans)
A = cat(2, trans.xcoef', trans.ycoef');
%A = cat(1, A, zeros(4, 2));

ptype = trans.dim;
if ptype == 2
    A([2 3], 2) = A([3 2], 2);
end
% switch ptype
%     case 2
%         A([2 3], 2) = A([3 2], 2);
%     case 4
%         A([4 5], :) = A([5 4], :);
% end

T = A;
trans.type = @taylorMat;

%monomial order
morder = ...
    [ 0 0;
    1 0;
    0 1;
    1 1;
    2 0;
    0 2];

if ptype > 1
    order = 2;
else
    order = 1;
    morder = morder(1:3,:);
    T = T(1:3,:);
end
[T morder] = rectifyOrder(T, morder, order);


trans.T = T;
[trans.doTrans trans.createTrans] = taylorMat();
trans.matrixFun = @taylorMat;
trans.param = trans.createTrans();
trans.order = order;
trans.monomialOrder = morder;
trans.ndim = 2;

trans.data.n = 32;
[trans.data.u trans.data.v] = getTransformSupport(trans);

trans = fitInverseTransform(trans);
trans.inv.data = trans.data;
trans = invertTransStruct(trans);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [u v] = getTransformSupport(trans)

if ~isfield(trans, 'Image') || isempty(trans.Image)
    if isempty(trans.Contour)
        u = [];
        v = [];
        return;
    else
        pts = [];
        for ii = 1:numel(trans.Contour)
            if ~isempty(trans.Contour(ii).points)
                pts = cat(1, pts, trans.Contour(ii).points);
            end
        end        
    end
else
    cname = {trans.Contour.name};
    %     idom = regexp(cname, '^domain');
    %     while isempty(idom{1})
    %         idom = {idom{2:end}};
    %     end
    
    %Assume the this trans has only one Contour, and it corresponds to the
    %image.
    pts = trans.Contour(1).points * trans.Image.mag;
end

u = [min(pts(:,1)) max(pts(:,1))];
v = [min(pts(:,2)) max(pts(:,2))];

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function transform = processTransform(transform)

if ~isfield(transform, 'Image')
    transform.Image = '';
    keyboard
end

transform = convertTransform(transform);

for i_c = 1:numel(transform.Contour)
    transform.Contour(i_c).transPoints = ...
        applyTransform(transform.Contour(i_c).points, transform);
end

if isfield(transform, 'Image') && ~isempty(transform.Image)
    for i_c = 1:numel(transform.Contour)
        transform.Contour(i_c).imageDomainPoints = ...
            transform.Contour(i_c).points * transform.Image(1).mag;
        transform.Contour(i_c).imageDomainTransPoints = ...
            applyTransform(transform.Contour(i_c).imageDomainPoints, transform);
    end
end


end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% function ptsout = applyReconstructTransform(contour, ptsin)
% if ~isempty(ptsin)
%     A = taylorMat(ptsin(:,1), ptsin(:,2), contour.order);
%     ptsout = A * contour.T;
% else
%     ptsout = [];
% end
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function secStr = readSections(prefix)

dd = dir([prefix '.*']);
slash = find(prefix == '/', 1, 'last');
if isempty(slash)
    serPath = './';
else
    serPath = prefix(1:slash);
end

secStr = repmat(struct('index', 0, 'name', '', 'section', []), [1 numel(dd)]);

keepsel = true(1, numel(dd));

parfor i_f = 1:numel(dd)
    keepsel(i_f) = ~isempty(regexp(dd(i_f).name, '\.[0-9]+$', 'once'));
    if keepsel(i_f)
        name = dd(i_f).name;
        dot = find(name == '.', 1, 'last');
        secStr(i_f).index = str2double(name((dot+1):end));
        secStr(i_f).name = [serPath dd(i_f).name];
    end
end

secStr = secStr(keepsel);

secIndex = [secStr.index];
[junk sortIndex] = sort(secIndex); %#ok

secStr = secStr(sortIndex);

parfor i_sec = 1:numel(secStr)
    fprintf('Reading section %d\n', secStr(i_sec).index);
    secStr(i_sec).section = xmlToMatstruct(readXML(secStr(i_sec).name));
end

end
