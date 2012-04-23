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
    for i_trans = 1:numel(secdoc(i_sec).section.Transform)
        newtrans(i_trans) = ... 
            processTransform(secdoc(i_sec).section.Transform(i_trans));%#ok
        
        %Record the transform containing image information.
        if ~isempty(secdoc(i_sec).section.Transform(i_trans).Image)
            secdoc(i_sec).section.transImageIndex = i_trans;
            
        end
    end    
    secdoc(i_sec).section.Transform = newtrans;
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function trans = convertTransform(trans)
A = cat(2, trans.xcoef', trans.ycoef');
%A = cat(1, A, zeros(4, 2));


switch trans.dim
    case 2
        A([2 3], 2) = A([3 2], 2);
    case 4
        A([4 5], :) = A([5 4], :);
end

%Tinv = Tinv([7:10 4:6 2:3 1],:);

Tinv = A;

% if ~isempty(trans.Image)
%     k = max(trans.Contour.points(:,1) * trans.Image.mag);
%     %Expect to have coefficients only for order up to one, plus the
%     %coefficient corresponding to xy
% 
%     Tinv(1, 1) = A(1, 1) + k * (A(3, 1) + A(5, 1));
%     Tinv(1, 2) = k - k * (A(3, 2) + A(5, 2)) - A(1, 2);
%     Tinv(3, 1) = -A(3, 1);
%     Tinv(2, 2) = -A(2, 2);
%     Tinv(5, 1) = -A(5, 1);
% end

trans.type = @taylorMat;

%T = invertTransform(Tinv, @taylorMat, 32, [0 10]);

trans.T = [];
trans.Tinv = Tinv;

trans.order = 2;
trans.iDim = 2;
trans.oDim = 2;

trans.data.n = 32;
[trans.data.u trans.data.v] = getTransformSupport(trans);

trans = populateTransInverse(trans);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [u v] = getTransformSupport(trans)

if isempty(trans.Image)
    pts = cat(1, trans.Contour.points);
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

transform = convertTransform(transform);

for i_c = 1:numel(transform.Contour)
    transform.Contour(i_c).transPoints = ...
        applyTransform(transform, transform.Contour(i_c).points);    
end

if ~isempty(transform.Image)
    for i_c = 1:numel(transform.Contour)
        transform.Contour(i_c).imageDomainPoints = ...
            transform.Contour(i_c).points * transform.Image(1).mag;
        transform.Contour(i_c).imageDomainTransPoints = ...
            applyTransform(transform, ...
            transform.Contour(i_c).imageDomainPoints);
    end
end
    

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function ptsout = applyTransform(contour, ptsin)
if ~isempty(ptsin)
    A = taylorMat(ptsin(:,1), ptsin(:,2), contour.order);
    ptsout = A * contour.T;
else
    ptsout = [];
end
end

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
