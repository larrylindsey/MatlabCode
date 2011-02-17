function pb = applyPBmodelOEOBri( model, img, rFitParab )
% OE + Oriented-smoothed Brightness
% optional parabolic fit

sigma=model.sigma;
elong=model.elong;
if isfield(model,'ntheta')
    ntheta=model.ntheta;
else
    ntheta=numel(model.thetalist);
end

oe = detOE( img, sigma,elong,ntheta);
[oemax,oeangle]=max( oe, [], 3 );

smoothOri = detOE( img, sigma/2,elong*2,ntheta,0); %half width,same length,0=deriv
smooth = zeros(size(oeangle));
for ia=1:ntheta 
    sma = smoothOri(:,:,ia);
    smooth(oeangle==ia) = sma(oeangle==ia);
end
    
if length(model.beta)==length(model.fstd)
    disp('Legacy model: same size fstd and beta'); % as if 1's for intercept explicit in data
    fstd = model.fstd(2:end);
else
    fstd = model.fstd;
end
pb = makePBori( {oe smooth}, ntheta, fstd, model.beta );

if nargin>=3
    pb = fitParabOri( pb, model.thetalist, rFitParab );
end
