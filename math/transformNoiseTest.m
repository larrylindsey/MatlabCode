function [e tpts_noisefit tpts_dupfit] = transformNoiseTest(fpts, tr, noise,...
    applyNoise)

% set up noise defaults
if nargin < 4
    applyNoise = 'add';
    if nargin < 3
        noise = .01 * rand(size(fpts));
    end    
end

% If applyNoise is textual, convert it to the appropriate function
if ischar(applyNoise)
    if strcmpi(applyNoise, 'add')
        applyNoise = @(x, y) x + y;
    elseif strcmpi(applyNoise(1:3), 'mul')
        applyNoise = @(x, y) x .* y;
    end
end

tpts = applyTransform(fpts, tr);

% If noise is a function handle, generate some noise.
if ~isnumeric(noise)
    % Assume noise is a function handle if it isn't numeric
    noise = noise(tpts);
end


tpts_noisy = applyNoise(tpts, noise);

tr_dup = refitTransform(fpts, tpts, tr);
tr_noisy = refitTransform(fpts, tpts_noisy, tr);

tpts_dupfit = applyTransform(fpts, tr_dup);
tpts_noisefit = applyTransform(fpts, tr_noisy);

e(1) = trError(tpts, tpts_noisefit);
e(2) = trError(tpts_dupfit, tpts_noisefit);
e(3) = trError(tpts, tpts_dupfit);

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function e = trError(pts1, pts2)
e = rms(sqrt(sum((pts1 - pts2).^2, 2)));
end
