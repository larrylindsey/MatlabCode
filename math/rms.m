function r = rms(in, varargin)

in = in .^ 2;
r = mean(in, varargin{:});
r = sqrt(r);