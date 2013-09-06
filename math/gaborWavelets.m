function GW = gaborWavelets(parm)

retparm = nargin > 0 && ischar(parm);
defaultparm = nargin < 1 || retparm;

if defaultparm
    parm.R = 128;
    parm.C = 128;
    parm.Kmax = pi / 2;
    parm.f = sqrt( 2 );
    parm.Delt = 2 * pi;
    parm.vs = 0:4;
    parm.us = 1:8;   
end

if retparm
    GW = parm;
    return;
end

R = parm.R;
C = parm.C;
Kmax = parm.Kmax;
f = parm.f;
Delt = parm.Delt;
us = parm.us;
vs = parm.vs;
Delt2 = Delt * Delt;


GW = zeros(R, C, numel(vs), numel(us));


for iv = 1:numel(vs)
    for iu = 1:numel(us)
        v = vs(iv);
        u = us(iu);
        GW(:,:, iv, iu) = GaborWavelet ( R, C, Kmax, f, u, v, Delt2 ); % Create the Gabor wavelets
    end
end


end

function GW = GaborWavelet (R, C, Kmax, f, u, v, Delt2)
% Create the Gabor Wavelet Filter
% Author : Chai Zhi  
% e-mail : zh_chai@yahoo.cn

k = ( Kmax / ( f ^ v ) ) * exp( 1i * u * pi / 8 );% Wave Vector

kn2 = ( abs( k ) ) ^ 2;

GW = zeros ( R , C );

for m = -R/2 + 1 : R/2
    
    for n = -C/2 + 1 : C/2
        
        GW(m+R/2,n+C/2) = ...
            ( kn2 / Delt2 ) * exp( -0.5 * kn2 * ( m ^ 2 + n ^ 2 ) / Delt2) * ...
            ( exp( 1i * ( real( k ) * m + imag ( k ) * n ) ) - ...
            exp ( -0.5 * Delt2 ) );
    
    end

end
end
