function c = leastSquaresFit(A, z)

ATA = A' * A;

[U S V] = svd(ATA);

invS = S;
invS(logical(eye(size(S)))) = 1./invS(logical(eye(size(S))));

ATAInv = U * invS * V';

eyemat = eye(size(ATAInv));
eyetest = ATAInv * ATA;

invErr = rms(eyemat(:) - eyetest(:));

s = warning('off', 'MATLAB:nearlySingularMatrix');
if invErr > 1e-3
    ATAInvTest = eyemat / ATA;
    eyetest = ATAInv * ATA;
    if rms(eyemat(:) - eyetest(:)) < invErr
        ATAInv = ATAInvTest;
    end
end
warning(s.state, 'MATLAB:nearlySingularMatrix');

c = (ATAInv * A' * z);