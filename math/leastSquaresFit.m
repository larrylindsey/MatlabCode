function c = leastSquaresFit(A, z)

ATA = A' * A;

[U S V] = svd(ATA);

invS = S;
invS(logical(eye(size(S)))) = 1./invS(logical(eye(size(S))));

ATAInv = U * invS * V';

eyemat = eye(size(ATAInv));
eyetest = ATAInv * ATA;

invErr = rms(eyemat(:) - eyetest(:));

if invErr > 1e-3
    ATAInv_alt = inv(ATA);
    eyetest = ATAInv_alt * ATA;
    if rms(eyemat(:) - eyetest(:)) < invErr
        ATAInv = ATAInv_alt;
    end
end

c = (ATAInv * A' * z);