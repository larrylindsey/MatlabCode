function imout = compositeImage(im1, im2)

sz1 = size(im1);
sz2 = size(im2);

minsz = min(sz1, sz2);

selr0 = 1:minsz(1);
selc0 = 1:minsz(2);

selr1 = adjustSel(sz1(1), minsz(1), selr0);
selc1 = adjustSel(sz1(2), minsz(2), selc0);
selr2 = adjustSel(sz2(1), minsz(1), selr0);
selc2 = adjustSel(sz2(2), minsz(2), selc0);

imout = zeros([minsz 3]);

imout(:,:,1) = im1(selr1, selc1);
imout(:,:,3) = im2(selr2, selc2);
imout(:,:,2) = mean(imout(:,:,[1 3]), 3);

end

function selOut = adjustSel(sz, minsz, sel)
if sz > minsz
    selOut = sel + ceil((sz - minsz)/2) - 1;
else
    selOut = sel;
end
end



    