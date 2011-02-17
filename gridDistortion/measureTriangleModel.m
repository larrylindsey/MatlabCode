function measure = measureTriangleModel(m, s, extra)

lmodel = repmat(m(1,:), [size(s, 1), 1]);
rmodel = repmat(m(2,:), [size(s, 1), 1]);

cvect = extra(s(:,5), :);
nbd_ind = s(:,6:end);
nbd_ind = unique(nbd_ind(nbd_ind > 0));
nbd_vect = extra(nbd_ind(:), :);

lexp = cvect + lmodel;
rexp = cvect + rmodel;
uexp = cvect - lmodel;
dexp = cvect - rmodel;

dmapL = min(dist2(lexp, nbd_vect), [], 2);
dmapR = min(dist2(rexp, nbd_vect), [], 2);
dmapU = min(dist2(uexp, nbd_vect), [], 2);
dmapD = min(dist2(dexp, nbd_vect), [], 2);

dmap = min(cat(2, dmapL, dmapR, dmapU, dmapD), [], 2);

measure = sqrt(max(dmap)) / mean(rms(m, 2) * sqrt(2));