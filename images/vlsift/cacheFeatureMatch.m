function cacheFeatureMatch(siftCache, rparm, force)

if nargin < 3
    force = false;    
    if nargin < 2
        rparm = ransacRegressionTransform;
    end
end

cache = load(siftCache);

if ~and(isfield(cache, 'f'), isfield(cache, 'd'))
    error('%s is not a sift cache file', siftCache);
end

if or(force, ~isfield(cache, 'm'))
    m = cell(1, numel(cache.f) - 1);
    s = m;
    d = cache.d;
    d_pp = d(2:end);
    parfor ii = 1:numel(m)
        [m{ii} s{ii}] = thirdpartyfunction('vl_ubcmatch', d{ii}, d_pp{ii});
    end
    cache.m = m;
    cache.s = s;
    save(siftCache, '-v7.3', '-struct', 'cache');
end

mr = cell(size(cache.m));
sr = mr;
feat1 = mr;
feat2 = mr;

f = cache.f;
f_pp = f(2:end);
m = cache.m;
s = cache.s;
sz = cache.sz;

parfor ii = 1:numel(mr)
    f1 = f{ii};
    f2 = f_pp{ii};
    mm = m{ii};
    ss = s{ii};
    f1 = f1(1:2, mm(1,:))';
    f2 = f2(1:2, mm(2,:))';
    f1 = f1 * 2 / max(sz(ii,:)) - 1;
    f2 = f2 * 2 / max(sz(ii,:)) - 1;
    [~, matchSel] = ransacRegressionTransform(rparm, f1, f2, 1);
    mr{ii} = mm(:, matchSel);
    sr{ii} = ss(:, matchSel);
    feat1{ii} = f1(matchSel,:);
    feat2{ii} = f2(matchSel,:);
end

cache.rparm = rparm;
cache.mr = mr;
cache.sr = sr;
cache.feat1 = feat1;
cache.feat2 = feat2;

save(siftCache, '-v7.3', '-struct', 'cache');
