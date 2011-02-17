function coordsout = try_iter(im, template, PP, coordsin, step)
%crap.  make it work in 2d.
%close(1);
%i_f = figure(1);
%image(template);
%i_a = gca;
%coordsin in dim x coord, ie, 2x6

xydim = size(coordsin, 1);
ncoord = size(coordsin, 2);
dirstep = [-step step];


%Distort the image as according to the input distortion coordinates
[im_d mask] = distort_image(im, PP, coordsin);

%Calculate the correlation given the current distortion coordinates
corr0 = correlation(im_d.*mask, template.*mask);
coord_results = repmat(struct('coords', zeros(2, ncoord), 'corr', 0),...
    [xydim + 2, 1]);

coord_results(1).coords = coordsin;
coord_results(1).corr = corr0;

%Treat x, y dimensions as independent.
for i_dim = 1:xydim
    %Test results will be stored in a struct for simplicity
    %This struct is indexed by [direction (-/+), dimension]
    coord_results_dim = repmat(struct('coords', zeros(2, ncoord), 'corr', 0),...
        [2 ncoord]);

    for i_step = 1:2
        for i_c = 1:ncoord
            trycoords = coordsin;
            trycoords(i_dim, i_c) = trycoords(i_dim, i_c) + dirstep(i_step);
            coord_results_dim(i_step, i_c).coords = trycoords;
            [im_d mask] = distort_image(im, PP, trycoords);
            coord_results_dim(i_step, i_c).corr = correlation(im_d.*mask, ...
                template.*mask);
        end
    end

    corr_results = [coord_results_dim.corr];
    corr_results = corr_results(:);
    
    bestcorr = max(corr_results(:));
    ibest = find(corr_results(:) == bestcorr, 1, 'first');
    coord_results(i_dim + 1) = coord_results_dim(ibest);
end

trycoords = zeros(size(coordsin));
for i_dim = 1:xydim
    trycoords(i_dim,:) = coord_results(i_dim).coords(i_dim,:);
end

coord_results(end).coords = trycoords;
[im_d mask] = distort_image(im, PP, trycoords);
coord_results(end).corr = correlation(im_d.*mask, template.*mask);

corr_results = [coord_results.corr];
bestcorr = max(corr_results(:));
ibest = find(corr_results == bestcorr, 1, 'first');

coordsout = coord_results(ibest).coords;