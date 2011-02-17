function mov =  volume2avi(Vol, filename, xyz, frame_axis, params, matched_volumes)

if nargin < 3,
    xyz{3} = 1:size(Vol,3);
    xyz{1} = 1:size(Vol,1);
    xyz{2} = 1:size(Vol,2);
end
if nargin < 4,
    frame_axis = 3;
end
image_axes = setdiff([1 2 3], frame_axis);
if nargin < 5,
    params.fps = 3;
end

if ~isfield(params, 'fps'),
    params.fps=3;
end

if nargin < 6,
    matched_volumes = {};
end
if ~isfield(params, 'paste_dim'), 
    params.paste_dim = 2;
end

original_mov_size = [numel(xyz{image_axes(1)}) numel(xyz{image_axes(2)}) 3 ...
                                                    numel(xyz{frame_axis})];
mov_size = original_mov_size;
mov_size(params.paste_dim) = original_mov_size(params.paste_dim) * ...
                        (numel(matched_volumes)+1) + numel(matched_volumes);

mov = zeros(mov_size);

Vol = Vol(xyz{1},xyz{2},xyz{3});
for i=1:numel(matched_volumes),
    matched_volumes{i} = matched_volumes{i}(xyz{1},xyz{2},xyz{3});
end

for i=1:numel(xyz{frame_axis}),
    switch(frame_axis)
        case 1
            xind = i;
            yind = 1:size(Vol, 2);
            zind = 1:size(Vol, 3);
        case 2
            xind = 1:size(Vol, 1);
            yind = i;
            zind = 1:size(Vol, 3);
        case 3
            xind = 1:size(Vol, 1);
            yind = 1:size(Vol, 2);
            zind = i;
    end
    for j=1:3,
        if ndims(Vol) == 4, % the image is color
            mov(1:original_mov_size(1),1:original_mov_size(2),j,i) = ...
                                                        Vol(xind,yind,zind,j);
        else
            mov(1:original_mov_size(1),1:original_mov_size(2),j,i) = ...
                                                        Vol(xind,yind,zind);
        end
        for k=1:numel(matched_volumes),
            switch params.paste_dim
                case 1
                    mov_xind = (k*original_mov_size(1)+1+k): ...
                                            ((k+1)*original_mov_size(1)+k);
                    mov_yind = 1:original_mov_size(2);
                case 2
                    mov_xind = 1:original_mov_size(1);
                    mov_yind = (k*original_mov_size(2)+1+k): ...
                                            ((k+1)*original_mov_size(2)+k);
            end
            if ndims(matched_volumes{k}) == 4, % the image is color
                mov(mov_xind, mov_yind, j, i) = ...
                                matched_volumes{k}(xind, yind, zind, j);
            else
                mov(mov_xind, mov_yind, j, i) = ...
                                matched_volumes{k}(xind, yind, zind);
            end
        end
    end
end
mov = immovie(double(mov)/255);
movie2avi(mov, filename, 'compression', 'None', 'fps', params.fps);

end