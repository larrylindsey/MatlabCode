function square_images
setappdata(0,'UseNativeSystemDialogs',0)
[files path] = uigetfile('*.*', 'Select Images', 'MultiSelect', 'on');

if ischar(files)
    files = {files};
end

for ii = 1:numel(files)
    file = files{ii};
    try
        disp(['Reading ', file]);
        imdata = imread([path file]);
        goflag = 1;
    catch
        le = lasterror;
        disp(['Error reading ' file]);
        disp(['Error message ' le.message]);
        disp(['Error ID ' le.identifier]);
        fprintf('Skipping...\n\n');
        goflag = 0;
    end
    
    if goflag == 1
        imsz = size(imdata);
        bigdim = max(imsz(1:2));
        if numel(imsz) < 3
            imsz(3) = 1;
        end
        newimdata = zeros(bigdim, bigdim, imsz(3));

        if isinteger(imdata)
            imclass = class(imdata);
            maxval = intmax(imclass);
            newimdata = cast(newimdata, imclass);
        else
            maxval = 1.0;
        end

        newimdata = newimdata + maxval;

        dim1lim = [1 bigdim];
        dim2lim = [1 bigdim];

        if imsz(1) > imsz(2)
            szdiff = imsz(1) - imsz(2);
            fstep = floor(szdiff / 2);
            bstep = ceil(szdiff / 2);
            dim2lim(1) = 1 + fstep;
            dim2lim(2) = dim2lim(2) - bstep;
        elseif imsz(2) > imsz(1)
            szdiff = imsz(2) - imsz(1);
            fstep = floor(szdiff / 2);
            bstep = ceil(szdiff / 2);
            dim1lim(1) = 1 + fstep;
            dim1lim(2) = dim1lim(2) - bstep;
        end
        dim1sel = dim1lim(1):dim1lim(2);
        dim2sel = dim2lim(1):dim2lim(2);
        newimdata(dim1sel,dim2sel,:) = imdata;

        %Figure out the new filename
        finddot = findstr(file, '.');
        if numel(finddot) > 0
            lastdot = finddot(end);
            newfile = [file(1:lastdot) 'sqr' file(lastdot:end)];
        else
            newfile = [file 'sqr.png'];
        end
        imwrite(newimdata, [path newfile]);
        disp(['Wrote ' newfile]);
        fprintf('\n\n');
    end
end