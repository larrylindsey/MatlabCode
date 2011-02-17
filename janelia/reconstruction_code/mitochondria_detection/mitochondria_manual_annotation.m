function mitochondria_manual_annotation(config)
% mitochondria_manual_annotation(config)
% for manually annotating mitochondria using GUI
%
% Shiv N. Vitaladevuni
% Janelia Farm Research Campus, HHMI
%
% v0  ~Feb.2008  init code
% v1  04112008  modified for reconstruction pipeline
%

fprintf('Manual annotation of mitochondria for training ..\n');

is_debugging = config.DEBUG;

mitochondria_config = config.mitochondria;
train_config = mitochondria_config.train;

image_dir = [get_reconstruction_dir(config), mitochondria_config.dir, train_config.dir];


message_1 = 'Delineate a mitochondria by clicking on vertices of polygon.';
message_2 = 'Press p to cancel previous vertex, d to stop tracing the mitochondria.';

for case_id  = train_config.case_ids
  fprintf('%d ', case_id);
  if(mod(case_id, 20)==0)
    fprintf('\n');
  end;
  
  img = im2double( imread( sprintf([image_dir, train_config.image_prefix, ...
    train_config.image_suffix], case_id) ) );
  img = img/max(img(:));
  img = histeq(img);
  figure(1);
  set(gcf, 'Menu', 'none');
  set(gcf, 'Toolbar', 'none');

  flag_last_zone = 0;

  zone = [];
  n_zone = 0;
  while(flag_last_zone==0)
    figure(1);
    imshow(img);
    set(gcf, 'name', message_1);
    title(message_2);
    hold on;
    for z = 1:length(zone)
      plot(zone{z}(1,:), zone{z}(2,:), 'go-');
    end;
    n=0;
    but=1;
    xy=[];
    set(gcf, 'currentcharacter', '0');
    hPoint = [];
    while(but==1)
      press_flag = waitforbuttonpress();
      if(press_flag==0)
        c_p = get(gca, 'currentpoint');
        xi = c_p(1,1); yi = c_p(1,2);
        n=n+1;
        xy(:,n) = [xi; yi];
        figure(1);
        imshow(img);
        set(gcf, 'name', message_1);
        title(message_2);
        hold on;
        for z = 1:length(zone)
          plot(zone{z}(1,:), zone{z}(2,:), 'go-');
        end;
        plot(xy(1,:), xy(2,:), 'ro-');
        hold off;
      else
        i = get(gcf, 'currentcharacter');
        switch(i)
          case 'd'
            but=0;
          case 'p'
            if(n>0)
              n=n-1;
              xy = xy(:, 1:end-1);
              if(~isempty(hPoint))
                imshow(img);
                set(gcf, 'name', message_1);
                title(message_2);
              end;
              hold on;
              for z = 1:length(zone)
                plot(zone{z}(1,:), zone{z}(2,:), 'go-');
              end;
              hPoint = plot(xy(1,:), xy(2,:), 'ro-');
              hold off;
            end;
        end;
        set(gcf, 'currentcharacter', '0');
      end;
    end;

    set(gcf, 'name', 'Save this mitochondria?');
    title('Save & continue (y)? Cancel & continue (c)? Save & go to next image (n)? Cancel & go to next image (a)?');
    press_flag = waitforbuttonpress();
    while(press_flag~=1)
      press_flag = waitforbuttonpress();
    end;
    i = get(gcf, 'currentcharacter');
    switch(i)
      case 'y'
        n_zone = n_zone + 1;
        zone{n_zone} = xy;
      case 'c'
      case 'n'
        n_zone = n_zone + 1;
        zone{n_zone} = xy;
        break;
      case 'a'
        break;
    end;

  end;

  zone_mask = zeros(size(img));
  for i = 1:n_zone
    zone_mask = max(zone_mask, poly2mask(zone{i}(1,:), zone{i}(2,:), size(zone_mask,1), size(zone_mask,2)) );
  end;
  figure(2); imshow(zone_mask);
  title('Mask of the mitochondria');

  zone_mask = 255*(zone_mask>0);

  save2( sprintf([image_dir, train_config.image_prefix, train_config.annot_suffix, '.mat'], ...
    case_id), 'zone_mask', 'zone');

end;

end

