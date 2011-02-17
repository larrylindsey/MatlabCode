function [result, segment_map_ws, boundary_map] = ...
      testbed_segmentation_from_superpixels_mean_boundary_value_RAG(case_id)

  result = [];
  switch(case_id)
    case 1
      segment_map_ws = uint32([1 0 2; 1 0 2; 1 0 2]);
      boundary_map = uint8([1 1 1; 1 1 1; 1 1 1]);
      f_threshold_seq = uint8(2);
      boundary_length_seq = uint8(0);
      run_test();
    case 2
      segment_map_ws = uint32([1 0 2; 1 0 2; 1 0 2]);
      boundary_map = uint8([1 2 3; 4 5 6; 7 8 9]);
      f_threshold_seq = uint8(2);
      boundary_length_seq = uint8(0);
      run_test();

    case {3,4}
      test_dir = '../test_data/';
      switch(case_id)
        case 3
          image_prefix = 'lamina.OTO.003';
        case 4
          image_prefix = 'medulla.leginon.681.153';
      end
      image = imread([test_dir, image_prefix, '.tif']);
      boundary_map = 255 - medfilt2(image, [3,3]);
      %  figure(1);
      %imshow(boundary_map);
      %title('Boundary map');
      
      segment_ws_name = [test_dir, image_prefix, '.ws.mat'];
      if(exist(segment_ws_name, 'file')~=2)
        fprintf('Computing initial watershed ...');
        tic
        segment_map_ws = uint32(watershed(boundary_map));
        save2(segment_ws_name, 'segment_map_ws');
        toc
        fprintf('done.\n');
      else
        load2(segment_ws_name, 'segment_map_ws');
      end
      %figure(2);
      %[py, px] = find(segment_map_ws==0);
      %imshow(boundary_map);
      %hold on; plot(px, py, '.'); hold off;
      
      f_threshold_seq = double(190); % double(190:2:200);
      boundary_length_seq = uint32(zeros(1, length(f_threshold_seq)));
      
      fprintf('Computing segmentation ...\n');
      tic
      result = segmentation_from_superpixels_mean_boundary_value_RAG(...
          segment_map_ws, boundary_map, f_threshold_seq, ...
          boundary_length_seq);
      toc
      fprintf('done.\n');
    case 100
      test_dir = '../test_data/';
      image_prefix = 'a.007';
      image = imread([test_dir, image_prefix, '.tif']);
      figure(1);
      imshow(image);
      title('Image');

      boundary_suffix = '.BEL.c_5_mf7';
      boundary = load2([test_dir, image_prefix, boundary_suffix, '.mat']);
      boundary_map = uint8(255*(1-boundary.boundary));
      boundary_map = imerode(imdilate(boundary_map, strel('disk',4)), ...
        strel('disk', 4));
      figure(2);
      imshow(boundary_map);
      title('Boundary map');
      
      segment_ws_name = [test_dir, image_prefix, boundary_suffix, ...
        '.ws.mat'];
      if(exist(segment_ws_name, 'file')~=2)
        fprintf('Computing initial watershed ...');
        tic
        segment_map_ws = uint32(watershed(boundary_map));
        save2(segment_ws_name, 'segment_map_ws');
        toc
        fprintf('done.\n');
      else
        load2(segment_ws_name, 'segment_map_ws');
      end
      figure(3);
      [py, px] = find(segment_map_ws==0);
      imshow(boundary_map);
      hold on; plot(px, py, '.'); hold off;
      
      f_threshold_seq = double(50:1:60);
      boundary_length_seq = uint32(0*ones(1, length(f_threshold_seq)));
      
      fprintf('Computing segmentation ...\n');
      tic
      result = segmentation_from_superpixels_mean_boundary_value_RAG(...
          segment_map_ws, boundary_map, f_threshold_seq, ...
          boundary_length_seq);
      toc
      fprintf('done.\n');
    otherwise
      error('Could not recognize test case id');
  end

  fprintf('Test OK\n');
  return;
  function run_test()
    display(segment_map_ws)
    display(boundary_map)
    display(f_threshold_seq)
    display(boundary_length_seq)
    segmentation_from_superpixels_mean_boundary_value_RAG(...
        segment_map_ws, boundary_map, f_threshold_seq, boundary_length_seq);
    return;
  end
end
