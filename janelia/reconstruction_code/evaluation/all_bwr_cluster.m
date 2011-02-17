function all_bwr_cluster(output_fn, gt_fn, ws_fns, amb_fns, ts, index_sets, ...
                         colname_base, trim, min_body_size, ...
                         markers_per_body, submatrix_multiplier, d)
    if ~isdeployed,
        fprintf(['this function is not designed to be run from Matlab, ' ...
                 'but from the linux command line after compilation.']);
        return;
    end
    spec = struct();
    spec.gt = readVTK(gt_fn);
    fprintf('Read gt correctly.\n');
    spec.ws_filenames = eval(ws_fns);
    fprintf('Read ws filenames correctly.\n');
    spec.amb_file_patterns = eval(amb_fns);
    fprintf('Read amb filenames correctly.\n');
    spec.ts = eval(ts);
    fprintf('Read ts correctly.\n');
    spec.index_sets = eval(index_sets);
    fprintf('Read index_sets correctly.\n');
    spec.colname_base = eval(colname_base);
    fprintf('Read colnames correctly.\n');
    spec.trim = eval(trim);
    fprintf('Read trim correctly.\n');
    spec.min_body_size = eval(min_body_size);
    fprintf('Read min_body_size correctly.\n');
    spec.markers_per_body = eval(markers_per_body);
    fprintf('Read markers_per_body correctly.\n');
    spec.submatrix_multiplier = eval(submatrix_multiplier);
    fprintf('Read submatrix_multiplier correctly.\n');
    spec.d = eval(d);
    [prtable, colnames] = all_bwr(spec);
    dlmwrite(output_fn, prtable, 'delimiter', ',', 'precision', 12);
    header = fopen([output_fn '.header.csv'], 'w');
    for i=1:(numel(colnames)-1),
        fprintf(header, [colnames{i} ',']);
    end
    fprintf(header, [colnames{end} '\n']);
end
