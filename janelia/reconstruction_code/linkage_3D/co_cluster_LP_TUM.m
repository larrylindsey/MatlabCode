function P = co_cluster_LP_TUM(F, A, file_name_prefix)
% F(F==0) = -1;

fprintf('Dumping linear program.\n');
lp_file_name = [file_name_prefix, '.lp'];
A_TUM = get_ineq_TUM_2LDH_MIAS_v4(double(A));
dump_linear_program_p3_v2(F, double(A_TUM),lp_file_name);

fprintf('Optimizing.\n');
lp_sol_file_name = [file_name_prefix, '.sol'];
lp_solv_file_name = [file_name_prefix, '.solv'];
tic
system(['lp_solve -lp ', lp_file_name, ' > ', lp_sol_file_name]);
system(['cat ', lp_sol_file_name, ' | grep "x_" | ', ...
  'awk ''{print $0;}'' > ', lp_solv_file_name]);
toc

P = eye(size(F));
fin = fopen(lp_solv_file_name, 'rt');
while(true)
  v_name = fscanf(fin, '%s', [1 1]);
  if(feof(fin)==1)
    break;
  end
  v_value = 1-fscanf(fin, '%g', 1);
 
  v_name_brkt = strfind(v_name, '_');
  node_i = str2double(v_name(v_name_brkt(1)+1:v_name_brkt(2)-1)) + 1;
  node_j = str2double(v_name(v_name_brkt(2)+1:end)) + 1;
  P(node_i, node_j) = v_value;
  P(node_j, node_i) = v_value;
end
fclose(fin);

delete([file_name_prefix, '.lp']);
delete([file_name_prefix, '.sol']);
delete([file_name_prefix, '.solv']);
return;
end
