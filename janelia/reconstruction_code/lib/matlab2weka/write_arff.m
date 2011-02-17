function write_arff(dataset_name, features, labels, file_name)

if(nargin<4)
  file_name = [dataset_name, '.arff'];
end

fout = fopen(file_name, 'wt');
fprintf(fout, '@relation %s\n', dataset_name);

for i = 1:size(features,2)
  fprintf(fout, '@attribute feature_%d real\n', i);
end
fprintf(fout, '@attribute label {p, n}\n');

fprintf(fout, '@data\n');
for i = 1:size(features,1)
  fprintf(fout, '%g, ', features(i, :));
  if(labels(i)==1)
    fprintf(fout, 'p\n');
  else
    fprintf(fout, 'n\n');
  end
end

fclose(fout);

return
end
