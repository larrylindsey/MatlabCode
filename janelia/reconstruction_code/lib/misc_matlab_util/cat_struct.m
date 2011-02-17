function u = cat_struct(s, t)
% concatenate the fields of structures s and t into u
u = cell2struct(cat(1,struct2cell(s),struct2cell(t)), ...
  cat(1,fieldnames(s),fieldnames(t)), 1);
return
end
