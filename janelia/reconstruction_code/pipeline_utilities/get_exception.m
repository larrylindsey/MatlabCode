function value = get_exception(exceptions, plane_id, parameter_name)
value = [];
for i = 1:size(exceptions, 1)
  if(ismember(plane_id, exceptions{i,1}))
    j = strfind(exceptions{i,2}, ['-', parameter_name]);
    if(isempty(j))
      continue;
    end
    j = j + length(parameter_name) + 1;
    sub_str = exceptions{i,2}(j:end);
    j = strfind(sub_str, ' ');
    if(~isempty(j))
      sub_str = sub_str(1:j-1);
    end
    eval(['value.', parameter_name, sub_str]);
  end
end
return
end
