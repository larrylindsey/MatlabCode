function current_dir = pwd2()

d = pwd();
current_dir = strrep(d, '\', '/');

return;
end

