function e = stderr(v)
e = std(v(:)) / sqrt(numel(v));
end