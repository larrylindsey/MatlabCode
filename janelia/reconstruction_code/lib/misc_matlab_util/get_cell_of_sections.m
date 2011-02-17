function c = get_cell_of_sections(start_sections, n_section)
% get_cell_of_sections(start_sections, n_section)
% Returns a cell array, one cell element for each start_section equal to
% [start_section:start_section+n_section-1]

c = {};
for i = 1:length(start_sections)
  c{i} = start_sections(i):start_sections(i)+n_section-1; %#ok<AGROW>
end
return
end
