function sub_strings = strsplit2(s, split_str, remove_esc_char)
% sub_strings = strsplit2(s, split_str)
% Split string with escape character \

  if(nargin==1)
    sub_strings = strsplit(s);
    return;
  end

  split_loc = regexp(s, ['[^\\]', split_str]);
  if(~isempty(split_loc))
    split_loc = split_loc + 1;
  end

  if(isempty(split_loc))
    sub_strings = {s};
  else
    split_loc = [0, split_loc, length(s)+1];
    sub_strings = {};
    for i = 1:length(split_loc)-1
      sub_strings{i} = s(split_loc(i)+1:split_loc(i+1)-1);
    end
  end

  if(nargin>=3 && remove_esc_char)
    for i = 1:length(sub_strings)
      sub_strings{i} = regexprep(sub_strings{i}, '(?<char>[^\\])\', '$<char>');
    end
  end
  return;
end
