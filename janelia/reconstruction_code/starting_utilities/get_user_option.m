function o = get_user_option(choices)
% o = get_user_option(choices)
% Get a numerical option from user from a defined set of choices
o = 0;
while(ismember(o, choices)==0)
  fprintf('Option:');
  o = input('');
  if(ismember(o, choices)==0)
    fprintf('Incorrect option - choose again\n');
  end;
end;
return