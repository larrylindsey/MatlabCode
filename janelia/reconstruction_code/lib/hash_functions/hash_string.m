function h = hash_string(s, hash_config)
% h = hash_string(s)
% Compute a hash value "h" for a variable length string "s"
%
% Inputs:
% 1. s      variable length string
% Outputs:
% 1. h      variable length hash value: string
%
% Shiv N. Vitaladevuni
% Janelia Farm Research Campus, HHMI
%
% v0  12152008    init code
%

h1 = num2str( hash_multiplication_method(s(1:2:end), hash_config) );
h2 = num2str( hash_multiplication_method(s(2:2:end), hash_config) );

hash_version = [num2str(hash_config.string_index_cycle), ...
  num2str(hash_config.multiplication_factor), num2str(hash_config.hash_range)];
hash_config.hash_range = 2^6;
hash_version_hash = ...
  ['_', num2str(hash_multiplication_method(hash_version, hash_config))];

h = [h1, h2, hash_version_hash];

return
end
