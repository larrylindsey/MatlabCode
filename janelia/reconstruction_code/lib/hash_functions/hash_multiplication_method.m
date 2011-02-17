function h = hash_multiplication_method(s, hash_config)
% Generate hash value from a character string using Multiplication method,
% see Introduction to Algorithms, Cormen et al., Page 233.
%
% Inputs:
% 1. s      variable length string
% Outputs:
% 1. h      integer
%
% Steps:
% 1. A given string is partitioned into substrings of length
% STRING_INDEX_CYCLE = hash_config.string_index_cycle.
% 2. A radix-notation integer is computed for each substring.
% 3. The radix-notations values of all substrings are added together.
% 4. A hash-value is computed for the sum using multiplication method,
% where the real-valued factor is A = hash_config.multiplication_factor,
% and the range for the hash-value is defined as HASH_RANGE = ...
% hash_config.hash_range. 
%
% Shiv N. Vitaladevuni
% Janelia Farm Research Campus, HHMI
%
% v0  12152008    init code
%

STRING_INDEX_CYCLE = hash_config.string_index_cycle;
A = hash_config.multiplication_factor;
HASH_RANGE = hash_config.hash_range;

p = double(1);
for i = 1:length(s)
  p = p + double(int16(s(i)))*(2^mod(i,STRING_INDEX_CYCLE));
end

% multiplication method for hashing
h = floor(HASH_RANGE*mod(p*A,1));

return
end
