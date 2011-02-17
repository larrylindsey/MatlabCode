function is_hash_correct = check_hash_compatibility()
% Check to make hash function is compatibile
% Important: Don't change the following values in
% code/lib/hash_functions/hash_multiplication_method.m
% STRING_INDEX_CYCLE = 4;
% A = (sqrt(5)-1)/2;
% HASH_RANGE = 2^20;

global config_global

hash_config = config_global.storage.hash;

hash_check_0 = strcmp(hash_string('Ticket to Ride', hash_config), '716847682745_40')==1;
hash_check_1 = strcmp(hash_string('Eleanor Rigby', hash_config), '170623585413_40')==1;
hash_check_2 = strcmp(hash_string('Norwegian Wood', hash_config), '1024830632835_40')==1;

is_hash_correct = hash_check_0 && hash_check_1 && hash_check_2;

return
end
