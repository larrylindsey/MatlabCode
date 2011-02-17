function testbed_get_connected_components_c()
message = 'test case 1: disconnected graph';
fprintf([message, '\n']);
A = eye(10);
threshold = 0;
expected_answer = 1:10;
process_test_case();

message = 'test case 2: connected graph but high threshold';
fprintf([message, '\n']);
A = ones(10);
threshold = 2;
expected_answer = 1:10;
process_test_case();

message = 'test case 3: fully connected graph';
fprintf([message, '\n']);
A = ones(10);
threshold = 0;
expected_answer = ones(1, 10);
process_test_case();

message = 'test case 4: cliques';
fprintf([message, '\n']);
A = zeros(10);
A(1:3, 1:3) = 1;
A(4:7, 4:7) = 1;
A(8:10, 8:10) = 1;
threshold = 0;
expected_answer = [1 1 1 2 2 2 2 3 3 3];
process_test_case();

message = 'test case 5: chains';
fprintf([message, '\n']);
A = zeros(10);
A(1,2)=1; A(2,3)=1;
A(4,5)=1; A(5,6)=1; A(6,7)=1;
A(8,9)=1; A(9,10)=1;
threshold = 0;
expected_answer = [1 1 1 2 2 2 2 3 3 3];
process_test_case();

message = 'test case 6: stars';
fprintf([message, '\n']);
A = zeros(10);
A(1,2)=1; A(2,3)=1;
A(4,6)=1; A(5,6)=1; A(6,7)=1;
A(8,6)=1; A(9,6)=1;
A(10,10)=1;
threshold = 0;
expected_answer = [1 1 1 2 2 2 2 2 2 3];
process_test_case();

return

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Sub routines
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  function process_test_case()
    A = max(A, A');
    fprintf('A:\n');
    disp(A);
    fprintf('threshold:\n');
    disp(threshold);
    fprintf('Expected answer:\n');
    disp(expected_answer);
    A = sparse(A);
    if(~issparse(A))
      error('Could not make A sparse');
    end
    c = double(get_connected_components_c(sparse(A), threshold)');
    fprintf('Obtained answer:\n');
    disp(c);
    if(sum(abs(c-expected_answer))~=0)
      error(['Error in ', message]);
    end
    return
  end
end
