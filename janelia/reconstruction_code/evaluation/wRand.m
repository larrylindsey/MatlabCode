function [oversegm,undersegm] = wRand( test, gt )
% weighted Rand coefficient
% elements weighted by inverse size of gt segment they belong to
% treats all labels equally, including zero
n = length(test(:));
if n~=length(gt(:)), error('Partitions should be of the same size'); end

[gt_sort,gt_idx] = sort(gt(:));
test_by_gt = test(gt_idx);
test_by_gt = test_by_gt(:); %should be column
% weights: count gt class sizes
steps = [find(diff(gt_sort)); length(gt_sort)];
classSizes = [steps(1); diff(steps)];
weights = 1./classSizes;

% overlaps & oversegm
oversegm = 0;
ovlpSizes=[];
ovlpWeights=[];
ovlpTestLabels=[];
from=1;
for k=1:length(classSizes)
    to = steps(k);
    testChunk = sort( test_by_gt(from:to) );
    from = to+1;

    chunkSteps = [find(diff(testChunk)); length(testChunk)];
    chunkOvlpSizes = [chunkSteps(1); diff(chunkSteps)];
    chunkOvlpWeights = chunkOvlpSizes * weights(k);
    chunkOvlpTestLabels = testChunk( chunkSteps );
    
    ovlpSizes = [ovlpSizes; chunkOvlpSizes];
    ovlpWeights = [ovlpWeights; chunkOvlpWeights];
    ovlpTestLabels = [ovlpTestLabels; chunkOvlpTestLabels];
    
    oversegm = oversegm + (1 - weights(k)*weights(k)*sum( chunkOvlpSizes.^2 ) )/2;

end

% re-sort to compute undersegm
[ovlpTestLabels,ovlpIdx] = sort( ovlpTestLabels );
ovlpWeights = ovlpWeights(ovlpIdx);
steps = [find(diff(ovlpTestLabels)); length(ovlpTestLabels)];
undersegm = 0;
from=1;
for k=1:length(steps)
    to = steps(k);
    chunk = ovlpWeights(from:to);
    from = to+1;

    undersegm = undersegm + (sum(chunk)^2 - sum(chunk.^2))/2;
end

% normalize by total weighted sum of pairs
norm_oversegm = sum((classSizes-1)./classSizes)/2;
oversegm = oversegm / norm_oversegm;
norm_undersegm = length(classSizes)*(length(classSizes)-1)/2;
undersegm = undersegm / norm_undersegm;

end

%%
function test_wRand()
%%
a=[1 2];
b=[1 1];
A=[a a a];      %  1 2 1 2 1 2
B=[b 2*b 2*b];  %  1 1 2 2 2 2
%{
6 elem, 15 pairs
1)gt=B, cl.sizes 2,4; w=1/2,1/4; 
 oversegm: 1 pairs in gt==1, 4 pairs in gt==2; (1*(1/4)+4*(1/16)) = 1/2
  norm_oversegm=1/4+3/8=5/8;  result (1/2)/(5/8)=4/5
 undersegm: 2 pairs in each test class; 2*2*((1/2)*(1/4)) = 1/2
  norm_undersegm=1
2)gt=A, cl.sizes 3,3; w=const=1/3; 
 oversegm: 2 pairs in each class; 2*2*(1/9)=4/9
  norm_oversegm=2*2/(2*3)=2/3;  result (4/9)/(2/3)=2/3
 undersegm: 1 pair in test==1, 4 pairs in test==2; 5*(1/9)=5/9
  norm_undersegm=1
%}
perm=randperm(length(A));
[oversegm1,undersegm1] = wRand( A(perm), B(perm) );
[oversegm2,undersegm2] = wRand( B(perm), A(perm) );
disp('Test wRand...');good=true;
if abs(oversegm1-4/5)>eps,     disp(['oversegm1=' num2str(oversegm1)]);    good=false;end
if abs(undersegm1-1/2)>eps,    disp(['undersegm1=' num2str(undersegm1)]);  good=false;end
if abs(oversegm2-2/3)>eps,     disp(['oversegm2=' num2str(oversegm2)]);    good=false;end
if abs(undersegm2-5/9)>eps,     disp(['undersegm2=' num2str(undersegm2)]);  good=false;end
if good, disp('OK'); end
%%
a=[1 2 3];
b=[1 1 1];
A=[a a a a a a];
B=[b 2*b 2*b 3*b 3*b 3*b];
%%
perm=randperm(length(A));
[oversegm,undersegm] = wRand( A(perm), B(perm) )
[oversegm,undersegm] = wRand( B(perm), A(perm) )

%%
end
