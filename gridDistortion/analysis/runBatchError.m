function runBatchError(dirs, param)
wd = pwd;
for ii = 1:numel(dirs)
    cd(dirs{ii});
    load distortion.mat
    load exclude.mat
    trblock = makeTrBlock(rc_found, rc_aff, param);
    trconstrblock = constrainBlock(trblock(1,3:5));
    [eraw, eblock] = distortionTransformErrorBlock(trblock, 'siftCache.mat', exclude);
    [erawc, eblockc] = distortionTransformErrorBlock(trconstrblock, 'siftCache.mat', exclude);
    save error_analysis eraw eblock erawc eblockc trblock trconstrblock
    cd(wd)
end
