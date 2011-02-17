function adjmatr=linkPropagate( ovlpList, P )
% build adjacency matrix by propagating link propensities P

nseg=max([ [ovlpList.loIdx] [ovlpList.upIdx] ]);

adjmatr=zeros(nseg,nseg);
for k=1:length(ovlpList)
    ovlp=ovlpList(k);
    adjmatr( ovlp.loIdx , ovlp.upIdx  ) = P(k);
    adjmatr( ovlp.upIdx , ovlp.loIdx  ) = P(k); %symm!
end

% Transitive closure: Floyd-Warshall, but instead of min/+ use max/min
for k=1:nseg
    prev=adjmatr;
    for i1=1:nseg, 
        for i2=(i1+1):nseg
            adjmatr(i1,i2)=max( prev(i1,i2), min(prev(k,i1),prev(k,i2)) );
            adjmatr(i2,i1) = adjmatr(i1,i2);
        end
    end
    %if 0==mod(k,1000),  disp(['Propagate k=' num2str(k) ' time ' num2str(toc)]); end
end
