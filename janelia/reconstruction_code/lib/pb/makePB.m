%
% OE + Oriented Brightness
%
modelBest=model(3,1);
for islice=1:2
    pbRaw{islice} = applyPBmodelOEOBri( modelBest, al{islice}, rFitParab );
    pb{islice} = max( pbRaw{islice}, [], 3 );
    pbsoft{islice} = max( gentleNonmax(pbRaw{islice}), [], 3 );
end
save('pb-OEOBri-HRP1-ShivAnnot_1_2.mat','code','pb','pbsoft');

