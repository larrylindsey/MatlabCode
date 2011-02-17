function model = createTriangleModel(s, junk)%#ok
model = cat(1, mean(s(:,[1 2]), 1), mean(s(:,[3 4]), 1));