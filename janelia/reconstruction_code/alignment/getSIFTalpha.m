function adjacencies=getSIFTalpha(tforms,sizes)
% adjacencies=getSIFTalpha(tforms,sizes)
% Construct overlaps matrix from tforms for set of tiles

%Get bounding rectangles for tiles
X={}; Y={};
Xi={}; Yi={};
for i=1:size(tforms,1)
  for j=1:size(tforms,2)
    if(isempty(tforms{i,j}))
      continue;
    end

    x=[sizes{i,j}(1)*[0 0 1 1];sizes{i,j}(2)*[0 1 1 0]];
    x=tforms{i,j}(1:2,1:2)*x+repmat(tforms{i,j}(:,3),[1 size(x,2)]);

    X{i,j}=x(1,:);
    Y{i,j}=x(2,:);
    
    xb = 1:10:sizes{i,j}(1);
    yb = 1:10:sizes{i,j}(2);
    xm = repmat(xb, [length(yb) 1]);
    ym = repmat(yb', [1 length(xb)]);
    x = [xm(:)'; ym(:)'];
    x=tforms{i,j}(1:2,1:2)*x+repmat(tforms{i,j}(:,3),[1 size(x,2)]);

    Xi{i,j}=x(1,:);
    Yi{i,j}=x(2,:);
  end
end


%construct overlaps matrix
adjacencies=zeros(size(tforms,1),size(tforms,2),size(tforms,1),size(tforms,2));
for i=1:size(tforms,1)
  for j=1:size(tforms,2)
    if(isempty(tforms{i,j}))
      continue;
    end

    %check only adjacent slices
    for i1=max(1,i-1):i
      for j1=1:size(tforms,2)
        if(isempty(tforms{i1,j1}))
          continue;
        end
        %check overlap conditions - if any of the vertices
        %from one polygon is inside the other -- there is overlap
        IN = inpolygon(Xi{i,j},Yi{i,j},X{i1,j1},Y{i1,j1});
        IN1 = inpolygon(Xi{i1,j1},Yi{i1,j1},X{i,j},Y{i,j});

        if(max([IN,IN1])>0)
          adjacencies(i,j,i1,j1)=1;
          adjacencies(i1,j1,i,j)=1;
        end
      end
    end
  end
end
