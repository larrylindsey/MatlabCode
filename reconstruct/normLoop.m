function npts = normLoop(pts, c)
%pts in [m 2]
m = size(pts, 1);
npts = zeros(c, 2);

pts = pts - repmat(mean(pts, 1), [m 1]);

mag = mean(sqrt(pts(:,1).^2 + pts(:,2).^2));
pts = pts / mag;

pts = cat(1, pts, pts(1,:))';
%now, pts in [2 m + 1]

theta = linspace(0, 2 * pi, c + 1);
theta = theta(1:(end - 1));

l2(:,1) = 0;

for i_c = 1:c
    l2(1,2) = sin(theta(i_c)) * 100;
    l2(2,2) = cos(theta(i_c)) * 100;    
    
    raymag = [];
    
    for i_p = 1:m
        l1(:,1) = pts(:, i_p);
        l1(:,2) = pts(:, i_p + 1);
        [ipt, doesI] = intersection(l1, l2);
        if doesI
            raymag = [raymag sqrt(sum(ipt.^2))]; %#ok
        end
    end
    raymag = mean(raymag);
    npts(i_c, 1) = raymag * sin(theta(i_c));
    npts(i_c, 2) = raymag * cos(theta(i_c));
end

end