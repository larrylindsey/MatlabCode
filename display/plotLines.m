function plotLines(lines, ah)

if nargin < 2
    ah = gca;
end

hold(ah, 'on');

for k = 1:length(lines)
    xy = [lines(k).point1; lines(k).point2];
    
    plot(ah, xy(:,1), xy(:,2), 'LineWidth', 2, 'Color', 'green');
    plot(ah, xy(1,1), xy(1,2), 'x', 'LineWidth', 2, 'Color', 'yellow');
    plot(ah, xy(2,1), xy(2,2), 'x', 'LineWidth', 2, 'Color', 'red');
end