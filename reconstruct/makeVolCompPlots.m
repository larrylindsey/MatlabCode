function varargout = makeVolCompPlots(str, f1, f2, stat, opts)

if nargin < 5
    opts.color = 'blue';
    opts.fontsize = 14;
    opts.axis = [];
    opts.linewidth = 2;
    opts.font = 'Helvetica';
end

if nargin < 1
    varargout = {opts};
    return;
else
    varargout = {};
end

if isempty(str.(f1))
    varargout = {[]};
    return;
end

popCon = [str.(f1).(stat)];
popExp = [str.(f2).(stat)];

plotArgs = {'Color', opts.color, 'LineWidth', opts.linewidth};
axisArgs = {'FontSize', opts.fontsize, 'FontName',...
    opts.font};

% Cumulative distribution of percentage divergence
if isempty(opts.axis)
    figure(1);
    ax = gca;
else
    ax = opts.axis;
end
    
hold on;
plot(ax, sort((popExp - popCon) ./ popCon), linspace(0,1,numel(popCon)), ...
    '*-', plotArgs{:});
xlim([-1 1] * max(abs(xlim)));
grid on;
set(gca, axisArgs{:});
% 
% % Percentage divergence vs. size
% figure(2);
% plot(popCon, (popExp - popCon) ./ popCon, '.', plotArgs{:});
% hold on;
% grid on;
% set(gca, axisArgs{:});
% 
% % Direct comparison
% figure(3);
% plot(popCon, popExp, '.', plotArgs{:});
% hold on;
% grid on;
% set(gca, axisArgs{:});

if nargout > 0
    varargout = {(popExp - popCon) ./ popCon };
end

end
