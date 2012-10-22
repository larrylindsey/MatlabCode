function doSfNPlots(k24, oblique, volumejosef, jstats)


% Plot options
labelsize = 12;
titlesize = 14;

c_k24 = [1 0 0];
c_obl = [0 1 0];
c_vol = [0 0 1];
cmap = cat(2, c_k24', c_obl', c_vol');


opt_k24 = makeVolCompPlots();
opt_k24.axis = 0;
opt_k24.font = 'Arial';

opt_k24.fontsize = labelsize;
opt_k24.color = c_k24;

opt_obl = opt_k24;
opt_obl.color = c_obl;

opt_vol = opt_k24;
opt_vol.color = c_vol;

legnames = {'Spine', 'Oblique', 'Apical'};
xticks = -.5:0.1:.5;
yticks = 0:.1:1;

allax = [];

% Histogram parameters
binctr = linspace(-.5, .5, 50);

% Dendrites
% hold on;
% plot([-.1 -.1 nan .1 .1], [0 1 nan 0 1], 'k', 'LineWidth', 2);

objtype = {'dendrite', 'axon', 'cfa'};
titles = {'Dendrite Volume', 'Axon Volume', 'Synapse Area'};
mtype = {'Volume', 'Volume', 'Flatarea'};

table_label = {'n', '+/- 5%', '+/- 10%', 'mean'};

table_stat = zeros(4, 4, numel(objtype));

for ii = 1:numel(objtype)
    figure;
    firstsub();
    ax = gca;
    opt_k24.axis = ax; opt_obl.axis = ax; opt_vol.axis = ax;
    
    
    stat_k24 = makeVolCompPlots(k24, [objtype{ii} 'Orig'], ...
        [objtype{ii} 'Elastic'], mtype{ii},...
        opt_k24);
    stat_obl = makeVolCompPlots(oblique, [objtype{ii} 'Orig'], ...
        [objtype{ii} 'Elastic'], mtype{ii},...
        opt_obl);
    stat_vol = makeVolCompPlots(volumejosef, [objtype{ii} 'Orig'], ...
        [objtype{ii} 'Elastic'], mtype{ii},...
        opt_vol);
    
    xlim([-.5 .5]);
    
    title(sprintf('Percent Change in %s\nCumulative Distribution', ...
        titles{ii}), 'FontSize', titlesize);
    lh = legend(legnames{:}, 'Location', 'SouthEast');
    set(lh, 'FontSize', 12);
    
    h_k24 = hist(stat_k24, binctr);
    h_obl = hist(stat_obl, binctr);
    h_vol = hist(stat_vol, binctr);
    h_all = cat(1, h_k24, h_obl, h_vol)';
    
    allax = [allax gca];
    
    secondsub();
    bar(binctr, h_all, 1, 'stacked');
    set(gca, 'FontSize', labelsize);
    title('Histogram', 'FontSize', titlesize);
    colormap(cmap);
    xlim([-.5 .5]);
    allax = [allax gca];
    
    table_stat(2:4, 1, ii) = [numel(stat_vol) numel(stat_obl) numel(stat_k24)];
    table_stat(2:4, 2, ii) = [ sum(abs(stat_vol) < .05), sum(abs(stat_obl) < .05),...
            sum(abs(stat_k24) < .05)];
    table_stat(2:4, 3, ii) = [ sum(abs(stat_vol) < .1), sum(abs(stat_obl) < .1),...
            sum(abs(stat_k24) < .1)];
    table_stat(2:4, 4, ii) = [ mean(stat_vol), mean(stat_obl), mean(stat_k24) ];
    
    table_stat(1, 1:3, ii) = sum(table_stat(2:4, 1:3, ii), 1);
    table_stat(1, 4, ii) = mean(cat(2, stat_vol, stat_obl, stat_k24));
    
    
end

table_stat(:, 2, :) = table_stat(:, 2, :) ./ table_stat(:, 1, :);
table_stat(:, 3, :) = table_stat(:, 3, :) ./ table_stat(:, 1, :);

fprintf('Table output, Volumes\n');

for ii = 1:size(table_stat,3)
    fprintf('%s,All, Apical, Oblique, Spine\n', titles{ii});
    for jj = 1:numel(table_label)
        fprintf('%s, %g, %g, %g, %g\n', table_label{jj}, ...
            table_stat(1, jj, ii), table_stat(2, jj, ii), ...
            table_stat(3, jj, ii), table_stat(4, jj, ii));
    end
    fprintf('\n');
end


set(allax, 'XTickMode', 'manual');
set(allax, 'XTick', xticks);
set(allax(1:2:end), 'YTickMode', 'manual');
set(allax(1:2:end), 'YTick', yticks);




end

function firstsub
subplot(11, 1, 2:5);
end

function secondsub
subplot(11, 1, 8:10);
end



