function doSfNPlots(k24, oblique, volumejosef)

figure;
firstsub();
ax = gca;


% Plot options
labelsize = 12;
titlesize = 14;

c_k24 = [1 0 0];
c_obl = [0 1 0];
c_vol = [0 0 1];
cmap = cat(2, c_k24', c_obl', c_vol');


opt_k24 = makeVolCompPlots();
opt_k24.axis = ax;

opt_k24.fontsize = labelsize;
opt_k24.color = c_k24;

opt_obl = opt_k24;
opt_obl.color = c_obl;

opt_vol = opt_k24;
opt_vol.color = c_vol;

legnames = {'Spine', 'Oblique', 'Apical'};
xticks = -.5:0.05:.5;
yticks = 0:.1:1;

allax = [];

% Histogram parameters
binctr = linspace(-.5, .5, 50);

% Dendrites
% hold on;
% plot([-.1 -.1 nan .1 .1], [0 1 nan 0 1], 'k', 'LineWidth', 2);

stat_k24 = makeVolCompPlots(k24, 'dendriteOrig', 'dendriteElastic', 'Volume',...
    opt_k24);
stat_obl = makeVolCompPlots(oblique, 'dendriteAffine', 'dendriteElastic',...
    'Volume', opt_obl);
stat_vol = makeVolCompPlots(volumejosef, 'dendriteAffine', 'dendriteElastic',...
    'Volume', opt_vol);
xlim([-.5 .5]);

title(sprintf('Percent Change in Dendrite Volume\nCumulative Distribution'),...
    'FontSize', titlesize);
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

% Axons
figure;
firstsub();
ax = gca;

% hold on;
% plot([-.1 -.1 nan .1 .1], [0 1 nan 0 1], 'k', 'LineWidth', 2);

opt_k24.axis = ax; opt_obl.axis = ax; opt_vol.axis = ax;

stat_k24 = makeVolCompPlots(k24, 'axonOrig', 'axonElastic', 'Volume',...
    opt_k24);
stat_obl = makeVolCompPlots(oblique, 'axonAffine', 'axonElastic',...
    'Volume', opt_obl);
stat_vol = makeVolCompPlots(volumejosef, 'axonAffine', 'axonElastic',...
    'Volume', opt_vol);
xlim([-.5 .5]);
title(sprintf('Percent Change in Axon Volume\nCumulative Distribution'),...
    'FontSize', titlesize);

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

% CFAs
figure;
firstsub();
ax = gca;

% hold on;
% plot([-.1 -.1 nan .1 .1], [0 1 nan 0 1], 'k', 'LineWidth', 2);

opt_k24.axis = ax; opt_obl.axis = ax; opt_vol.axis = ax;

stat_k24 = makeVolCompPlots(k24, 'cfaOrig', 'cfaElastic', 'Flatarea',...
    opt_k24);
stat_obl = makeVolCompPlots(oblique, 'cfaAffine', 'cfaElastic',...
    'Flatarea', opt_obl);
stat_vol = makeVolCompPlots(volumejosef, 'cfaAffine', 'cfaElastic',...
    'Flatarea', opt_vol);
xlim([-.5 .5]);
title(sprintf('Percent Change in Contact Flat Area\nCumulative Distribution'),...
    'FontSize', titlesize);

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


set(allax, 'XTickMode', 'manual');
set(allax, 'XTick', xticks);
set(allax(1:2:end), 'YTickMode', 'manual');
set(allax(1:2:end), 'YTick', yticks);

end

function firstsub
subplot(11, 1, [2:5]);
end

function secondsub
subplot(11, 1, 8:10);
end



