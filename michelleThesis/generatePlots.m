function generatePlots(stat, kind)
% generatePlots(stat, kind)
% Generate plots on a stat struct as returned by launchMitoAnalysis
%
% stat - the stat struct
% kind - one of 'g', 'ag', 'hg', or 'hc'
%   g: histogram of all stats by group
%  ag: plots of stats by animal and annotation class (39 in total)
%  hg: histogram of area size by group and annotation class
%  hc: histogram of area size by annotation class only

if ~isfield(stat, 'group')
    error('Need group info. Did you run applyGroupType?');
end

% if ~isfield(stat, 'calibrated')
%     error('Did you run calibrateStat?');
% end

global name_k strfields_f name_title_k strfield_title_f name_id_map ...
    strfields_id_map f k doMito_k g minArea name_g label_font_size ...
    title_font_size hist_n axis_font_size subplot_reduction ...
    pbox_aspect ebar_line_width mbar_line_width marker_line_width ...
    marker_size mean_bar_width marker_color

% class names in the stat struct
name_k = {'terminal', 'glia', 'capillary', 'background', 'allMito', ...
    'terminal_s', 'terminal_m', 'terminal_l', 'glia_s', 'glia_m', 'glia_l'};
% stat names, referenced to computestat.m
strfields_f = {'nnorm', 'fraction', 'avgsize'};
% class title as it will appear in the plot title
name_title_k = {'Terminal', 'Glia', 'Capillary', 'Background', ...
    'Mitochondria',...
    'Small Terminals', 'Medium Terminals', 'Large Terminals', ...
    'Small Glia', 'Medium Glia', 'Large Glia'};
% stat name as it will appear in the title and y-axis
strfield_title_f = {'# / Area (um^{-2})', 'Fraction', 'Size (um^2)'};
% group names, as they appear on titles and axes
name_g = {'YV', 'YE', 'YEP', 'AV', 'AE', 'AEP'};

% number of histogram bins
hist_n = 32;

% The group histograms are six to a figure. Reduce the font size globally
% by this many points
subplot_reduction = 2;

% Marker size for plot kind ag
marker_size = 18;
% Aspect ratio for the ag plots
pbox_aspect = [4 5 1];
% error bar line weight
ebar_line_width = 2;
% mean bar line weight
mbar_line_width = 2;
% marker line weight (plot ag only)
marker_line_width = 3;
% The width (or length, really) of the mean bar
mean_bar_width = .5;
% Marker colors for ag plots
marker_color = [.25 .25 .25];
% Minimum area cutoff for annotations and mitochondria, in that order
minArea = [2500 256];

% Title font size
title_font_size = 20;
% ylabel, xlabel font size
label_font_size = 16;
% Font size for other axes text
axis_font_size = 14;

k = numel(name_k);
f = numel(strfields_f);

name_id_map = struct;
strfields_id_map = struct;

doMito_k = 1:12 < 3; % eventually, we'll have k = 12.

for i_k = 1:k
    name_id_map.(name_k{i_k}) = i_k;
end

for i_f = 1:f
    strfields_id_map.(strfields_f{i_f}) = i_f;
end

% cal_f = [1 1 .002 * .002];

g = numel(name_g);
groups = [stat.group];
g_idx = unique(groups);

stat_g = cell(1, g);

for i_g = 1:g
    gsel = groups == g_idx(i_g);
    stat_g{i_g} = stat(gsel);
end

switch kind
    case 'g'
        plotDataByGroup(stat_g);
    case 'ag'
        plotDataByAnimalAndGroup(stat_g);
    case 'hg'
        plotHistogramsByGroup(stat_g);
    case 'hc'
        plotAreaHistogramsByClass(stat_g);
    otherwise
        error('I don''t know how to make a plot of kind %s', kind);
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function plotAreaHistogramsByClass(stat_g)
global name_k

for i_k = 1:5
    cname = name_k{i_k};
    kclass = [stat.(cname)];
    summary_stat.(cname).a = [kclass.a];
    figure;
    histogramHelper(summary_stat, 'a', 'All', i_k, [0 0]);
    print(gcf, '-dpng', sprintf('all-histogram-%s.png', cname));
    print(gcf, '-depsc', sprintf('all-histogram-%s.eps', cname));
end

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function plotHistogramsByGroup(stat_g)
global k doMito_k g
cstat_g = collateByGroup(stat_g);

for i_k = 1:k
    groupHistogramPlotHelper(cstat_g, i_k, 'a');
end

for i_k = find(doMito_k)
    groupHistogramPlotHelper(cstat_g, i_k, 'm');
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function groupHistogramPlotHelper(cstat_g, i_k, type)
global name_k name_g hist_n g strfields_f minArea
i_f = 3;

figure;
cname = name_k{i_k};

stat_all = [cstat_g.(cname)];
data_all = [stat_all.(type)];

if type == 'a'
    data_all(data_all < minArea(1)) = [];
elseif type == 'm'
    data_all(data_all < minArea(2)) = [];
end

data_all = data_all * .002 * .002;

[~, hist_x] = hist(log10(data_all), hist_n);
dx = hist_x(2) - hist_x(1);
xmin = min(hist_x) - dx / 2;
xmax = max(hist_x) + dx / 2;
minXTick = ceil(2 * xmin) / 2;
maxXTick = floor(2 * xmax) / 2;
altXTick = minXTick:.5:maxXTick;

ymax = 0;
axes = zeros(1, g);

sc = linspace(.5, 0, g);

for i_g = 1:g    
    subplot(2, 3, i_g);
    hh = histogramHelper(cstat_g(i_g), type, name_g{i_g}, i_k, minArea, hist_x);
    set(hh, 'FaceColor', sc(i_g) * [1 1 1]);
    axes(i_g) = gca;
    g_ymax = max(ylim);
    if g_ymax > ymax
        ymax = g_ymax;
    end
end

set(gcf,'units','normalized','outerposition',[0 0 1 1])

drawnow;

for i_g = 1:g
    set(axes(i_g), 'YLim', [0 ymax], 'XLim', [xmin xmax]);
    xTick = get(axes(i_g), 'XTick');
    if numel(xTick) < numel(altXTick)
        set(axes(i_g), 'XTick', altXTick);  
    end
end

drawnow;

print(gcf, '-dpng', sprintf('histogram_%s_%s_%s.png',...
    name_k{i_k}, strfields_f{i_f}, type));
print(gcf, '-depsc', sprintf('histogram_%s_%s_%s.eps',...
    name_k{i_k}, strfields_f{i_f}, type));

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function hh = histogramHelper(stat, type, name, i_k, minArea, hist_x)
i_f = 3;
global strfield_title_f label_font_size title_font_size name_title_k ...
    hist_n name_k axis_font_size subplot_reduction

cname = name_k{i_k};
data = stat.(cname).(type);

if type == 'a'
    data(data < minArea(1)) = [];
elseif type == 'm'
    data(data < minArea(2)) = [];
end

data = data * .002 * .002;

if nargin > 5
    hist_arg = hist_x;
else
    hist_arg = hist_n;
end

hold on;
all_ch = get(gca, 'Children');

hist(log10(data), hist_arg);

hh = setdiff(get(gca, 'Children'), all_ch);

if type == 'a'
    title(sprintf('%s - %s', name, name_title_k{i_k}), ...
        'FontSize', title_font_size - subplot_reduction, 'FontWeight', 'bold');   
elseif type == 'm'
    title(sprintf('%s - %s\nMitochondria', name, name_title_k{i_k}), ...
        'FontSize', title_font_size - subplot_reduction, 'FontWeight', 'bold');
end

%     xTicks = get(gca, 'XTick');
%     xTickLabels = 10.^xTicks;
%     set(gca, 'XTickLabel', xTickLabels);
xlabel(['Log_{10} ' strfield_title_f{i_f}], 'FontSize', label_font_size - subplot_reduction);

set(gca, 'FontSize', axis_font_size - subplot_reduction);
set(gca, 'Box', 'off');
set(gca, 'TickLength', [0 0]);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function cstat_g = collateByGroup(stat_g)
global g

cstat_g = cell(g, 1);

for i_g = 1:6
    cstat_g{i_g} = computeSummary(stat_g{i_g});
end

cstat_g = cat(1, cstat_g{:});

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function plotDataByGroup(stat_g)
global k f doMito_k

cstat_g = collateByGroup(stat_g);

for i_k = 1:k
    for i_f = 1:f
        generateGroupBarPlotHelper(cstat_g, i_k, i_f, 'a');
    end
end

for i_k = find(doMito_k)
    for i_f = 1:f
        generateGroupBarPlotHelper(cstat_g, i_k, i_f, 'm');
    end
end

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function generateGroupBarPlotHelper(cstat_g, i_k, i_f, type)
global strfields_f minArea g name_k ebar_line_width
figure;
b = zeros(1, g);
e = b;

cname = name_k{i_k};

for i_g = 1:g
    [b(i_g), e(i_g)] = computestat(strfields_f{i_f}, cname,...
        cstat_g(i_g), type, minArea);
end
bar(b);
hold on;
if i_f == 3 % std err only means something for mean area
    for i_g = 1:g
        plot([i_g, i_g], [b(i_g) - e(i_g), b(i_g) + e(i_g)], ...
            'k', 'LineWidth', ebar_line_width);
    end
end

makeItPretty(type, i_k, i_f);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function plotDataByAnimalAndGroup(stat_g)
global k f doMito_k

for i_k = 1:k
    for i_f = 1:f
        dataByAnimalPlotHelper(stat_g, i_k, i_f, 'a');
    end
end

for i_k = find(doMito_k)
    for i_f = 1:f
        dataByAnimalPlotHelper(stat_g, i_k, i_f, 'm');
    end
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function dataByAnimalPlotHelper(stat_g, i_k, i_f, type)
global g strfields_f name_k minArea marker_line_width pbox_aspect ...
    mbar_line_width ebar_line_width mean_bar_width marker_color

% this will be referenced later as a char-array
style = 'oxvds^*';

hf = figure;
hold on;
cname = name_k{i_k};
mbr_2 = mean_bar_width / 2; % half mean bar width
% For each group...
for i_g = 1:g
    stat = stat_g{i_g};
    animal_idx_m = [stat.animalidx];
    animal_idx_a = unique(animal_idx_m);
    a = numel(animal_idx_a);
    stat_values = zeros(1, a);
    % Plot the per-animal (i_a) mean
    for i_a = 1:a
        sel = animal_idx_m == animal_idx_a(i_a);
        astat = computeSummary(stat(sel));
        % Computer our stat, then save it for later
        stat_values(i_a) = computestat(strfields_f{i_f}, cname,...
            astat, type, minArea);
        % Plot an individual symbol
        doPlot(i_g, stat_values(i_a), style(i_a), 'Color', marker_color, ...
            'LineWidth', marker_line_width);
    end
    e = stderr(stat_values);
    m = mean(stat_values);

    % Plot mean an error bars
    doPlot(i_g + [-mbr_2, mbr_2], m * [1, 1], 'k', 'LineWidth', ...
        mbar_line_width);
    doPlot(i_g * [1, 1], m + [-e, e], 'k', 'LineWidth', ebar_line_width);    
end

set(gca, 'PlotBoxAspectRatio', pbox_aspect);

makeItPretty(type, i_k, i_f);

print(hf, '-dpng', sprintf('animal_plot_%s_%s_%s.png',...
    name_k{i_k}, strfields_f{i_f}, type));
print(hf, '-depsc', sprintf('animal_plot_%s_%s_%s.eps',...
    name_k{i_k}, strfields_f{i_f}, type));
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function makeItPretty(type, i_k, i_f)
global name_title_k strfield_title_f name_g label_font_size ...
    title_font_size axis_font_size g

pfstring = '';
switch type
    case 'a'
        pfstring = '%s - %s';
    case 'm'
        pfstring = '%s Mitochondria - %s';
end

title(sprintf(pfstring, name_title_k{i_k},...
    strfield_title_f{i_f}), 'FontSize', title_font_size,...
    'FontWeight', 'bold');

ylabel(strfield_title_f{i_f}, 'FontSize', label_font_size);

set(gca, 'XTick', 1:g);
set(gca, 'XTickLabel', name_g);
set(gca, 'FontSize', axis_font_size);

xlim([0 g + 1]);
ylim([0 max(ylim)]);

set(gca, 'Box', 'off');
set(gca, 'TickLength', [0 0]);

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function p = doPlot(varargin)
global marker_size
ph = plot(varargin{:});
set(ph, 'Markers', marker_size);
if nargout > 0
    p = ph;
end
end