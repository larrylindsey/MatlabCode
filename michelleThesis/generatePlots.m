function generatePlots(stat)

if ~isfield(stat, 'group')
    error('Need group info. Did you run applyGroupType?');
end

% if ~isfield(stat, 'calibrated')
%     error('Did you run calibrateStat?');
% end

global name_k strfields_f name_title_k strfield_title_f name_id_map ...
    strfields_id_map f k doMito_k g minArea name_g fig_font_size ...
    title_font_size hist_n

name_k = {'terminal', 'glia', 'capillary', 'background', 'allMito', ...
    'terminal_s', 'terminal_m', 'terminal_l', 'glia_s', 'glia_m', 'glia_l'};
strfields_f = {'nnorm', 'fraction', 'avgsize'};

name_title_k = {'Terminal', 'Glia', 'Capillary', 'Background', ...
    'Mitochondria',...
    'Small Terminals', 'Medium Terminals', 'Large Terminals', ...
    'Small Glia', 'Medium Glia', 'Large Glia'};
strfield_title_f = {'# / Area (um^{-2})', 'Fraction', 'Size (um^2)'};

name_g = {'YV', 'YE', 'YEP', 'AV', 'AE', 'AEP'};

hist_n = 32;

k = numel(name_k);
f = numel(strfields_f);

minArea = 256;

name_id_map = struct;
strfields_id_map = struct;

doMito_k = [true, true, false, false, false, false, false, false, false, false, false, false];

for i_k = 1:k
    name_id_map.(name_k{i_k}) = i_k;
end

for i_f = 1:f
    strfields_id_map.(strfields_f{i_f}) = i_f;
end

% cal_f = [1 1 .002 * .002];

fig_font_size = 14;
title_font_size = 16;

g = 6;
groups = [stat.group];
g_idx = unique(groups);

stat_g = cell(1, g);

for i_g = 1:g
    gsel = groups == g_idx(i_g);
    stat_g{i_g} = stat(gsel);
end

% plotDataByGroup(stat_g);
% plotDataByAnimalAndGroup(stat_g);
plotHistogramsByGroup(stat_g);

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
global name_k name_g hist_n strfield_title_f name_title_k title_font_size...
    fig_font_size g
i_f = 3;

figure;
cname = name_k{i_k};

stat_all = [cstat_g.(cname)];
data_all = [stat_all.(type)];
[~, hist_x] = hist(log10(data_all), hist_n);

for i_g = 1:g
    data = cstat_g(i_g).(cname).(type);
    
    subplot(2, 3, i_g);
    
    hist(log10(data), hist_x);
    title(sprintf('%s - %s', name_g{i_g}, name_title_k{i_k}), ...
        'FontSize', title_font_size);
%     xTicks = get(gca, 'XTick');
%     xTickLabels = 10.^xTicks;
%     set(gca, 'XTickLabel', xTickLabels);
    xlabel(['Log_{10} ' strfield_title_f{i_f}], 'FontSize', fig_font_size);
    
    
    set(gca, 'FontSize', fig_font_size);
    
end
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
global strfields_f minArea g name_k
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
            'k', 'LineWidth', 2);
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
global g strfields_f name_k minArea

mean_bar_width = .25;

style = 'oxvds^*';

hf = figure;
hold on;
cname = name_k{i_k};
mbr_2 = mean_bar_width / 2; % half mean bar width
for i_g = 1:g
    stat = stat_g{i_g};
    animal_idx_m = [stat.animalidx];
    animal_idx_a = unique(animal_idx_m);
    a = numel(animal_idx_a);
    stat_values = zeros(1, a);
    for i_a = 1:a
        sel = animal_idx_m == animal_idx_a(i_a);
        astat = computeSummary(stat(sel));
        stat_values(i_a) = computestat(strfields_f{i_f}, cname,...
            astat, type, minArea);
        
        plot(i_g, stat_values(i_a), style(i_a), 'LineWidth', 2);
    end
    e = stderr(stat_values);
    m = mean(stat_values);

    plot(i_g + [-mbr_2, mbr_2], m * [1, 1], 'k', 'LineWidth', 2);
    plot(i_g * [1, 1], m + [-e, e], 'k', 'LineWidth', 2);    
end
makeItPretty(type, i_k, i_f);

print(hf, '-dpng', sprintf('animal_plot_%s_%s_%s.png',...
    name_k{i_k}, strfields_f{i_f}, type));
print(hf, '-depsc', sprintf('animal_plot_%s_%s_%s.eps',...
    name_k{i_k}, strfields_f{i_f}, type));
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function makeItPretty(type, i_k, i_f)
global name_title_k strfield_title_f name_g fig_font_size title_font_size

pfstring = '';
switch type
    case 'a'
        pfstring = '%s - %s';
    case 'm'
        pfstring = '%s Mitochondria - %s';
end

title(sprintf(pfstring, name_title_k{i_k},...
    strfield_title_f{i_f}), 'FontSize', title_font_size);

ylabel(strfield_title_f{i_f}, 'FontSize', fig_font_size);

set(gca, 'XTick', 1:6);
set(gca, 'XTickLabel', name_g);
set(gca, 'FontSize', fig_font_size);

xlim([0 7]);

end
