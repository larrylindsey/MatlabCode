function stat_n = calibrateStat(stat_n)

cal = .002 * .002; % .002 um/pixel, or that ^2 for pixels-squared.

n = numel(stat_n);

fn = fieldnames(stat_n);
f = numel(fn);

for i_n = 1:n
    for i_f = 1:f
        stat = stat_n(i_n);
        stat.totpix = stat.totpix * cal;
        if isstruct(stat.(fn{i_f}))
           stat.(fn{i_f}).a = stat.(fn{i_f}).a * cal;
           stat.(fn{i_f}).m = stat.(fn{i_f}).m * cal;
        end
        stat_n(i_n) = stat;
        stat_n(i_n).calibrated = true;
    end
end