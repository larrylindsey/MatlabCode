function plotError(eanalysis, ioki, nhx)


eraw = eanalysis.eraw;
etr = eanalysis.eblock(eanalysis.iok(ioki));

hx = linspace(0,.2, nhx);
hx = hx + .5 * (hx(2) - hx(1));

[hraw, hx] = hist(eraw.all, hx);
htr = hist(etr.all, hx);

hx(1) = [];
hraw(1) = [];
htr(1) = [];

plot(hx, hraw / sum(hraw), 'r', hx, htr / sum(htr), 'k', 'LineWidth', 2);
ylim([0 .1]);
set(gca, 'FontSize', 14);
grid on;

ylabel('Frequency', 'FontSize', 14);
xlabel('Alignment Error (pixels)', 'FontSize', 14);
set(legend('Raw Error', 'Model Reversed Error'), 'FontSize', 12);
