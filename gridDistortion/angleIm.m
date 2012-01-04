function im = angleIm(rc_grid1, rc_diff1, rc_grid2, rc_diff2)

[sel1 sel2 rc_grid1] = rcMatch(rc_grid1, rc_grid2);

rc_grid = rc_grid1(sel1,:);
rc_diff1 = rc_diff1(sel1,:);
rc_diff2 = rc_diff2(sel2,:);

rc_orientation1 = atan2(rc_diff1(:,2), rc_diff1(:,1));
rc_orientation2 = atan2(rc_diff2(:,2), rc_diff2(:,1));

rc_angle = rc_orientation2 - rc_orientation1;

rc_angle = mod(rc_angle + pi, 2 * pi) - pi;

im = gridData2im(rc_grid, rc_angle) * 180 / pi;

end
