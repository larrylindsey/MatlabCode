subplot(2,2,1);
pp = plot(pmax);
fs = 16;
grid on
set([xlabel('V_d_c [V]') ylabel('P_m_a_x [W]')],'FontSize', fs);
subplot(2,2,2)
pp(2) = plot([optstructsm.maxh]);
set([ylabel('h_0 [m]') xlabel('V_d_c [V]')], 'FontSize', fs);
grid on;
subplot(2,2,3)
pp(3) = plot([optstructsm.maxn]);
ylim([min(ylim) max(ylim)*1.1])
grid on
set([xlabel('V_d_c [V]') ylabel('n')], 'FontSize', fs);
subplot(2,2,4)
pp(4) = plot([optstructsm.maxm]);
grid on;
set([xlabel('V_d_c [V]') ylabel('m [kg]')], 'FontSize', fs);
set(pp, 'LineWidth', 2);