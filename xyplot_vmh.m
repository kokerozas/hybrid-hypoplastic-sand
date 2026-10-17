% XYPLOT_VMH: plotting section (model in red, measured curves in blue when they are available)
%
% Strains are plotted in %, stresses in kPa; q = sigma_a - sigma_r (negative in extension).

CM = [0.80 0.15 0.15];   CE = [0.20 0.40 0.70];
ea = 100*EE(:,3);   ev = 100*INV_E(:,1);   p = INV_S(:,1);   q = SS(:,3) - SS(:,1);
has = ~isempty(DATA);
figure('Color', 'w', 'Position', [80 80 900 620]);
switch TEST
    case 'TCU'
        N = NLIQ.(['N_' NLIQ_RULE]);   ru = 1 - p/p0;
        if has, DATA.N = cycles_halfwave(DATA.q);   DATA.ru = 1 - DATA.p/DATA.p(1); end   % tests: half-wave count
        subplot(2,2,1), hold on, box on, grid on
        if has, plot(DATA.e1, DATA.q, 'Color', CE), end
        plot(ea, q, 'Color', CM), xlabel('axial strain \epsilon_1 (%)'), ylabel('deviator stress q (kPa)')
        subplot(2,2,2), hold on, box on, grid on
        if has, plot(DATA.p, DATA.q, 'Color', CE), end
        plot(p, q, 'Color', CM), xlabel('mean effective stress p'' (kPa)'), ylabel('q (kPa)'), xlim([0 inf])
        subplot(2,2,3), hold on, box on, grid on
        if has, plot(DATA.N, DATA.ru, 'Color', CE), end
        plot(N, ru, 'Color', CM), yline(0.95, ':k');
        xlabel('N (cycles)'), ylabel('r_u = 1 - p''/p''_0'), ylim([-0.1 1.05])
        subplot(2,2,4), hold on, box on, grid on
        if has, plot(DATA.N, DATA.e1, 'Color', CE), end
        plot(N, ea, 'Color', CM), xlabel('N (cycles)'), ylabel('\epsilon_1 (%)')
        tit = sprintf('%s TCU %s   q_{amp} = %.2f kPa   N_{liq} = %.3f', SAND, num2str(ID), q_amp, Nliq);
    case {'TMD', 'TMU', 'TMU2'}
        subplot(2,2,1), hold on, box on, grid on
        if has, plot(DATA.e1, DATA.q, 'Color', CE), end
        plot(ea, q, 'Color', CM), xlabel('\epsilon_1 (%)'), ylabel('q (kPa)')
        subplot(2,2,2), hold on, box on, grid on
        if has, plot(DATA.p, DATA.q, 'Color', CE), end
        plot(p, q, 'Color', CM), xlabel('p'' (kPa)'), ylabel('q (kPa)')
        subplot(2,2,3), hold on, box on, grid on
        if has && isfield(DATA, 'ev'), plot(DATA.e1, DATA.ev, 'Color', CE), end
        plot(ea, ev, 'Color', CM), xlabel('\epsilon_1 (%)'), ylabel('\epsilon_v (%)')
        subplot(2,2,4), hold on, box on, grid on
        plot(ea, HARD(:,1), 'Color', CM), xlabel('\epsilon_1 (%)'), ylabel('void ratio e')
        tit = sprintf('%s %s %s   p_0 = %.1f kPa   e_0 = %.3f', SAND, TEST, num2str(ID), p0, e0);
    case 'OED'
        e = HARD(:,1);   s1 = SS(:,3);
        if has                                       % unloading drawn down to the last measured stress
            [~, jm] = max(s1);   s_end = max(15, min(DATA.s1(find(DATA.s1 == max(DATA.s1), 1):end)));
            keep = (1:numel(s1))' <= jm | s1 >= s_end;   e = e(keep);   s1 = s1(keep);
        end
        hold on, box on, grid on
        if has, plot(DATA.s1, DATA.e, 'Color', CE), end
        plot(s1, e, 'Color', CM), set(gca, 'XScale', 'log')
        xlabel('vertical stress \sigma_1 (kPa)'), ylabel('void ratio e')
        tit = sprintf('%s OED %s   e(%g kPa) = %.4f', SAND, num2str(ID), p0, e0);
end
sgtitle([tit '   (red: model, blue: test)'], 'FontSize', 10)
