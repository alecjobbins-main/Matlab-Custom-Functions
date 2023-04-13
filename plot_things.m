function [] = plot_things(x, y, signal_units, type)
    if type == "PSD"
        loglog(x, y);
        xlabel('f [Hz]', 'FontSize', 16, 'Interpreter','latex');
        ylabel("$S_{xx}$ $[\frac{(" + num2str(signal_units) + ")^2}{Hz}]$", 'Interpreter','latex','FontSize',16);
        title("Power Spectral Density", 'FontSize', 16);
        grid on;
    end
end

