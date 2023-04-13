function [p,CDF,Edges] = probability(x, N, dx, edges)
    figure();
    histogram(x,edges);
    hold on;
    Ni = floor(1.87*(N^(0.4)));
    plot(edges, Ni*ones(size(edges)))
    hold off;
    grid on;
    ylabel('Number of Occurrences');
    xlabel('y [m/s]');
    title('Histogram');
    legend('Histogram of Data', "Rule of Thumb Number of Occurrences (in 1+ bin(s)) = " + num2str(Ni));
    

    [Values, Edges] = histcounts(x, [min(x):dx:max(x)]);
    Edges = 0.5*(Edges(1:end-1) + Edges(2:end));
    figure();
    p = Values / (N*dx);
    trapz(Edges, p)
    plot(Edges, p);
    grid on;
    ylabel('p(x)');
    xlabel('x [m/s]');
    title('p(x) Approximation');

    figure();
    CDF = cumsum(Values)/N;
    plot(Edges, CDF);
    ylabel('CDF(y)');
    xlabel('y [m/s]');
    title('CDF');
    grid on;



end