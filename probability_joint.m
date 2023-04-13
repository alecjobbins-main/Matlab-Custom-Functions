function [] = probability_joint(x, y, dx, dy, edgesX, edgesY)
    figure();
    histogram2(x,y,edgesX, edgesY);
    grid on;
    xlabel('x [m/s]', 'FontSize', 16);
    ylabel('y [m/s]', 'FontSize', 16);
    title('Joint Histogram', 'FontSize', 14);
    

    [Counts, X2, Y2] = histcounts2(x, y, edgesX, edgesY);
    X2=X2(1:end-1) + dx/2;
    Y2=Y2(1:end-1) + dy/2;
    figure();
    [C, h2] = contour(X2,Y2,Counts'/(sum(Counts(:))*dx*dy), [10:10:100],'-b');
    clabel(C,h2);
    axis equal;
    grid on;
    xlabel('x [m/s]', 'FontSize', 16);
    ylabel('y [m/s]', 'FontSize', 16);
    title('p(x,y)', 'FontSize', 14);

end