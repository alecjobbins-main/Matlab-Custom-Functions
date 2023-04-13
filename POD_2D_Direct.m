function [a, Phi] = POD_2D_Direct(xloc, yloc, u, compression)
    My = size(u,1);
    Mx = size(u,2);
    N = size(u,3);
    
    u2 = reshape(u, My*Mx, N); % N (time) columns
    M = size(u2,1); % rows of u
    
    R = (u2 * u2') /  N; % Autocorrelation Matrix
    ds = (xloc(2)-xloc(1))*(yloc(2)-yloc(1)); % uniform grid spacing
    Rw = R*ds; % Rw matrix
    
    [eigvect, eigvalues] = eig(Rw);
    
    Phi(:,:) = eigvect(:,M:-1:1); % ASSUMES AN ORDER FROM MATLAB's eig function
    Phi = Phi / sqrt(ds);
    lambda = diag(eigvalues);
    Lambda = lambda(M:-1:1); % ASSUMES AN ORDER FROM MATLAB's eig function
    
    Energy = Lambda/sum(Lambda); % array of energies
    CumEnergy = cumsum(Energy); % array of cumulative energies
    
    figure();
    semilogx(CumEnergy, 'square');
    hold on;
    semilogx(Energy, 'o');
    hold off;
    grid on;
    xlabel('Mode #', 'FontSize', 16);
    ylabel('Fraction of Energy', 'FontSize', 16);
    title('Energy', 'FontSize', 16);
    legend('Cumulative Energy', 'Mode Energy', 'FontSize', 14);
    
    
    Phi2 = reshape(Phi, My, Mx, M);
    
    figure();
    contourf(xloc, yloc, Phi2(:,:,1), [-10:0.1:10], LineColor="none"); 
    colormap('jet');
    hold on;
    [C1, h1] = contour(xloc,yloc,Phi2(:,:,1),[-10:1:10], LineColor="k"); clabel(C1,h1);
    axis equal;
    xlabel('x [m]', 'FontSize',16);
    ylabel('y [m]', 'FontSize',16);
    title("Mode 1, (" + num2str(Energy(1)*100) + "%, " + num2str(CumEnergy(1)*100) + "%)", 'FontSize', 16);
    
    figure();
    contourf(xloc, yloc, Phi2(:,:,2),  [-10:0.1:10], LineColor="none"); 
    colormap('jet');
    hold on;
    [C2,h2] = contour(xloc, yloc, Phi2(:,:,2),  [-10:1:10], LineColor="k"); clabel(C2,h2);
    axis equal;
    xlabel('x [m]', 'FontSize',16);
    ylabel('y [m]', 'FontSize',16);
    title("Mode 2, (" + num2str(Energy(2)*100) + "%, " + num2str(CumEnergy(2)*100) + "%)", 'FontSize', 16);
    
    %% Check 1
    Orthonormality_mode1 = sum((Phi(:,1).^2)*ds)
    
    Orthonormality_mode2 = sum((Phi(:,2).^2)*ds)
    
    %% Check 2
    a = u2' * Phi * ds;
    
    figure();
    plot(a(:,1));
    hold on;
    plot(a(:,2));
    xlim([0 200]);
    grid on;
    legend('a1(t)', 'a2(t)');
    
    %% Check 3, reconstruct signal
    
    u_reconstructed = (a*(Phi'))';
    u_reconstructed_2D = reshape(u_reconstructed, My, Mx, N);
    
    figure();
    time_series_loc1 = u(20,3,:);
    time_series_loc1 = reshape(time_series_loc1, [1,N]);
    plot(time_series_loc1);
    hold on;
    time_series_loc1_reconstructed = u_reconstructed_2D(20,3,:);
    time_series_loc1_reconstructed = reshape(time_series_loc1_reconstructed, [1,N]);
    plot(time_series_loc1_reconstructed, '--');
    hold off;
    grid on;
    xlim([0 200]);
    xlabel('Time Index');
    ylabel('u');
    legend('Original Signal', 'POD Reconstruction');
    
    if compression ~= 0
        a = a(:, 1:compression);
        Phi = Phi(:, 1:compression);
        u_compressed = (a*(Phi'))';
        u_compressed_2D = reshape(u_compressed, My, Mx, N);
        u_compressed_x = u_compressed_2D(20,3,:);
        u_compressed_x = reshape(u_compressed_x, [1,N]);
    
        figure();
        plot(time_series_loc1);
        hold on;
        plot(u_compressed_x, '--');
        grid on;
        xlim([0 200]);
        xlabel('Index (not time)', 'FontSize', 16);
        ylabel('u', 'FontSize', 16);
        legend("Original Signal", "POD Compression Up To Mode " + num2str(compression), 'FontSize', 14);
        title("u(x,t) at x = #", 'FontSize', 16);
    end

end