function [a, Phi] = POD_1D_NonUniformGrid_Snapshot(xloc, w, u, x, compression)
    %% Uses snapshot method, not direct method.  Good for M>>N (large space small time)
% xloc is a 1D array of xlocations of each 'sensor'
% u is the signal.  It is MxN where rows 1-M correspond to each position in
% xloc and the columns are time
% x is the location where you would like to fully reconstruct the signal as
% a check at the POD worked
% compression = 0 or the mode # you would like to compress up to
% a is the temporal coefficients
% phi is the eigen modes
    
    M = size(u,1); % rows of u
    N = length(u); % columns of u

    w1 = w*ones(1,N);
    u_tilde = sqrt(w1).*u;

    R = (u_tilde' * u_tilde) /  N; % Autocorrelation matrix
    
    [eigvect, eigvalues] = eig(R);
    
    b(:,:) = eigvect(:,N:-1:1); % ASSUMES AN ORDER FROM MATLAB's eig function

    Psi = (1/N)*u_tilde*b;
    
    lambda = diag(eigvalues);
    Lambda = lambda(N:-1:1); % ASSUMES AN ORDER FROM MATLAB's eig function
    Psi_tilde = Psi ./ sqrt(Lambda'/N);
    Phi = Psi_tilde ./ sqrt(w1);

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
    
    figure();
    plot(xloc, Phi(:,1), '.-');
    hold on;
    plot(xloc, Phi(:,2), 'square-');
    hold on;
    plot(xloc, Phi(:,3), '*-');
    hold on;
    plot(xloc, Phi(:,4), 'o-');
    grid on;
    hold off;
    xlabel('x [m]', 'FontSize',16);
    title('Modes', 'FontSize', 16);
    legend("Mode 1 (" + num2str(Energy(1)) + ", " + num2str(CumEnergy(1)) + ")", ...
        "Mode 2 (" + num2str(Energy(2)) + ", " + num2str(CumEnergy(2)) + ")", ...
        "Mode 3 (" + num2str(Energy(3)) + ", " + num2str(CumEnergy(3)) + ")", ...
        "Mode 4 (" + num2str(Energy(4)) + ", " + num2str(CumEnergy(4)) + ")", ...
        'FontSize', 14);
    
    a = sqrt(Lambda'*N).*b;
    figure();
    plot(a(:,1));
    hold on;
    plot(a(:,2));
    legend('a1', 'a2')
    xlim([0 200]);
    grid on;
    hold off;
    title('Temporal Coefficients');
    
    % Reconstruct signal at specified location
    index = find(xloc == x); 
    u_x = u(index,:); % Actual signal at this location
    u_reconstructed = a*(Phi');
    u_reconstructed_x = u_reconstructed(:,index);
    
    figure();
    plot(u_x);
    hold on;
    plot(u_reconstructed_x, '--');
    hold off;
    grid on;
    xlim([0 200]);
    xlabel('Index (not time)', 'FontSize', 16);
    ylabel('u', 'FontSize', 16);
    legend('Original Signal', 'POD Full Reconstruction', 'FontSize', 14);
    title("u(x,t) at x = " + num2str(x), 'FontSize', 16);

    if compression ~= 0
        a = a(:, 1:compression);
        Phi = Phi(:, 1:compression);
        u_compressed = a*(Phi');
        u_compressed_x = u_compressed(:,index);

        figure();
        plot(u_x);
        hold on;
        plot(u_compressed_x, '--');
        grid on;
        xlim([0 200]);
        xlabel('Index (not time)', 'FontSize', 16);
        ylabel('u', 'FontSize', 16);
        legend("Original Signal", "POD Compression Up To Mode " + num2str(compression), 'FontSize', 14);
        title("u(x,t) at x = " + num2str(x), 'FontSize', 16);   
    end

end