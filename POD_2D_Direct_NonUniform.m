function [a, Phi, Phi2, Energy, CumEnergy] = POD_2D_Direct_NonUniform(w, u, compression, fsamp)
    My = size(u,1);
    Mx = size(u,2);
    N = size(u,3);

    u_tilde = sqrt(w).*u;
    u2 = reshape(u_tilde, My*Mx, N); % N (time) columns
    M = size(u2,1); % rows of u
    
    R = (u2 * u2') /  N; % Autocorrelation Matrix  
    [eigvect, eigvalues] = eig(R);
    
    Psi(:,:) = eigvect(:,M:-1:1);
    w1 = reshape(w, My*Mx, 1);
    w2 = w1 * ones(1,M);
    Phi = Psi ./ sqrt(w2);

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
    
  
    
    %% Check 1
    Orthonormality_mode1 = sum((Psi(:,1).^2))
    
    Orthonormality_mode2 = sum((Psi(:,2).^2))

    Orthonormality_mode3 = sum((Psi(:,3).^2))

    NormalCheck1 = sum((Psi(:,1).*Psi(:,2)))

    NormalCheck2 = sum((Psi(:,1).*Psi(:,3)))

    NormalCheck3 = sum((Psi(:,1).*Psi(:,4)))


    
    %% Check 2
    a = (u2' * Psi);
    dt = 1/fsamp;
    T = (N-1)*dt;
    t = 0:dt:T;
    figure();
    plot(t, a(:,1), 'r--');
    hold on;
    plot(t, a(:,2), 'b.-');
    hold on;
    plot(t, a(:,3));
    %xlim([0 200]);
    grid on;
    legend('a_1(t)', 'a_2(t)', 'a_3(t)', 'FontSize', 16);
    xlabel('Time [s]', 'FontSize', 16);
    ylabel('Coeff. Value', 'FontSize',16);
    title('Temporal Coefficients', 'FontSize', 16);

    
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
    %% Check 4: Plot Dominant Modes
    
end