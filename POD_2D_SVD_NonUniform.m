function [a, Phi, Phi2, Energy, CumEnergy] = POD_2D_SVD_NonUniform(w, u, compression, fsamp)
    My = size(u,1);
    Mx = size(u,2);
    N = size(u,3);

    u_tilde = sqrt(w).*u;
    u2 = reshape(u_tilde, My*Mx, N); % N (time) columns
    M = size(u2,1); % rows of u2
    [Psi, Sigma, V] = svd(u2, 'econ');
    
    lambda = (1/N)*Sigma*Sigma;
    Lambda = diag(lambda);
    Energy = Lambda/sum(Lambda); % array of energies
    CumEnergy = cumsum(Energy); % array of cumulative energies
    
    
    w1 = reshape(w, My*Mx, 1);
    w2 = w1 * ones(1,N);
    Phi = Psi ./ sqrt(w2);
    Phi2 = reshape(Phi, My, Mx, N);


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
  
  
    
    %% Check 1
    Orthonormality_mode1 = sum((Psi(:,1).^2))
    
    Orthonormality_mode2 = sum((Psi(:,2).^2))

    Orthonormality_mode3 = sum((Psi(:,3).^2))

    NormalCheck1 = sum((Psi(:,1).*Psi(:,2)))

    NormalCheck2 = sum((Psi(:,1).*Psi(:,3)))

    NormalCheck3 = sum((Psi(:,1).*Psi(:,4)))


    
    %% Check 2
    a = (V*Sigma);
    dt = 1/fsamp;
    T = (N-1)*dt;
    t = 0:dt:T;
    figure();
    plot(t, a(:,1), 'r--');
    hold on;
    plot(t, a(:,2), 'b.-');

    %xlim([0 200]);
    grid on;
    legend('a_1(t)', 'a_2(t)', 'FontSize', 16);
    xlabel('Time [s]', 'FontSize', 16);
    ylabel('Coeff. Value', 'FontSize',16);
    title('Temporal Coefficients', 'FontSize', 16);

    %% Check 2.5 Plot Spectra of temporal coefficients
    a1 = a(:,1);
    a2 = a(:,2);

    window = 'hanning';
    Num_blocks = 1;
    points = 0;
    overlapping = 'no';
    phase_correction = 'no';
    energyLB = 0;
    energyUB = 0;
    truncateAfterNyquist = 'yes';

    [f_k_a1, c_k, c_k_real, c_k_imag, A, C_continuous_approximation, ...
    C_continuous_approximation_real, C_continuous_approximation_imag, Sxx_a1,...
    sigma_squared_based_on_x, sigma_squared_based_on_Sxx, ...
    N_block, Num_blocks] ...
    = spectra(a1, fsamp, Num_blocks, overlapping, points, window, ...
    phase_correction, energyLB, energyUB, truncateAfterNyquist);

    [f_k_a2, c_k, c_k_real, c_k_imag, A, C_continuous_approximation, ...
    C_continuous_approximation_real, C_continuous_approximation_imag, Sxx_a2,...
    sigma_squared_based_on_x, sigma_squared_based_on_Sxx, ...
    N_block, Num_blocks] ...
    = spectra(a2, fsamp, Num_blocks, overlapping, points, window, ...
    phase_correction, energyLB, energyUB, truncateAfterNyquist);

    figure();
    loglog(f_k_a1, Sxx_a1);
    hold on;
    loglog(f_k_a2, Sxx_a2);
    grid on;
    xlabel('f [Hz]', 'FontSize', 16);
    ylabel('$S_{xx}$ $[\frac{[x^2]}{Hz}]$', 'Interpreter','latex','FontSize',16);
    title('S_{xx} of Temporal Coefficients', 'FontSize', 16);
    legend('PSD of a_1(t)', 'PSD of a_2(t)', 'FontSize', 14);


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