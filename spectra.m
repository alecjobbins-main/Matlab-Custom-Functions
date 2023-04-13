function [f_k, c_k, c_k_real, c_k_imag, A, C_continuous_approximation, ...
    C_continuous_approximation_real, C_continuous_approximation_imag, Sxx,...
    sigma_squared_based_on_x, sigma_squared_based_on_Sxx, ...
    N_block, Num_blocks] ...
    = spectra(x, fsamp, Num_blocks, overlapping, points, window, ...
    phase_correction, energyLB, energyUB, truncateAfterNyquist)
    %% Returns DFT, Approximation to continuous FT, Power Spectrum, and Energy
    % based on a signal
    % [f_k, c_k, c_k_real, c_k_imag, A, C_continuous_approximation, ...
    % C_continuous_approximation_real, C_continuous_approximation_imag, Sxx,...
    % sigma_squared_based_on_x, sigma_squared_based_on_Sxx, ...
    % N_block, Num_blocks] ...
    % = spectra(x, fsamp, Num_blocks, points, window, phase_correction, energyLB, energyUB)

    % Based on inputs, determine important parameters
    if points == 0 % If care more about how many blocks there are than how many points
        dt = 1/fsamp;
        N_block = floor(length(x) / Num_blocks); % Number of points per block
        k = 0:N_block-1; % mode number vector
        TB = (N_block)*dt; % Time of each block
        t=k*dt; % time vector for each block
    else % If care more about how many points per block
        dt = 1/fsamp;
        N_block = points; % Number of points per block
        Num_blocks = floor(length(x) / N_block);
        k = 0:N_block-1; % mode number vector
        TB = (N_block)*dt; % Time of each block
        t=k*dt; % time vector for each block
    end
    
    % Truncate signal if data overflows blocks, reshape signal into blocks
    if overlapping == "yes"
        Num_blocks = Num_blocks + (Num_blocks-1);
        x_new = [];
        for i = 1:Num_blocks
            if i == 1
                x_new = x(1:N_block);
            else 
                x_new = [x_new; [x( (((i-1)/2)*N_block)+1 : ((i+1)/2)*N_block )]];
            end
        end
        x = x_new;
        blocks = reshape(x,N_block, Num_blocks);
    elseif overlapping == "no"
        x = x(1:Num_blocks*N_block);
        blocks = reshape(x,N_block, Num_blocks);
    end
    
    % Define Window to use
    if window == "no"
        scalingConstant = 1;
        win = ones(Num_blocks, N_block);
        win = win';
    elseif window == "hanning"
        scalingConstant = sqrt(8/3);
        hanningWindow = 0.5*(1-cos(2*pi*t/(TB)));
        win = hanningWindow.*ones(Num_blocks, N_block);
        win = win';
    end
    
    % Compute DFT, Continuous Approximation of FT, Power Spectra, and Energy
    f_k = fsamp*(k/N_block); 

    if phase_correction == "no"
        c_k = scalingConstant*fft(blocks.*win)/(N_block);
        c_k_real = real(c_k);
        c_k_imag = imag(c_k);
        C_continuous_approximation = abs(c_k) * (TB);
        C_continuous_approximation_real = real(c_k * (TB));
        C_continuous_approximation_imag = imag(c_k * (TB));
    elseif phase_correction == "yes"
        c_k = (scalingConstant.*fft(blocks.*win)/(N_block));
        phase_corrector = exp(2*pi*1i*f_k*(TB+dt)/2)';
        c_k = phase_corrector.*c_k;
        c_k_real = real(c_k);
        c_k_imag = imag(c_k);
        C_continuous_approximation = abs(c_k) * (TB);
        C_continuous_approximation_real = real(c_k * (TB));
        C_continuous_approximation_imag = imag(c_k * (TB));
    end

    Sxx = (abs(c_k)).^2*(N_block / fsamp);
    Sxx = mean(Sxx,2);
    
    % Return arrays only up to nyquist frequency
    if truncateAfterNyquist == "yes"
         f_nyquist = fsamp / 2;
         f_k(f_k > f_nyquist) = []; %drop all freq greater than nyquist
    
         c_k = abs(c_k);
         c_k = c_k(1:length(f_k));
         A = 2*c_k;
         c_k_real = c_k_real(1:length(f_k));
         c_k_imag = c_k_imag(1:length(f_k));
 
         C_continuous_approximation = C_continuous_approximation(1:length(f_k));
         C_continuous_approximation_real = C_continuous_approximation_real(1:length(f_k));
         C_continuous_approximation_imag = C_continuous_approximation_imag(1:length(f_k));
 
         Sxx = Sxx(1:length(f_k));
    elseif truncateAfterNyquist == "no" 
         c_k = abs(c_k);
         A = 2*c_k; 
    end

    if energyLB == energyUB % return energy under entire spectra
        % Energy based on entire (truncated) signal
        sigma_squared_based_on_x = std(x)^2; % includes negative energy
        %sigma_squared_based_on_x = (1/N)*sum(x.^2);
        delta_f = 1 / TB;
        % Energy based on integral under power spectra (rectange rule)
        sigma_squared_based_on_Sxx = 2*delta_f*sum(Sxx); % includes imaginary energy
        %sigma_squared_based_on_Sxx = 2*delta_f*sum(Sxx(1:floor(N_block/2))); % less accurate
    else % return energy in desired frequency range
        % Energy based on entire (truncated) signal
        sigma_squared_based_on_x = 0; % irrelevant statistic in this case
        delta_f = 1 / TB;
        % Energy based on integral under power spectra (rectange rule)
        lower_index = find(f_k>=energyLB,1);
        upper_index = find(f_k>=energyUB,1);
        sigma_squared_based_on_Sxx = 2*delta_f*sum(Sxx(lower_index:upper_index)); % includes imag. energy
        %sigma_squared_based_on_Sxx = 2*delta_f*sum(Sxx(1:floor(N_block/2))); % less accurate
    end
    
end

