function [t,x] = example_signal_from_spectra(Sxx, fsamp, N)
    % Assumes provided Spectra only up to f_Nyq = fsamp/2
    % Need to provide N (number of points in signal), ie. T
    T = N/fsamp;
    a_squiggle = 1/sqrt(2)*randn(1,N/2);
    b_squiggle = 1/sqrt(2)*randn(1,N/2);
    c_k_squiggle = a_squiggle + 1i*b_squiggle;

    ck = c_k_squiggle .* sqrt(Sxx / T);
    ck(N/2 + 1: N) = 0;

    x = 2*N*real(ifft(ck));

    j = 1:N;
    t = ((j-1)/N)*T;
end