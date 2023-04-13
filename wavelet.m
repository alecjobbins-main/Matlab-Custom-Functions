function [t, wavelet_transform] = wavelet(x, fsamp, freq_array, b)
    %% Assumes complex morlet wavelet
    dt = 1/fsamp;
    N = length(x);
    T = N*dt;
    t = 0:dt:T-dt;
    k = [0:N-1];
    fk = fsamp*k/N;
    w = 2*pi*fk;
    kap_array = 2*pi*freq_array/b;
    G = zeros(length(kap_array), N);

    for i = 1:length(kap_array)
        k = kap_array(i);
        g_star = sqrt(2*pi)*(1-exp(-b*w/k)).*exp((-1/2)*((w/k)-b).^2);
        G(i,:) = ifft(fft(x).*g_star);
    end

    wavelet_transform = abs(G);
end