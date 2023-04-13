function [filtered_signal] = signal_filtering(x, fsamp, fc, type, order)
    % x = unfiltered signal
    % fsamp = sampling frequency of unfiltered signal
    % fc = cutoff frequency in Hz.  [f_low f_high] if bandpass
    % type = type of filter: "RC Low Pass" or "RC High Pass" or "Ideal Low
    % Pass" or "Ideal High Pass" or "Ideal Band Pass"
    % order = order of filter
    N = length(x);
    k = 0:1:N-1;
    fk = fsamp*k/N;
    wk = 2*pi*fk;
    
    %% RC Filters
    if order == 1
        if type == "RC Low Pass"
            w0 = 2*pi*fc;
            G = 1./(1+(1i*wk/w0));
            filtered_signal = 2*real(ifft(G.*fft(x)));
        elseif type == "RC High Pass"
            w0 = 2*pi*fc;
            G = (1i.*(wk/w0))./(1+(1i*wk/w0));
            filtered_signal = 2*real(ifft(G.*fft(x)));
        end
    end

    %% Ideal Filters
    if type == "Ideal Low Pass"
        K = find(fk == fc);
        G = [];
        G(1:K) = 1;
        G(K+1:N) = 0;
        filtered_signal = 2*real(ifft(G.*fft(x)));
    elseif type == "Ideal High Pass"
        K = find(fk == fc);
        G = [];
        G(1:K-1) = 0;
        G(K:N) = 1;
        filtered_signal = 2*real(ifft(G.*fft(x)));
    elseif type == "Ideal Band Pass"
        K1 = find(fk == fc(1));
        K2 = find(fk == fc(2));
        G = [];
        G(1:K1-1) = 0;
        G(K1:K2) = 1;
        G(K2+1:N) = 0;
        filtered_signal = 2*real(ifft(G.*fft(x)));
    end

    %% Butter filters
    if type == "butter_bandpass" || type == "butter_low" || type == "butter_high"
        if type == "butter_bandpass"
            % fc should be in form [fc_low fc_high]
            Wn = fc./(fsamp/2);
            [b,a] = butter(order,Wn);
            filtered_signal = filter(b,a,x);
        elseif type == "butter_low"
            Wn = fc/(fsamp/2);
            [b,a] = butter(order,Wn,'low');
            filtered_signal = filter(b,a,x);
        elseif type == "butter_high"
            Wn = fc/(fsamp/2);
            [b,a] = butter(order,Wn,'high');
            filtered_signal = filter(b,a,x);
        end
    end
end