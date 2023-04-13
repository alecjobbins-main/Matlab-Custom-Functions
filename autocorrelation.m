function [R_tau, rho_tau, tau] = autocorrelation(x, fsamp, subtract_mean, biased)
    % Returns the autorrelation and autocorrelation coeff, unbiased or
    % biased
    N = length(x);
    if subtract_mean == "yes"
        if biased == "yes"
            R_tau = xcorr(x-mean(x))/N;
        elseif biased == "no"
            R_tau = xcorr(x-mean(x), 'unbiased');
        end
    elseif subtract_mean == "no"
        if biased == "yes"
            R_tau = xcorr(x)/N;
        elseif biased == "no"
             R_tau =xcorr(x, 'unbiased');
        end
    end
    m = -(N-1):1:(N-1);
    tau = m*(1/fsamp);
    index = find(tau==0);
    rho_tau = R_tau / R_tau(index);
end