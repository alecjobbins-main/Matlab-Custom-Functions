function [R_xy, rho_xy, tau] = cross_correlation(x, y, fsamp, time_delay, biased)
    N = length(x);
    if time_delay == "no"
        R_xy = mean(x.*y);
        rho_xy = R_xy / (std(x)*std(y)); % pearson coefficient
        tau = 0;
    elseif time_delay == "yes"
        if biased == "yes"
            R_xy = conj(xcorr(y,x,'biased'));
            rho_xy = R_xy / (std(x)*std(y)); 
        elseif biased == "no"
            R_xy = conj(xcorr(y,x,'unbiased'));
            rho_xy = R_xy / (std(x)*std(y)); 
        end
        m = -(N-1):1:(N-1);
        tau = m*(1/fsamp);
    end
end