function [t,x] = inverse_spectra(c,fsamp,type,spectrum_domain)
    % Returns a signal given its fourier transform, ie. resturns the
    % inverse fourier transform.  
    % type = the type of fourier transform supplied ('discrete' or
    % 'continuous')
    % spectrum_domain = if a continuous FT is supplied, was the entire
    % continuous FT supplied (for positive and negative frequencies) or only
    % the positive domain?  If a discrete was supplied, was it supplied
    % with padded zeros at the end, or was it a normal dft?
    dt = 1/fsamp;
    if type == "discrete"
        if spectrum_domain == "both"
            N = length(c);
            j = 0:1:N-1;
            t = j*dt;
            x = N*ifft(c);
        elseif spectrum_domain == "positiveOnly" % assumes padded zeros at end...rarely use this
            N = length(c);
            j = 0:1:N-1;
            t = j*dt;
            x = 2*N*real(ifft(c));
        end
    elseif type == "continuous"
        if spectrum_domain == "both"
            % assumes you have a bell curve about 0
            ck = ifftshift(c);
            x = fsamp*ifft(ck);
            N = length(x);
            j = 0:1:N-1;
            t = j*dt;
        elseif spectrum_domain == "positiveOnly" 
            % assumes zeros are on negative axis and you only have the
            % positive axis information (ie. truncated half the info).
            ck = ifftshift(c);
            x = 2*fsamp*real(ifft(ck));
            N = length(x);
            j = 0:1:N-1;
            t = j*dt;
        elseif spectrum_domain == "positiveOnly_DFT_format"
            % assumes continuous FT that is already in DFT format, such as
            % Sxx.  ie. no need to perform ifftshift()
            x = fsamp*(ifft(c));
            N = length(x);
            j = 0:1:N-1;
            t = j*dt;
        end
    end
end