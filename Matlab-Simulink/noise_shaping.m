% noise shaping quantisation 

function [C_NSD] = noise_shaping(Nb, b_nsf, a_nsf, Qmodel, Qstep, Vmin, YQns, MLns, Xcs)

    % Code storage container
    C_NSD = zeros(1, length(Xcs));
    
    % Buffer storage for error
    error_buffer = zeros(size(b_nsf));
    
    % Quantiser output
    u_ns = zeros(size(Xcs));
    
    for i = 1:length(Xcs)
    
        desired = Xcs(i) - b_nsf*error_buffer';
    
        q = floor(desired/Qstep +1/2);
    
        c = q- floor(Vmin/Qstep);
    
        c(c>2^Nb-1) = 2^Nb-1;
    
        c(c<0) = 0;

        C_NSD(1,i) = c;
    
        switch Qmodel
            case 1
                u_ns(1,i) = YQns(1, c+1);
            case 2
                u_ns(1,i) = MLns(1, c+1);
        end
    
        error = u_ns(i) - desired;
        error_buffer(1)  = - error;
        error_buffer(1)  = - a_nsf*error_buffer';
        error_buffer  = circshift(error_buffer,1);
    end
end