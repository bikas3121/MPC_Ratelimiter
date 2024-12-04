% Simulate the static DAC 

function [DAC_out] = get_dac_output(C, Qlevels)
    % Get DAC output
    C_size = size(C);
    DAC_out = zeros(C_size(1),C_size(2));

    for j = 1:C_size(1)
        for i = 1:C_size(2)
            DAC_out(j,i) = Qlevels(C(i)+1);
        end
    end
end