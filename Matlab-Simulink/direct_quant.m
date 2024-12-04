% Perform the direct quantisation of the provided signall 
function [C] = direct_quant(Xcs, Qstep, Vmin, Vmax)
    % Return the vector of directly quantiser signal
    dq_Xcs = floor(Xcs/Qstep+1/2)*Qstep;
    dq_Xcs(dq_Xcs > Vmax) = Vmax;
    dq_Xcs(dq_Xcs < Vmin) = Vmin;
    dq_Xcs = dq_Xcs/Qstep;
    C = dq_Xcs - floor(Vmin/Qstep);
end