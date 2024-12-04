% Quantiser configuration 
% Only midtread quantiser are considered
function [Nb, Mq, Vmin, Vmax, Rng, Qstep, YQ]  = quant_config(Qconfig)

switch Qconfig
    case 1
        Nb = 4;
        Mq = 2^Nb -1;
    case 2
        Nb = 6;
        Mq = 2^Nb -1;
    case 3
        Nb = 8;
        Mq = 2^Nb -1;
    case 4
        Nb = 12;
        Mq = 2^Nb -1;
    case 5
        Nb = 16;
        Mq = 2^Nb -1;
end
Vmax = 2;
Vmin = 0;
Rng = Vmax-Vmin; 
Qstep = Rng/Mq;
YQ = linspace(Vmin, Vmax, Mq+1);
end