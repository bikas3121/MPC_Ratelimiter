function [SINAD, ENOB] =sinad_enob(XCS)
% Calculate the SINAD and ENOB
    SINAD = sinad(XCS);
    ENOB = (SINAD - 1.76)/6.02;
end