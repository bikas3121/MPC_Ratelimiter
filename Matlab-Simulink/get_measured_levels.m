% Dac measurements;  Simulates the  measuremed of the DAC levels

function [INL] = get_measured_levels(Qconfig)
    switch Qconfig
        case 1
            if exist("Measured_INL/INL_levels_04bit.csv", 'file')
                INL = readmatrix("Measured_INL/INL_levels_04bit.csv");
            else
                disp('File not found')
            end
        case 2
            if exist("Measured_INL/INL_levels_06bit.csv", 'file')
                INL = readmatrix("Measured_INL/INL_levels_06bit.csv");
            else
                disp('File not found')
            end
        case 3
            if exist("Measured_INL/INL_levels_08bit.csv", 'file')
                INL = readmatrix("Measured_INL/INL_levels_08bit.csv");
            else
                disp('File not found')
            end
        case 4
            if exist("Measured_INL/INL_levels_12bit.csv", 'file')
                INL = readmatrix("Measured_INL/INL_levels_12bit.csv");
            else
                disp('File not found')
            end
        case 5
            if exist("Measured_INL/INL_levels_16bit.csv", 'file')
                INL = readmatrix("Measured_INL/INL_levels_16bit.csv");
            else
                disp('File not found')
            end
        otherwise
            disp('Invalid quantiser configuration')
    end 
end