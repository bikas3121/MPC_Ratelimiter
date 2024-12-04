yalmip('clear')
close all
clear
clc
format long

%% Sampling rate and frequency
Fs = 1e6;       % Sampling frequency
Ts = 1/Fs;  % Sampling rate

%% Quantiser model
% 1 = Ideal
% 2 = Non ideal (with INL)
Qmodel = 1;

%% Qunatization Parameters 

Qconfig = 1;
[Nb, Mq, Vmin, Vmax, Rng, Qstep, YQ]  = quant_config(Qconfig);
INL = get_measured_levels(Qconfig)*Qstep;

YQns = YQ; % Ideal levels;
MLns = YQ + INL; % Measured levels

%% Reconstruction  Filter Parameters

Fc = 1e5;       % Cutoff frequency
Wn = Fc/(Fs/2); % Normalised frequency
switch 2
    case 1
        [b, a] = butter(2, Wn);     % Butterworth Filter

    case 2   % Optimal filter H2-Hinf optimal 
        b = [1.000000000000000e+00 ,   -1.892962366058354e-01 ,    1.390289638645996e-01];
        a = [ 1.000000000000000e+00 ,   -1.142985437884713e+00  ,   4.128062411723096e-01];
end

[A,B,C,D]= tf2ss(b,a);      % Transfer function to state space

%% Reference Signal
Xcs_Scale = 100;
Xcs_freq = 1e3;
Xcs_maxamp = Rng/2;
Xcs_offset = -Qstep/2;

switch 2
    case 1
        Nts = 1e5;
        Np = int(ceil(Xcs_freq*Ts*Nts));
    case 2
        Np = 1;
end

Npt = 1;
Npt = Np + Npt;
t_end = Np/Xcs_freq;
t = 0:Ts:t_end-Ts;
Xcs = (Xcs_Scale/100)*Xcs_maxamp*sin(2*pi*Xcs_freq*t) + Rng/2;


%% Direct Qunatization
C_DQ = floor(Xcs/Qstep+1/2);

% Static Model
switch Qmodel
    case 1
        Xcs_DQ1 = get_dac_output(C_DQ, YQns);
    case 2
        Xcs_DQ1 = get_dac_output(C_DQ, MLns);
end

sr = 1e4; % slew rate


IN_DQ = [t', Xcs_DQ1'];
OUT_DAC_DQ = sim("DACwithSlewRate.slx", t(end));
Xcs_DQ2 = OUT_DAC_DQ.out_dac_dq.Data;
%% Resampling the signal
Xcs_DQ1_re = datasample(Xcs, length(Xcs)*10);
Ts_re = 1/1e7;
t_re = 0:Ts_re:t_end-Ts_re;

figure
% plot (t, Xcs)
% hold on 
plot (t_re, Xcs_DQ1_re')
%%
figure
plot (t, Xcs)
hold on 
plot (t, Xcs_DQ1)
hold on 
plot (t, Xcs_DQ2)
legend("Ref","DQ1", "DQ2")
grid on 

figure
sinad(Xcs_DQ1)
figure
sinad(Xcs_DQ2)
%% Noise shaping
% b_nsf = b-a;
% a_nsf = b;
% C_NSD = noise_shaping(Nb, b_nsf, a_nsf, Qmodel, Qstep, Vmin, YQns, MLns, Xcs);
% switch Qmodel
%     case 1
%         Xcs_NSD1 = get_dac_output(C_NSD, YQns);
%     case 2
%         Xcs_NSD1 = get_dac_output(C_NSD, MLns);
% end
% 
% %% MPC 
% N = 1;
% x_init = [0,0]';
% 
% getControl = getControlMPC_BV(N, Nb, A, B, C, D, YQns, 'mosek');
% C_MHOQ = MHOQ_BV(Xcs, N, getControl, x_init, Qmodel, YQns, MLns, A, B);
% switch Qmodel
%     case 1
%         Xcs_MHOQ1 = get_dac_output(C_MHOQ, YQns);
%     case 2
%         Xcs_MHOQ1 = get_dac_output(C_MHOQ, MLns);
% end
% 
% 
% %% MPC with switching minimisation
% N = 1;
% x_init = [0,0]';
% beta = 0.01;
% getControl_Fmin = getControlMHOQ_Fmin(N, Nb, A, B, C, D, YQns, beta, 'mosek');
% C_MHOQ2 = MHOQ_FMIN(Xcs, N, Nb, getControl_Fmin, x_init, Qmodel, YQns, MLns, A, B);
% 
% switch Qmodel
%     case 1
%         Xcs_MHOQ_Fmin1 = get_dac_output(C_MHOQ2, YQns);
%     case 2
%         Xcs_MHOQ_Fmin1 = get_dac_output(C_MHOQ2, MLns);
% end
% 
% 
% %% Simulink Models
% 
% sr = 1e5; % slew rate
% 
% 
% IN_DQ = [t', Xcs_DQ1'];
% IN_NSD = [t', Xcs_NSD1'];
% IN_MHOQ = [t(1:length(Xcs_MHOQ1))', Xcs_MHOQ1'];
% IN_MHOQ_Fmin = [t(1:length(Xcs_MHOQ_Fmin1))', Xcs_MHOQ_Fmin1'];
% 
% 
% 
% OUT_DAC_DQ = sim("DACwithSlewRate.slx", t(end));
% % OUT_DAC_NSD = sim("DACwithSlewRate.slx", t(end));
% OUT_DAC_MHOQ = sim("DACwithSlewRate.slx", t(end));
% OUT_DAC_MHOQ_Fmin = sim("DACwithSlewRate.slx", t(end));
% 
% 
% Xcs_DQ2 = OUT_DAC_DQ.out_dac_dq.Data;
% % Xcs_NSD2 = OUT_DAC_NSD.out_dac_nsd.Data;
% Xcs_MHOQ2 = OUT_DAC_MHOQ.out_dac_mhoq.Data;
% Xcs_MHOQ_Fmin2 = OUT_DAC_MHOQ_Fmin.out_dac_mhoq_Fmin.Data;
% 
% 
% %% Filtering and SINAD
% [b,a] = butter(2,Fc/(Fs/2));
% FXcs_DQ1 = filter(b,a,Xcs_DQ1);
% FXcs_DQ2 = filter(b,a,Xcs_DQ2);
% 
% % FXcs_NSD1 = filter(b,a,Xcs_NSD1);
% % FXcs_NSD2 = filter(b,a,Xcs_NSD2);
% % 
% FXcs_MHOQ1 = filter(b,a,Xcs_MHOQ1);
% FXcs_MHOQ2 = filter(b,a,Xcs_MHOQ2);
% 
% FXcs_MHOQ_FMin1 = filter(b,a,Xcs_MHOQ_Fmin1);
% FXcs_MHOQ_Fmin2 = filter(b,a,Xcs_MHOQ_Fmin2);
% 
% figure
% sinad(FXcs_MHOQ1)
% figure
% sinad(FXcs_MHOQ2)
% 
% figure
% sinad(FXcs_MHOQ_FMin1)
% figure
% sinad(FXcs_MHOQ_Fmin2)
% 
% 
% 
% %% Performance analysis
% sl = 1000;
% figure
% plot(t(1:sl), Xcs(1:sl))
% hold on 
% plot(t(1:sl), Xcs_MHOQ_Fmin1(1:sl))
% hold on 
% plot(t(1:sl), Xcs_MHOQ_Fmin1(1:sl))
% grid on 
% legend('Ref','FMin1', 'Fmin2')
% 
% 
