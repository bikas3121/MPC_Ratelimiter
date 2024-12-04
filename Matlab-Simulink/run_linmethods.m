% main_mpc
% MPC Setup in YALMIP
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
Xcs_freq = 999;
Xcs_maxamp = Rng/2;
Xcs_offset = -Qstep/2;

switch 2
    case 1
        Nts = 1e5;
        Np = int(ceil(Xcs_freq*Ts*Nts));
    case 2
        Np = 3;
end

Npt = 1;
Npt = Np + Npt;
t_end = Np/Xcs_freq;
t = 0:Ts:t_end;
Xcs = (Xcs_Scale/100)*Xcs_maxamp*sin(2*pi*Xcs_freq*t) + Rng/2;


%% Direct Qunatization
C_DQ = floor(Xcs/Qstep+1/2);
switch Qmodel
    case 1
        Xcs_DQ = get_dac_output(C_DQ, YQns);
    case 2
        Xcs_DQ = get_dac_output(C_DQ, MLns);
end
  
%% Noise-shaping quantisation
b_nsf = b-a;
a_nsf = b;
C_NSD = noise_shaping(Nb, b_nsf, a_nsf, Qmodel, Qstep, Vmin, YQns, MLns, Xcs);
switch Qmodel
    case 1
        Xcs_NSD = get_dac_output(C_NSD, YQns);
    case 2
        Xcs_NSD = get_dac_output(C_NSD, MLns);
end




%% MPC steup 

N = 1;
x_init = [0,0]';
getControl = getControlMPC_BV(N, Nb, A, B, C, D, YQns, 'mosek');
C_MHOQ = MHOQ_BV(Xcs, N, getControl, x_init, Qmodel, YQns, MLns, A, B);

switch Qmodel
    case 1
        Xcs_MHOQ = get_dac_output(C_MHOQ, YQns);
    case 2
        Xcs_MHOQ = get_dac_output(C_MHOQ, MLns);
end


%% Switching minimisation
N = 1;
x_init = [0,0]';
beta = 0.01;
getControl_Fmin = getControlMHOQ_Fmin(N, Nb, A, B, C, D, YQns, beta, 'mosek');

C_MHOQ2 = MHOQ_FMIN(Xcs, N, Nb, getControl_Fmin, x_init, Qmodel, YQns, MLns, A, B);

switch Qmodel
    case 1
        Xcs_MHOQ2 = get_dac_output(C_MHOQ2, YQns);
    case 2
        Xcs_MHOQ2 = get_dac_output(C_MHOQ2, MLns);
end

%% Figure
sl = 100;
figure
plot(t(1:sl), Xcs_MHOQ(1:sl))
hold on 
plot(t(1:sl), Xcs_MHOQ2(1:sl))
grid on
legend("MHOQ1", "MHOQ2")



%% Post processing
[b1, a1] = butter(2, Wn);     % Butterworth Filter
f_ref = filter(b1,a1,Xcs);
f_ds = filter(b1, a1, Xcs_DQ);
f_nsd = filter(b1, a1, Xcs_NSD);
f_mhoq = filter(b1, a1, Xcs_MHOQ);
f_mhoq_fmin = filter(b1, a1, Xcs_MHOQ2);

var_dir = var(f_ref - f_ds)
var_nsd = var(f_ref - f_nsd )
var_mhoq = var(f_ref(1:length(f_mhoq)) - f_mhoq )
var_mhoq_fmin = var(f_ref(1:length(f_mhoq)) - f_mhoq_fmin )

figure
sinad(f_ds)
figure
sinad(f_nsd)
figure
sinad(f_mhoq)
figure
sinad(f_mhoq_fmin)
