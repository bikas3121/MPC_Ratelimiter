clearvars; close all; clc;
format long

%% Parameters
Fs = 1e6;
Ts = 1/Fs;      

%% Quantiser configurations
Nb = 16;     % DAC bits
Vmax = 1;
Vmin = -1;
Qstep = (Vmax - Vmin)/(2^Nb-1);
YQ = Vmin:Qstep:Vmax;
Rng = Vmax - Vmin;

%% Reconstruction filter
Fc = 1e5;        %10000 without learning, 2500 with 1-bit learning, 125000 with 16-bit learning
Wc = Fc/(Fs/2);
[b,a] = butter(2, Wc);
[A, B, C, D] = tf2ss(b, a);


%% Reference signal 
t = 0:Ts:1;
XCS_maxamp = Rng/2- Qstep;
XCS_offset = -Qstep/2;
XCS_freq = 999;
XCS = XCS_maxamp*cos(2*pi*XCS_freq*t) + XCS_offset ;

IN = [t', XCS'];
%% Direct Quantisation
OUT_DIR = sim("direct.slx", t(end));
XCS_DIR = OUT_DIR.direct_out.Data;
%% Noise shaping quantiser

sr = 1e4; % slew rate
OUT_NSQ = sim("nsq_rate_lim.slx", t(end));
XCS_NSQ = OUT_NSQ.nsq_out.Data;
% Q_NSQ = OUT_NSQ.nsq_qout.Data;
%%
figure 
plot(t, XCS)
hold on 
plot(t, XCS_DIR)
hold on 
plot(t, XCS_NSQ)
ylim([Vmin-0.1, Vmax+0.1])
legend('Reference', ' Quantised')

%% Performance analysis
 
[SINAD_DIR, ENOB_DIR] = sinad_enob(XCS_DIR)
[SINAD_NSQ, ENOB_NSQ] = sinad_enob(XCS_NSQ)



