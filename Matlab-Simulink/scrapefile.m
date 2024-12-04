% scrapefile
Fs = 1e4;
Ts = 1/Fs;
t_end = 0.001;
t = 0:Ts:t_end-Ts;

x = sin(2*pi*1e3*t);



IN = x;
OUT_upsampler = sim("upsampler.slx", t(end));
X_up = OUT_upsampler.out.Data;

X_us = []; 
for i = 1:length(x)
    X_us = [X_us; X_up(:,:,i)];
end

%%
Ts_up = 1/1e5;
t_up = 0:Ts_up:t_end-Ts_up;
figure
% plot(t, x)
% hold on 
plot(1:1:length(X_us), X_us)