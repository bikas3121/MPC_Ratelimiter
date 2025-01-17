% function [getControl]  = getControlMHOQ_RL(N, Nb, A, B, C, D, YQns, beta, solverstr)

%% Inputs and Outputs
% Inputs: 
    % N - Prediction horizon in optimization 
    % Q - Quantization levels (Constrained set for optimization variables)
    % A, B, C, D - System Matrices
    % solverstr - Solver name.  Can be chosen as follows, 
        % 'gurobi' - Gurobi Optimizer
        % 'mosek' - MOSEK Solver.

 % Outputs: 
    % getControl - returns optimal control (quantization levels) at each
    % prediction horizon ,
%% Problem Setup 
%% Sampling rate and frequency
Fs = 1e6;       % Sampling frequency
Ts = 1/Fs;  % Sampling rate
Qmodel = 1;
Qconfig = 1;
[Nb, Mq, Vmin, Vmax, Rng, Qstep, YQ]  = quant_config(Qconfig);
YQns = YQ;
b = [1.000000000000000e+00 ,   -1.892962366058354e-01 ,    1.390289638645996e-01];
a = [ 1.000000000000000e+00 ,   -1.142985437884713e+00  ,   4.128062411723096e-01];
Fc = 1e5;       % Cutoff frequency
Wn = Fc/(Fs/2); % Normalised frequency
[A,B,C,D]= tf2ss(b,a);      % Transfer function to state space

Xcs_Scale = 100;
Xcs_freq = 1e3;
Xcs_maxamp = Rng/2;
Xcs_offset = -Qstep/2;

Np = 1;
Npt = 1;
Npt = Np + Npt;
t_end = Np/Xcs_freq;
t = 0:Ts:t_end-Ts;
Xcs = (Xcs_Scale/100)*Xcs_maxamp*sin(2*pi*Xcs_freq*t) + Rng/2;
solverstr = "gurobi";

%% MPC Setup
Q = 0:1:2^Nb-1;
nu = 1; % dimension of the control
nx = 2; % dimension of the state
ny = 1; % dimension of the reference signal
N = 1;
%%
% Optimization Variables
% u = intvar(repmat(nu,1,N), repmat(1,1,N));      % control 

bv  = binvar(repmat(length(Q),1,N+1),repmat(1,1,N+1));
x = sdpvar(repmat(nx,1,N+1), repmat(1,1,N+1));  % state
r = sdpvar(repmat(ny,1,N+1),repmat(1,1,N+1));       % reference
bv2  = binvar(repmat(length(Q),1,N+1),repmat(1,1,N+1));


%% 
%%
%Constraints and Objective Functions

Constraints = [];
Objective = 0;

% Sum over prediction horizon 
for i = 1:N
    % Control
    con  = YQns * bv{i};
    % Error output
    et = C*x{i} + D* (con-r{i});   % output/error
    % Objective
    Objective = Objective + et'*et;  % Objective function 

    Constraints = [Constraints, x{i+1} == A*x{i} + B* (con -r{i})]; % State evolution constraints


    Constraints = [Constraints, sum(bv{i}) == 1];          % Input constraints
    
    Sw2 = bv{i} - bv2{1};
    Ns = Sw2'*Sw2;
    % Objective = Objective + beta*Ns;
end
% Parameters input to the optimizer
% Initial state and reference state at each prediction horizon
params_in = {x{1},[r{1:N}], bv2{1}};  

% Output of the optimizer
params_out = bv;

% Optimization setting
options = sdpsettings('verbose',0, 'solver', solverstr , 'showprogress',1);
getControl = optimizer(Constraints, Objective, [], params_in , params_out);

% end