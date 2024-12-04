function [getControl] = getControlMPC_BV(N, Nb, A, B, C, D, YQns, solverstr)
% This function forumlates the optimization problem using the binary
% variables to solve the mixed integer programming problem. 
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

%% MPC Setup
Q = 0:1:2^Nb-1;
nu = 1; % dimension of the control
nx = 2; % dimension of the state
ny = 1; % dimension of the reference signal

%%
% Optimization Variables
% u = intvar(repmat(nu,1,N), repmat(1,1,N));      % control 

bv  = binvar(repmat(length(Q),1,N+1),repmat(1,1,N+1));
x = sdpvar(repmat(nx,1,N+1), repmat(1,1,N+1));  % state
r = sdpvar(repmat(ny,1,N+1),repmat(1,1,N+1));       % reference

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
end
% Parameters input to the optimizer
% Initial state and reference state at each prediction horizon
params_in = {x{1},[r{1:N}]};  

% Output of the optimizer
params_out = bv;

% Optimization setting
options = sdpsettings('verbose',0, 'solver', solverstr , 'showprogress',1);
getControl = optimizer(Constraints, Objective, [], params_in , params_out);


end


