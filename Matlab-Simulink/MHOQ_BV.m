% MHOQ impelementation using binary variables 

function [getCodes] = MHOQ_BV(Xcs, N, getControl, x_init, Qmodel, YQns, MLns, A, B)
    f = waitbar(0, 'Starting');
    C_MHOQ = [];
    len_MPC = length(Xcs)-N;
    for k = 1:len_MPC
        ref_in = Xcs(k:k+N-1) ;                          % intance t = k to t = k+N
        inputs = {x_init, ref_in};       % Inputs: Initial state and reference value at t = k
                                       
        var_opt1 = getControl(inputs);           % Optimal quantization levels at t = k
        var_opt = round(var_opt1{1});
    
        % get index of the variable
        indx = find(var_opt);
        C_MHOQ = [C_MHOQ indx-1];
    
        switch Qmodel
            case 1
                QL = YQns;
            case 2
                QL = MLns;
        end
    
        u_opt = QL*var_opt;
        x_new = A*x_init + B*(u_opt-ref_in(1));   % State prediction using obtained 
        
        x_init = x_new;  
        % u_st = [u_st; u_opt];
        waitbar(k/len_MPC, f, sprintf('Progress: %d %%', floor(k/len_MPC*100)));
    end
    close(f)
    getCodes = C_MHOQ;
end