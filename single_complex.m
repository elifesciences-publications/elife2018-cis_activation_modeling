function out = single_complex(var_vals)

    global alpha_N gamma_N alpha_D gamma_D gamma_C k_minus k_plus 
   
    p = num2cell(var_vals);
    [alpha_N, gamma_N, gamma_D, gamma_C, k_minus, k_plus] = deal(p{:}); 
        
    % assume alpha_D has a 100-fold range, maxing at alpha_N
    alpha_D_vals = logspace(log10(var_vals(1)/10), log10(10*var_vals(1)), 10);
    
    out = struct;
    for k = 1:length(alpha_D_vals)
        alpha_D  = alpha_D_vals(k);
        options = optimoptions('fsolve', 'Display', 'off');
        init = [alpha_N/gamma_N alpha_D/gamma_D 1 1];
        X = fsolve(@eqn_generate, init, options);
        out(k).N = X(1); out(k).D = X(2); out(k).C = X(3); out(k).alpha_D = alpha_D;
    end
end

function eqns = eqn_generate(X)
    
    global alpha_N gamma_N alpha_D gamma_D gamma_C k_minus k_plus 
    N = X(1); D = X(2); C = X(3); 
    
    eqn_N = alpha_N - gamma_N*N + (k_minus*C - k_plus*N*D);
    eqn_D = alpha_D - gamma_D*D + (k_minus*C - k_plus*N*D);
    eqn_C = - gamma_C*C + (k_plus*N*D - k_minus*C);
     
    eqns = [eqn_N; eqn_D; eqn_C];

end
