function out = cis_and_trans(var_vals, trans_cis_ratio)

    global alpha_N gamma_N alpha_D gamma_D gamma_Cact gamma_Cinh kact_minus kact_plus kinh_minus kinh_plus 
    global T_kact_minus D_trans
   
    p = num2cell(var_vals);
    [alpha_N, gamma_N, gamma_D, gamma_Cact, gamma_Cinh, kact_minus, kact_plus, kinh_minus, kinh_plus] = deal(p{:}); 
            
    % assume alpha_D has a 100-fold range, maxing at alpha_N
    alpha_D_vals = logspace(log10(var_vals(1)/10), log10(10*var_vals(1)), 10);
    
    % First solve model without trans interactions to calculate range of
    % cis-Delta (D_cis) values for particular choice of parameters
    temp = struct;
    for k = 1:length(alpha_D_vals)
        alpha_D  = alpha_D_vals(k);
        options = optimoptions('fsolve', 'Display', 'off');
        init = [alpha_N/gamma_N alpha_D/gamma_D 1 1];
        X = fsolve(@cis_eqn_generate, init, options);
        temp(k).N = X(1); temp(k).D = X(2); temp(k).Cact = X(3); temp(k).Cinh = X(4); temp(k).alpha_D = alpha_D;
    end
    
    % relative dissociation constant of trans Notch-Delta complex, compared to cis
    % complex
    T_kact_minus = kact_minus/trans_cis_ratio;
    
    % Next, assume trans-Delta values (D_trans) spans the same range as D_cis
    D_cis = [temp.D]; 
    D_trans_vals = logspace(log10(D_cis(1)), log10(D_cis(end)), 10);
    
    % Solve model, now including both cis- and trans- interactions
    out = struct;
    for k = 1:length(D_trans_vals)
        for p = 1:length(alpha_D_vals)
            alpha_D = alpha_D_vals(p);
            D_trans = D_trans_vals(k);
            options = optimoptions('fsolve', 'Display', 'off');
            init = [alpha_N/gamma_N alpha_D/gamma_D 1 1 0];
            X = fsolve(@ct_eqn_generate, init, options);
            out(k, p).N = X(1); out(k, p).D = X(2); out(k, p).Cact = X(3); out(k, p).Cinh = X(4); out(k, p).T = X(5); out(k, p).alpha_D = alpha_D;      
        end
    end
    
end


function eqns = cis_eqn_generate(X)
    
    global alpha_N gamma_N alpha_D gamma_D gamma_Cact gamma_Cinh kact_minus kact_plus kinh_minus kinh_plus
    N = X(1); D = X(2); Cact = X(3); Cinh = X(4);
    
    %% 3) Model 2c: non-stoichiometric ND interactions for clustering
    eqn_N = alpha_N - gamma_N*N + (kact_minus*Cact - kact_plus*N*D);
    eqn_D = alpha_D - gamma_D*D + (kact_minus*Cact - kact_plus*N*D) + (kinh_minus*Cinh - kinh_plus*Cact*D);
    eqn_Cact =  - gamma_Cact*Cact + (kact_plus*N*D - kact_minus*Cact)  + (kinh_minus*Cinh - kinh_plus*Cact*D);
    eqn_Cinh =  - gamma_Cinh*Cinh + kinh_plus*Cact*D - kinh_minus*Cinh;

    eqns = [eqn_N; eqn_D; eqn_Cact; eqn_Cinh];
end

function eqns = ct_eqn_generate(X)
    
    global alpha_N gamma_N alpha_D gamma_D gamma_Cact gamma_Cinh kact_minus kact_plus kinh_minus kinh_plus
    global T_kact_minus D_trans
    
    N = X(1); D = X(2); Cact = X(3); Cinh = X(4); T = X(5);
    
    %% 3) Model 2c with trans interactions
    eqn_N = alpha_N - gamma_N*N + (kact_minus*Cact - kact_plus*N*D) + (T_kact_minus*T - kact_plus*N*D_trans);
    eqn_D = alpha_D - gamma_D*D + (kact_minus*Cact - kact_plus*N*D) + (kinh_minus*Cinh - kinh_plus*Cact*D);
    eqn_Cact =  - gamma_Cact*Cact + (kact_plus*N*D - kact_minus*Cact)  + (kinh_minus*Cinh - kinh_plus*Cact*D);
    eqn_Cinh =  - gamma_Cinh*Cinh + kinh_plus*Cact*D - kinh_minus*Cinh;
    eqn_T = - gamma_Cact*T + (kact_plus*N*D_trans - T_kact_minus*T);
 
    eqns = [eqn_N; eqn_D; eqn_Cact; eqn_Cinh; eqn_T];
end

