function out = double_complex(var_vals)

    global alpha_N gamma_N alpha_D gamma_D gamma_Cact gamma_Cinh kact_minus kact_plus kinh_minus kinh_plus
   
    p = num2cell(var_vals);
    [alpha_N, gamma_N, gamma_D, gamma_Cact, gamma_Cinh, kact_minus, kact_plus, kinh_minus, kinh_plus] = deal(p{:}); 
        
    % assume alpha_D has a 100-fold range, maxing at alpha_N
    alpha_D_vals = logspace(log10(var_vals(1)/10), log10(10*var_vals(1)), 10);
    
    out = struct;
    for k = 1:length(alpha_D_vals)
        alpha_D  = alpha_D_vals(k);
        options = optimoptions('fsolve', 'Display', 'off');
        init = [alpha_N/gamma_N alpha_D/gamma_D 1 1];
        X = fsolve(@eqn_generate, init, options);
        out(k).N = X(1); out(k).D = X(2); out(k).Cact = X(3); out(k).Cinh = X(4); out(k).alpha_D = alpha_D;

    end

end

% IMPORTANT: uncomment equations corresponding to model of interest
function eqns = eqn_generate(X)
    
    global alpha_N gamma_N alpha_D gamma_D gamma_Cact gamma_Cinh kact_minus kact_plus kinh_minus kinh_plus
    N = X(1); D = X(2); Cact = X(3); Cinh = X(4);

    
%% 0) Model 1: no clustering, parallel tracks for inhibitory and activating complexes
%     eqn_N = alpha_N - gamma_N*N + (kact_minus*Cact - kact_plus*N*D) + (kinh_minus*Cinh - kinh_plus*N*D);
%     eqn_D = alpha_D - gamma_D*D + (kact_minus*Cact - kact_plus*N*D) + (kinh_minus*Cinh - kinh_plus*N*D);
%     eqn_Cact = - gamma_Cact*Cact + (kact_plus*N*D - kact_minus*Cact);
%     eqn_Cinh = - gamma_Cinh*Cinh + (kinh_plus*N*D - kinh_minus*Cinh);
    
%% 1) Model 2a: trimolecular rxn for clustering
%     eqn_N = alpha_N - gamma_N*N + (kact_minus*Cact - kact_plus*N*D) + (kinh_minus*Cinh - kinh_plus*N*D*Cact);
%     eqn_D = alpha_D - gamma_D*D + (kact_minus*Cact - kact_plus*N*D) + (kinh_minus*Cinh - kinh_plus*N*D*Cact);
%     eqn_Cact = - gamma_Cact*Cact + (kact_plus*N*D - kact_minus*Cact) + (kinh_minus*Cinh - kinh_plus*N*D*Cact);
%     eqn_Cinh = - gamma_Cinh*Cinh + (kinh_plus*N*D*Cact - kinh_minus*Cinh);
   
%% 2) Model 2b: bimolecular complex merging for clustering
%     eqn_N = alpha_N - gamma_N*N + (kact_minus*Cact - kact_plus*N*D);
%     eqn_D = alpha_D - gamma_D*D + (kact_minus*Cact - kact_plus*N*D);
%     eqn_Cact =  - gamma_Cact*Cact + (kact_plus*N*D - kact_minus*Cact)  + (2*kinh_minus*Cinh - kinh_plus*Cact^2);
%     eqn_Cinh =  - gamma_Cinh*Cinh + kinh_plus*Cact^2 - kinh_minus*Cinh;

%% 3) Model 2c: non-stoichiometric ND interactions for clustering
%     eqn_N = alpha_N - gamma_N*N + (kact_minus*Cact - kact_plus*N*D);
%     eqn_D = alpha_D - gamma_D*D + (kact_minus*Cact - kact_plus*N*D) + (kinh_minus*Cinh - kinh_plus*Cact*D);
%     eqn_Cact =  - gamma_Cact*Cact + (kact_plus*N*D - kact_minus*Cact)  + (kinh_minus*Cinh - kinh_plus*Cact*D);
%     eqn_Cinh =  - gamma_Cinh*Cinh + kinh_plus*Cact*D - kinh_minus*Cinh;

%% 4) Model 2d: higher-order (2N + 2D -> Cinh) ND interactions for clustering
%     n = 2;
%     eqn_N = alpha_N - gamma_N*N + (kact_minus*Cact - kact_plus*N*D) + (n*kinh_minus*Cinh - kinh_plus*N^n*D^n);
%     eqn_D = alpha_D - gamma_D*D + (kact_minus*Cact - kact_plus*N*D) + (n*kinh_minus*Cinh - kinh_plus*N^n*D^n);
%     eqn_Cact =  - gamma_Cact*Cact + (kact_plus*N*D - kact_minus*Cact);
%     eqn_Cinh =  - gamma_Cinh*Cinh + (kinh_plus*N^n*D^n - kinh_minus*Cinh);
  
     eqns = [eqn_N; eqn_D; eqn_Cact; eqn_Cinh];

end