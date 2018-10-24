clear all;

% generates num_param parameters using Latin Hypercube Sampling
num_param = 1e4;
param_names = {'alpha_N', 'gamma_N', 'gamma_D', 'gamma_C', 'k_minus', 'k_plus'};
param_comb = lhsdesign(num_param, length(param_names))*1e2;

h = waitbar(0);
all_solutions = cell(1, num_param);
warning('OFF', 'optim:fsolve:NonSquareSystem');

% solves for equilibrium distribution of molecular species in 
% Model 0 for each parameter set in param_comb
for k = 1:size(param_comb, 1)      
    waitbar(k/size(param_comb, 1), h, [num2str(100*k/size(param_comb, 1)), '%']);
    solution_single_s = single_complex(param_comb(k, :));
    all_solutions{k} = solution_single_s;
end

close(h)

% IMPORTANT: solutions are saved to out_filename; modify as required 
out_filename = 'single_complex.mat';
save(out_filename, 'param_comb', 'all_solutions');

