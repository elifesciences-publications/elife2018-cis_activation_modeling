clear all;

% generates num_param parameters using Latin Hypercube Sampling
num_param = 1e4;
param_names = {'alpha_N', 'gamma_N', 'gamma_D', 'gamma_Cact', 'gamma_Cinh', 'kact_minus', 'kact_plus', 'kinh_minus', 'kinh_plus'};
param_comb = lhsdesign(num_param, length(param_names))*1e2;

h = waitbar(0);
all_solutions = cell(1, num_param);
warning('OFF', 'optim:fsolve:NonSquareSystem');

% solves for equilibrium distribution of molecular species in the model of
% interest (Model 1, Models 2a-2d) for each parameter set in param_comb
%
% IMPORTANT: Uncomment equations corresponding to model of interest in
% the 'eqn_generate' function within the 'double_complex.m' file

for k = 1:size(param_comb, 1)      
    waitbar(k/size(param_comb, 1), h, [num2str(100*k/size(param_comb, 1)), '%']);
    solution_double_s = double_complex(param_comb(k, :));
    all_solutions{k} = solution_double_s;
end

close(h)

% IMPORTANT: solutions are saved to out_filename; modify as required 
out_filename = 'double_parallel.mat';
save(out_filename, 'param_comb', 'all_solutions');
