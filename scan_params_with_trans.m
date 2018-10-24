clear all;

% generates num_param parameters using Latin Hypercube Sampling
num_param = 1e4;
param_names = {'alpha_N', 'gamma_N', 'gamma_D', 'gamma_Cact', 'gamma_Cinh', 'kact_minus', 'kact_plus', 'kinh_minus', 'kinh_plus'};
param_comb = lhsdesign(num_param, length(param_names))*1e2;

with_trans_solutions = {};

% the variable trans_cis_ratio sets the relative dissociation constant of
% trans Notch-Delta complexes to cis Notch-Delta complexes (higher the
% ratio, weaker trans interactions compared to cis interactions)
trans_cis_ratio = 1;

h = waitbar(0);

successful_param_comb = param_comb; 
% IMPORTANT: in the paper, trans interactions are only considered for parameter
% combinations that give rise to non-monotonic cis-activation profiles,
% i.e.
% successful_param_comb = param_comb(successful_param_comb_indices, :);

for k = 1:size(successful_param_comb)
    waitbar(k/size(successful_param_comb, 1), h, [num2str(100*k/size(successful_param_comb, 1)), '%']);
    solution = cis_and_trans(successful_param_comb(k, :), trans_cis_ratio);
    with_trans_solutions{end + 1} = solution;
end
close(h);

% IMPORTANT: solutions are saved to out_filename; modify as required 
out_filename = ['with_trans_ratio', num2str(trans_cis_ratio), '.mat'];
save(out_filename, 'with_trans_solutions', 'successful_param_comb');
