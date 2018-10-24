# elife2017-cis_activation_modeling
MATLAB code for calculating equilibrium concentrations of molecular species in cis-activation models, for a range of parameters

For Model 0, run 'scan_single_params.m', which solves the model equations in 'single_complex.m'
For Models 1, and 2a-2d, run 'scan_params.m', which solves model equations in 'double_complex.m'
IMPORTANT: Equations corresponding to model of interest have to be uncommented within double_complex.m. Note that sampled parameters
are saved along with solutions in the output file, and can be reused for other models for consistency. 

Equations for models including both cis- and trans-interactions are contained in 'cis_and_trans.m'. Run 'scan_params_with_trans.m' 
to solve equations for a range of parameter sets. 
IMPORTANT: The variable 'trans_cis_ratio' in 'scan_params_with_trans.m' sets the relative strength of trans- and cis-interactions
