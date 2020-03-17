

function y = mat_vec_mult_PCG_precon(x)
% function y = mat_vec_mult_PCG_precon(x)

global reconStruct_small;

y = x ./ reconStruct_small.P;







