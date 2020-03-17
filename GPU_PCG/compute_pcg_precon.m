


function P = compute_pcg_precon()
% function P = compute_pcg_precon()


global reconStruct_small;


B1M = (reconStruct_small.b1m_re).^2 + (reconStruct_small.b1m_im).^2;
B1P = (reconStruct_small.b1p_re).^2 + (reconStruct_small.b1p_im).^2;

P = squeeze( sum( sum(B1M,3),2 ) ) .* sum( B1P,2 ) .* reconStruct_small.nsamples;

