

function y = pcg_fun_INTERF(x)
% function y = pcg_fun_INTERF(x)


y1 = forw_mod_fully3D_INTERF(x,'notransp');

y = forw_mod_fully3D_INTERF(y1,'transp');





