


function plot_precon()
% function plot_precon()

global reconStruct_small;

A = gpuArray( zeros(reconStruct_small.npix) );
for i=1:reconStruct_small.npix
    fprintf('\tComputing A matrix %d/%d ...\n',i,reconStruct_small.npix);
    
    x = zeros(reconStruct_small.npix,1);
    x(i) = 1.0;
    
    A(:,i) = pcg_fun_INTERF(x);
end

A = gather(A);

P_img = diag(reconStruct_small.P);
figure; imagesc(abs([ A P_img A-P_img ])); axis image off; colorbar;
title('A   P   A-P');


eval( sprintf('save A_%dslices.mat A Pimg;',reconStruct_small.nz) );
