


function img = RECON( reconStruct,niters,tol,display,freqs,freq_weights )
% function img = RECON( reconStruct,niters,tol,display,freqs,freq_weights )

global reconStruct_small;

% close all;

% % TEST -- REMOVE IMAGE PIXELS FOR QUICK PRECONDITIONER ANALYSIS
% us_fac = 2;
% reconStruct.b1_mat = reconStruct.b1_mat(1:us_fac:end,1:us_fac:end,:,:,:);
% reconStruct.b1_plus = reconStruct.b1_plus(1:us_fac:end,1:us_fac:end,:,:);
% reconStruct.SEM_mat = reconStruct.SEM_mat(1:us_fac:end,1:us_fac:end,:,:);
% 
% reconStruct.XX = reconStruct.XX(1:us_fac:end,1:us_fac:end);
% reconStruct.YY = reconStruct.YY(1:us_fac:end,1:us_fac:end);
% reconStruct.ZZ = reconStruct.ZZ(1:us_fac:end,1:us_fac:end);
% 
% reconStruct.reconSize = reconStruct.reconSize(1:us_fac:end,1:us_fac:end);
% 
% 
% % TEST -- REMOVE HALF THE SAMPLES TO FIT DATA IN RECON
% reconStruct.recondata = reconStruct.recondata(1:2:end,:,:,:);
% reconStruct.numsamples = 128;
% reconStruct.traj_mat = reconStruct.traj_mat(1:2:end);
% reconStruct.time = reconStruct.time(1:2:end);


% PREPARE GPU
reconStruct_small.display = display;
reconStruct_small.freqs = freqs;
reconStruct_small.freq_weights = freq_weights;
reconStruct_small.nfreqs = size(freqs,2);

prepare_gpu( reconStruct );


% OPTIMIZE
x0 = single( zeros(reconStruct_small.npix,1) );


% % FFT + 2D LSQR -- this is not needed as the fully 3D recon already
% converges quickly
% DATA_ss = separate_trase_slices_fft( reconStruct.recondata );
% 
% DATA_ss = permute( DATA_ss,[1 3 2 4] );  % permute so that the "coil" dimension is the 2nd fastest varying
%                                          % data dimensions are samples * coils * rots * echoes
% 
% 
% % TEST -- 1 slice
% DATA_ss = DATA_ss(:,:,:,1);
% 
%                                          
% fun = @(x,transp)forw_mod_2Dstack( x,transp, ...
%     reconStruct.traj_mat, ...
%     reconStruct.SEM_mat, ...
%     reconStruct.b1_mat, ...
%     reconStruct.b1_plus, ...
%     B1M_ind, ...
%     use_gpu );
% 
% x = lsqr_( fun,DATA_ss(:),tol,niters,[],[],x0 );


DATA = permute( single(reconStruct.recondata),[1 3 2 4] );  % permute so that the "coil" dimension is the 2nd fastest varying
                                                            % data dimensions are samples * coils * rots * echoes

                                                            
% FULLY 3D PCG (normal equations)
fprintf('Compute preconditioner ...\n');
reconStruct_small.P = compute_pcg_precon();  % diagonal preconditioner

% plot_precon();  % plot preconditioner to compare to the full A matrix

fprintf('Transform RHS (normal equations) ...\n');
DATA_allfreqs = zeros( reconStruct_small.ndata_points*reconStruct_small.nfreqs,1 );
for f = 1:reconStruct_small.nfreqs
    inds = (f-1)*reconStruct_small.ndata_points + 1: f*reconStruct_small.ndata_points;
    DATA_allfreqs( inds ) = DATA(:) * sqrt( reconStruct_small.freq_weights(f) );
end
b = forw_mod_fully3D_INTERF(DATA_allfreqs,'transp');  % needed for normal equation

if reconStruct_small.use_gpu == 1
    b = gpuArray(b);
    x0 = gpuArray(x0);
    reconStruct_small.P = gpuArray( reconStruct_small.P );
end

fprintf('\n********** FULLY 3D PCG SOLVE (%d FREQS.)  **********\n',reconStruct_small.nfreqs);

x = pcg_( @(x)pcg_fun_INTERF(x),b,tol,niters,@(x)mat_vec_mult_PCG_precon(x),[],x0 );  % with precon
% x = pcg_( @(x)pcg_fun_INTERF(x),b,tol,niters,[],[],x0 );  % without precon


img = gather( reshape(x,reconStruct_small.nx,reconStruct_small.ny,reconStruct_small.nz) );
















