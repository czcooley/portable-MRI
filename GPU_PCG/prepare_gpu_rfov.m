

function reconStruct_small = prepare_gpu( reconStruct )
% function reconStruct_small = prepare_gpu( reconStruct )

global reconStruct_small;

reconStruct_small.use_gpu = 0;

[nx ny nz ncoils nrots] = size(reconStruct.SEM_mat();
npix = nx*ny*nz;    

reconStruct_small.npix = npix;
reconStruct_small.ncoils = ncoils;
reconStruct_small.nrots = nrots;
reconStruct_small.nsamples = size(reconStruct.traj_mat,1);
reconStruct_small.nechoes = nz;

reconStruct_small.nx = nx;
reconStruct_small.ny = ny;
reconStruct_small.nz = nz;

reconStruct_small.ndata_points = numel(reconStruct.recondata);

if parallel.gpu.GPUDevice.isAvailable       
    
    reconStruct_small.use_gpu = 1;
    
    fprintf('Grab GPU ...\n');
    reconStruct_small.g = gpuDevice;
    reset( reconStruct_small.g);
           
    % A * x kernel
    fprintf('Preparing A*x kernel ...\n');
    reconStruct_small.forw_mod_fully3D_cudaKernel_v1 = parallel.gpu.CUDAKernel('forw_mod_fully3D_cudaKernel_v4.ptx','forw_mod_fully3D_cudaKernel_v4.cu');
    
    M = 128;    
    N=ceil( reconStruct_small.ndata_points * reconStruct_small.nfreqs / M );
    
    reconStruct_small.forw_mod_fully3D_cudaKernel_v1.ThreadBlockSize = M;
    reconStruct_small.forw_mod_fully3D_cudaKernel_v1.GridSize = N;
        
    reconStruct_small.traj_mat = gpuArray( single( reconStruct.traj_mat ) );
    reconStruct_small.SEM_mat = gpuArray( single( reshape(reconStruct.SEM_mat,npix,nrots) ) );
    
    reconStruct_small.b1m_re = gpuArray( single( real(reshape(reconStruct.b1_mat,npix,ncoils,nrots)) ) );
    reconStruct_small.b1m_im = gpuArray( single( imag(reshape(reconStruct.b1_mat,npix,ncoils,nrots)) ) );
    
    reconStruct_small.b1p_re = gpuArray( single( real(reshape(reconStruct.b1_plus,npix,nrots)) ) );
    reconStruct_small.b1p_im = gpuArray( single( imag(reshape(reconStruct.b1_plus,npix,nrots)) ) );
    
    reconStruct_small.b1_plus_echo_exp = gpuArray( single(reconStruct.b1_plus_echo_exp) );
    
    reconStruct_small.y_re = gpuArray( single(zeros(reconStruct_small.ndata_points * reconStruct_small.nfreqs,1)) );
    reconStruct_small.y_im = gpuArray( single(zeros(reconStruct_small.ndata_points * reconStruct_small.nfreqs,1)) );
    
    
    % PE field
    reconStruct_small.PE_field = gpuArray( single( reshape(reconStruct.grad_mat,npix,nrots) ) );

    
    % A' * y kernel
    fprintf('Preparing A''*y kernel ...\n');
    reconStruct_small.forw_mod_fully3D_cudaKernel_TRANSP_v1 = parallel.gpu.CUDAKernel('forw_mod_fully3D_cudaKernel_TRANSP_v4.ptx','forw_mod_fully3D_cudaKernel_TRANSP_v4.cu');
    
    M = 128;    
    N=ceil( npix / M );

    reconStruct_small.forw_mod_fully3D_cudaKernel_TRANSP_v1.ThreadBlockSize = M;
    reconStruct_small.forw_mod_fully3D_cudaKernel_TRANSP_v1.GridSize = N;

    reconStruct_small.x_re = gpuArray( single(zeros(npix,1)) );
    reconStruct_small.x_im = gpuArray( single(zeros(npix,1)) );

    
else
    
    reconStruct_small.traj_mat = single( reconStruct.traj_mat );
    reconStruct_small.SEM_mat = single( reshape(reconStruct.SEM_mat,npix,nrots) );    
    reconStruct_small.b1m = single( reshape(reconStruct.b1_mat,npix,ncoils,nrots) );    
    reconStruct_small.b1p = single( reshape(reconStruct.b1_plus,npix,nrots) );
    reconStruct_small.b1_plus_echo_exp = single( b1_plus_echo_exp );
    
end













