

function y = forw_mod_2Dstack_v1( x,transp,traj_mat,SEM_mat,b1_mat,b1_plus,B1M_ind,use_gpu )
% function y = forw_mod_2Dstack_v1( x,transp,traj_mat,SEM_mat,b1_mat,b1_plus,B1M_ind,use_gpu )

numsamples = size(traj_mat,1);
[nx ny nz nechoes] = size(b1_plus);
[nx ny nz ncoils nrots] = size(b1_mat);

npix = nx*ny*nz;
npix_slice = nx*ny;

% flag for computing A*x ('no transp') or A'*y ('transp')
% this required for LSQR
if strcmp(transp,'transp')
    y = zeros( npix,1 );
else
    y = zeros( numsamples*nrots*ncoils*nz,1 );
end


% GPU flag
if use_gpu == 1
    x = gpuArray(x);
    y = gpuArray(y);
end


for z = 1:nz  % since slices are separated, more efficient to do a for loop at this level
              % than to use a block diagonal A matrix
              
    ind_pix_slice = 1 + (z-1)*npix_slice : z*npix_slice;

    % t1 = 0; t2 = 0; t3 = 0; t4 = 0;
    
    for rots = 1:nrots
        % tic;
        EXP = repmat( exp( -1j*traj_mat*reshape( SEM_mat(:,:,z,rots),1,npix_slice ) ),[ncoils 1] );
        % t1 = t1+ toc;

        % tic;
        B1M = repmat( reshape( b1_mat(:,:,z,:,rots),npix_slice,ncoils ).',[numsamples 1] );
        % t2 = t2 + toc;
        
        % tic;
        A = EXP .* B1M(B1M_ind,:);
        % t3 = t3 + toc;
        
        ind_start = 1 + (rots-1)*numsamples*ncoils + (z-1)*numsamples*ncoils*nrots;
        inds = ind_start : ind_start + numsamples*ncoils - 1;            
        
        % tic;
        if strcmp(transp,'transp')  % A'*y                
            y(ind_pix_slice) = y(ind_pix_slice) + A' * x(inds);
        else  % A*x
            y(inds) = A * x(ind_pix_slice);
        end        
        % t4 = t4 + toc;
        
    end
    
    % fprintf('T1 = %f  T2 = %f  T3 = %f  T4 = %f\n',t1/nrots,t2/nrots,t3/nrots,t4/nrots);
    
end


if use_gpu == 1
    y = gather(y);
end






















