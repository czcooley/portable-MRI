

function y = forw_mod_fully3D_v1( x,transp,traj_mat,SEM_mat,b1_mat,b1_plus,B1M_ind,use_gpu )
% function y = forw_mod_fully3D_v1( x,transp,traj_mat,SEM_mat,b1_mat,b1_plus,B1M_ind,use_gpu )


numsamples = size(traj_mat,1);
[nx ny nz nechoes] = size(b1_plus);
[nx ny nz ncoils nrots] = size(b1_mat);

npix = nx*ny*nz;
npix_slice = nx*ny;

% flag for computing A*x ('no transp') or A'*y ('transp')
% required for LSQR
if strcmp(transp,'transp')
    y = zeros( npix,1 );
else
    y = zeros( numsamples*nz*nrots*ncoils,1 );
end


% GPU flag
if use_gpu == 1
    x = gpuArray(x);
    y = gpuArray(y);
end

% for rots = 1:nrots   
%     fprintf('rot = %d\n',rots);
%     EXP = exp( -1j*traj_mat*reshape( SEM_mat(:,:,:,rots),1,npix ) );
%     
%     for coil = 1:ncoils
%         B1M = repmat( reshape(b1_mat(:,:,:,coil,rots),1,npix),[numsamples 1] );
%         
%         for ee=1:nz
%             B1P = repmat( reshape(b1_plus(:,:,:,ee),1,npix),[numsamples 1] );
%             
%             ind_start = 1 + (rots-1)*numsamples + (coil-1)*numsamples*nrots + (ee-1)*numsamples*nrots*ncoils;
%             inds = ind_start : ind_start + numsamples-1;
%             
%             A = EXP .* B1P .* B1M;
%             
%             if strcmp(transp,'transp')  % A'*y                
%                 y = y + A' * x(inds);
%             else  % A*x
%                 y(inds) = A * x;
%             end
%             
%         end
%     end
% end


for rots = 1:nrots   
    EXP = repmat( exp( -1j*traj_mat*reshape( SEM_mat(:,:,:,rots),1,npix ) ),[ncoils 1] );    
    B1M = repmat( reshape(b1_mat(:,:,:,:,rots),npix,ncoils).',[numsamples 1] );    
    B1M = B1M(B1M_ind,:);
    
    for ee=1:nz
        B1P = repmat( reshape(b1_plus(:,:,:,ee),1,npix),[numsamples*ncoils 1] );
        
        % B1P_B1M = reshape(b1_mat(:,:,:,:,rots),npix,ncoils) .* repmat( reshape(b1_plus(:,:,:,ee),npix,1),[1 ncoils] );
        % B1P_B1M = repmat( B1P_B1M.',[numsamples 1] );
        
        ind_start = 1 + (rots-1)*numsamples*ncoils + (ee-1)*numsamples*ncoils*nrots;
        inds = ind_start : ind_start + numsamples*ncoils -1;
        
        A = EXP .* B1P .* B1M;
        % A = EXP .* B1P_B1M(B1M_ind,:);
        
        if strcmp(transp,'transp')  % A'*y                
            y = y + A' * x(inds);
        else  % A*x
            y(inds) = A * x;
        end
        
    end
    
end

if use_gpu == 1
    y = gather(y);
end






















