

function y = forw_mod_fully3D_v2( x,transp,traj_mat,SEM_mat,b1_mat,b1_plus,B1M_ind,use_gpu )
% function y = forw_mod_fully3D_v2( x,transp,traj_mat,SEM_mat,b1_mat,b1_plus,B1M_ind,use_gpu )

global reconStruct_small;

use_gpu = 0;

% flag for computing A*x ('no transp') or A'*y ('transp')
% required for LSQR
if strcmp(transp,'transp')
    y = zeros( npix,1 );
else
    y = zeros( reconStruct_small.nfreqs*reconStruct_small.nsamples*reconStruct_small.nechoes*reconStruct_small.nrots*reconStruct_small.ncoils,1 );
end


% GPU flag
if use_gpu == 1
    x = gpuArray(x);
    y = gpuArray(y);
end


for f = 1:reconStruct_small.nfreqs

    for rots = 1:reconStruct_small.nrots  
        EXP = exp( -1j*traj_mat*( SEM_mat(:,rots).' + reconStruct_small.freqs(f) ) );    

        for c = 1:reconStruct_small.ncoils
            B1M = repmat( b1_mat(:,c,rots).',[reconStruct_small.nsamples 1] );

            for ee=1:reconStruct_small.nechoes
                B1P = repmat( b1_plus(:,ee).',[reconStruct_small.nsamples 1] );

                off_1 = (c-1)*reconStruct_small.nsamples + (rots-1)*reconStruct_small.nsamples*reconStruct_small.ncoils + (ee-1)*reconStruct_small.nsamples*reconStruct_small.ncoils*reconStruct_small.nrots;
                off_2 = (f-1)*reconStruct_small.nsamples*reconStruct_small.nechoes*reconStruct_small.nrots*reconStruct_small.ncoils;
                ind_start = 1 + off_1 + off_2;
                inds = ind_start : ind_start + reconStruct_small.nsamples - 1;

                A = EXP .* B1P .* B1M * sqrt(reconStruct_small.freq_weights(f));

                if strcmp(transp,'transp')  % A'*y                
                    y = y + A' * x(inds);
                else  % A*x
                    y(inds) = A * x;
                end

            end

        end

    end
end


if use_gpu == 1
    y = gather(y);
end

% reconStruct_small.g.FreeMemory






















