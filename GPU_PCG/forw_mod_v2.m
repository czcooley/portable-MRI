

function y = forw_mod_v2( x,transp,traj_mat,SEM_mat,b1_mat,B1P,num_slices,use_gpu )


numsamples = size(traj_mat,1);
[nP nSEM] = size(SEM_mat);
[nP nC nSEM] = size(b1_mat);

% flag for computing A*x ('no transp') or A'*y ('transp')
% this required for LSQR
if strcmp(transp,'transp')
    y = zeros( nP,1 );
else
    y = zeros( numsamples*num_slices*nSEM*nC,1 );
end
 
N = numsamples*num_slices;

% GPU flag
if use_gpu == 1
    x2 = gpuArray(x);
end


% start matrix multiplication
j = 1;
for rots = 1:nSEM    
    
    EXP = repmat( exp( -1j*traj_mat*SEM_mat(:,rots).' ),[num_slices 1] );
    
    for coil = 1:nC
        
        B1M = repmat( b1_mat(:,coil,rots).',[N 1] );
        
        fprintf('%d/%d\n',j,nSEM*nC);
        j = j+1;
        
               
        if strcmp(transp,'transp')  % compute A*x
            
            ind_start = 1 + (rots-1)*N + (coil-1)*nSEM*N;
            inds = ind_start : ind_start + N - 1;

            if use_gpu ==1
                y(inds) = ( EXP .* B1P .* B1M ) * x2;
            else
                y(inds) = ( EXP .* B1P .* B1M ) * x;
            end
            
        else  % compute A'*y
            
            ind_start = 1 + (rots-1)*N + (coil-1)*nSEM*N;
            inds = ind_start : ind_start + N - 1;

            if use_gpu ==1
                y(inds) = ( EXP .* B1P .* B1M ) * x2;
            else
                y(inds) = ( EXP .* B1P .* B1M ) * x;
            end
        end
        
    end
end























