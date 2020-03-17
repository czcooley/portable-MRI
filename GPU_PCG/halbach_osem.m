


function x = halbach_osem( x0,data,traj_mat,SEM_mat,b1_mat,b1_plus,use_gpu,nsubsets,niters,display )


numsamples = size(traj_mat,1);
[nP nSEM] = size(SEM_mat);
[nP num_slices] = size(b1_plus);
[nP nC nSEM] = size(b1_mat);

data = reshape( data,numsamples,nSEM,num_slices,nC );

% create subsets
fprintf('\tCreating subsets ...\n')

rots_subs = cell(nsubsets,1);
all_subs = [];
j = 1;
for i = 1:nsubsets
    rots_subs{i} = [j:nsubsets:nSEM]';
    all_subs = [ all_subs;rots_subs{i} ];
    j = j+1;
end

c = setxor([1:nSEM]',all_subs);
rots_subs{nsubsets} = [ rots_subs{nsubsets};c ];  % append the last rotations to the last subset to make sure we don't forget any data



% compute sensitivity images
fprintf('\tComputing sensitivitiy images');
x_sens = cell(nsubsets,1);
for s = 1:nsubsets
    fprintf('  #%d',s)
    nrots = size(rots_subs{s},1);
    proj = gpuArray(ones(numsamples*nrots*num_slices*nC,1));
    x_sens{s} = forw_mod_subsets_v1( proj,'transp',traj_mat,SEM_mat,b1_mat,b1_plus,use_gpu,rots_subs{s} );
    
    ind = find(abs(x_sens{s})==0);
    x_sens{s}(ind) = 1.0;
    
end
fprintf('\n');


% OSEM
x = x0;
if display == 1
    figure;
    nfigs_x = ceil(sqrt(num_slices));
    nfigs_y = nfigs_x;
    if (nfigs_x-1)*nfigs_y > num_slices
        nfigs_x = nfigs_x-1;
    end
end
for n = 1:niters
    
    fprintf('\tOSEM iteration %d/%d  subset',n,niters);
    
    for s = 1:nsubsets
        
        fprintf(' #%d',s);
        
        % projection
        proj = forw_mod_subsets_v1( x,'notransp',traj_mat,SEM_mat,b1_mat,b1_plus,use_gpu,rots_subs{s} );
        
        % correction ratio
        data2 = gpuArray( data( :,rots_subs{s},:,: ) );
        % corr_ratio = data2(:) ./ proj;
        corr_ratio = data2(:) - proj;
        
        % back-projection
        x_corr = forw_mod_subsets_v1( corr_ratio,'transp',traj_mat,SEM_mat,b1_mat,b1_plus,use_gpu,rots_subs{s} );
        
        % update
        % x = x .* x_corr;
        % x = x + x_corr./x_sens{s};
        x = x + x_corr;
        
        if display==1
            to_disp = conv_3Dstack_2_2Dimg( reshape(x,96,96,num_slices),nfigs_x,nfigs_y );
            imagesc(abs(to_disp)); axis image off; colorbar; caxis([0 1e-8]);
            title( sprintf('Iteration #%d / subset #%d',n,s) );
            drawnow;
        end
        
    end    
    fprintf('\n');
    
end
        
        









 