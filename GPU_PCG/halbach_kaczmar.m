

function x = halbach_kaczmar( x0,data,traj_mat,SEM_mat,b1_mat,b1_plus,use_gpu,niters,display )


numsamples = size(traj_mat,1);
[nP nSEM] = size(SEM_mat);
[nP num_slices] = size(b1_plus);
[nP nC nSEM] = size(b1_mat);

data = reshape( data,numsamples,nSEM,num_slices,nC );

if display == 1
    figure;
    nfigs_x = ceil(sqrt(num_slices));
    nfigs_y = nfigs_x;
    if (nfigs_x-1)*nfigs_y >= num_slices
        nfigs_x = nfigs_x-1;
    end
end

x = zeros(size(x0));  % Kaczmar must start @ 0

if use_gpu == 1
    x = gpuArray(x);
end

sbuf = [];

for n = 1:niters
    
    fprintf('\tKaczmar iteration %d/%d ',n,niters);

    for rots = 1:nSEM    
                
        if mod(rots,10) == 0
            for j=1:size(sbuf,2)
                fprintf('\b');
            end
            sbuf = sprintf('[%d/%d]',rots,nSEM);
            fprintf('%s',sbuf);
        end        
        
        EXP = exp( -1j*traj_mat*SEM_mat(:,rots).' );

        for coil = 1:nC
            B1M = repmat( b1_mat(:,coil,rots).',[numsamples 1] );

            for ee=1:num_slices
                B1P = repmat( b1_plus(:,ee).',[numsamples 1] );            
                
                A = EXP .* B1P .* B1M;
                
                for pts = 1:numsamples                

                    % projection
                    proj = A(pts,:)*x;

                    % backprojection
                    if use_gpu ==1
                        data2 = gpuArray(data(pts,rots,num_slices,coil));
                    else
                        data2 = data(pts,rots,num_slices,coil);
                    end

                    % norm_fac = sum(A,1).';               
                    % norm_fac = sum(A,2);
                    % norm_fac = sum( abs(A(:)).^2 );
                    norm_fac = real( A(pts,:) * A(pts,:)' );

                    % x = x + 0.1*A'*( ( data2-proj )./norm_fac );
                    x = x + 0.05 * A(pts,:)'*(data2-proj ) / norm_fac;
                    % x = 0.1 * A'*(data2-proj );
                    
                end
                               
            end
            
        end
        
        if display==1
            to_disp = conv_3Dstack_2_2Dimg( reshape(x,96,96,num_slices),nfigs_x,nfigs_y );
            imagesc(abs(to_disp)); axis image off; colorbar; caxis([0 5]);
            title( sprintf('Iteration #%d / rotation #%d',n,rots) );
            drawnow;
        end        
        
        
    end
    
    fprintf('\n');
    sbuf = [];
    
end

    
    





