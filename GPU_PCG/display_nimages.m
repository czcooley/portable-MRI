
function display_nimages(vol)
% function display_nimages(vol)

[nx ny nz]=size(vol);

% fig_handle=figure;

nfigs=ceil(sqrt(nz));
if nfigs*(nfigs-1)>=nz
    nfigsx=nfigs-1;
    nfigsy=nfigs;
else
    nfigsx=nfigs;
    nfigsy=nfigs;
end
    
to_disp=[];
k=1;
for i=1:nfigsx
    row=[];
    for j=1:nfigsy
        if k<=nz
            row=[row vol(:,:,k)];
        else
            row=[row zeros(nx,ny)];
        end
        k=k+1;
    end
    to_disp=[to_disp;row];
end

imagesc(to_disp);