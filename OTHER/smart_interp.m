% interpolates the input matrix out to a size specified by N

function out = smart_interp(in, N)

size_in=size(in);

if numel(in) == size_in(1);
    
    nx=numel(in);
    
    out = interp(linspace(1,nx,nx),in,linspace(1,nx,N));
    
else

    nc=numel(in(1,1,:));
    ny=numel(in(:,1,1));
    nx=numel(in(1,:,1));

    for cc=1:nc


%         out(:,:,cc) = interp2(linspace(1,ny,ny), linspace(1,nx,nx), in(:,:,cc),...
%             linspace(1,ny,N),linspace(1,nx,N).');
%         
        
        
        out(:,:,cc) = interp2(linspace(1,nx,nx),linspace(1,ny,ny).', in(:,:,cc), linspace(1,nx,N), linspace(1,ny,N).');


    end


end