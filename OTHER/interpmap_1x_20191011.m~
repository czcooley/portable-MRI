function [B0, Gx, Gz] = interpmap_1x_20191011( FOVx, FOVy, FOVz, Nx, Ny, Nz, B0name, Gxname, Gzname )


load(B0name); load(Gxname); load(Gzname);
load('seq16_coord.mat');


%% CREATE GRID FOR FIELD MAP RECONSTRUCTION

xlinq = linspace(FOVx,FOVx,1);
ylinq = linspace(-FOVy/2,FOVy/2,Ny);
zlinq = linspace(-FOVz/2,FOVz/2,Nz);

[Z, Y, X] = meshgrid(interp(zlin,2), ylin, xlin);
[Zq, Yq, Xq] = meshgrid(zlinq, ylinq, xlinq);

Gxx = interp3(Z, Y, X, Gx_field_rotated_permuted(:,:,:,1), Zq, Yq, Xq);
Gxy = interp3(Z, Y, X, Gx_field_rotated_permuted(:,:,:,2), Zq, Yq, Xq);
Gxz = interp3(Z, Y, X, Gx_field_rotated_permuted(:,:,:,3), Zq, Yq, Xq);
Gx(:,:,:,1) = Gxx;
Gx(:,:,:,2) = Gxy;
Gx(:,:,:,3) = Gxz;


Gzx = interp3(Z, Y, X, Gz_field_rotated_permuted(:,:,:,1), Zq, Yq, Xq);
Gzy = interp3(Z, Y, X, Gz_field_rotated_permuted(:,:,:,2), Zq, Yq, Xq);
Gzz = interp3(Z, Y, X, Gz_field_rotated_permuted(:,:,:,3), Zq, Yq, Xq);
Gz(:,:,:,1) = Gzx;
Gz(:,:,:,2) = Gzy;
Gz(:,:,:,3) = Gzz;

B0x = interp3(Z, Y, X, B0_field_rotated_permuted(:,:,:,1), Zq, Yq, Xq);
B0y = interp3(Z, Y, X, B0_field_rotated_permuted(:,:,:,2), Zq, Yq, Xq);
B0z = interp3(Z, Y, X, B0_field_rotated_permuted(:,:,:,3), Zq, Yq, Xq);
B0(:,:,:,1) = B0x;
B0(:,:,:,2) = B0y;
B0(:,:,:,3) = B0z;
