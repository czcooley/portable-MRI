function [B0, Gx, Gz] = interpmap_1x_20190807( FOVx, FOVy, FOVz, Nx, Ny, Nz )

%interp from load('B0mat_seq16_interp2'); load('Gx_Gz_seq16_interp2');


% FOVx = 0.18;
% FOVy = 0.18;
% FOVz = 0.18;
% Nx = 10;
% Ny = 64;
% Nz = 64;

load('B0mat_seq16_interp2'); load('Gx_Gz_seq16_interp2');
load('seq16_coord.mat');


%% CREATE GRID FOR FIELD MAP RECONSTRUCTION

xlinq = linspace(FOVx,FOVx,1);

ylinq = linspace(-FOVy/2,FOVy/2,Ny);
zlinq = linspace(-FOVz/2,FOVz/2,Nz);

[Z, Y, X] = meshgrid(interp(zlin,2), ylin, xlin);
[Zq, Yq, Xq] = meshgrid(zlinq, ylinq, xlinq);

Gx = interp3(Z, Y, X, Gx_field_rotated_permuted, Zq, Yq, Xq);
Gz = interp3(Z, Y, X, Gz_field_rotated_permuted, Zq, Yq, Xq);
B0 = interp3(Z, Y, X, B0_field_rotated_permuted, Zq, Yq, Xq) - B0_field_rotated_permuted(ceil(end/2),ceil(end/2), ceil(end/2));

% cvec = [-2e-4,2e-4];
% figure; mosaic1(Gz,5,5);  colormap jet; caxis(cvec);
% figure; mosaic1(Gz_field_rotated_permuted,5,5);  colormap jet;  caxis(cvec);
% 
% cvec = [80e-3,82e-3];
% figure; mosaic1(B0,5,5);  colormap jet; caxis(cvec);
% figure; mosaic1(B0_field_rotated_permuted,5,5);  colormap jet;  
% caxis(cvec);

