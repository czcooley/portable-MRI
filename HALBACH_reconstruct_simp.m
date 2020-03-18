%% Top level reconstruction script for portable MRI with gradient phase encoding in X and Z and built-in permanent readout gradient in Y
%% FFT is performed to partition in data in X (along the magnet cylinder bore)
%% Then generalized reconstruction is performed on 2D partitions with measured B0 and Gz maps

clear all;

%% LOAD DATA
%%T2 weighted spectral echo data
load('T2_data_example.mat');

%% FIELDMAP FILES
B0name = 'B0_allcomps_noHF_seq16_interp2_191011'; %% matrix vector components
Gzname = 'Gz_allcomps_noHF_seq16_interp2_191011'; %% matrix vector components
Gxname = 'Gx_allcomps_noHF_seq16_interp2_191011'; %% matrix vector components

%% USER  INPUTS FOR RECONSTRUCTION
% reconstruction settings - some vary per acquisition
BW= 100e3;                  %%acquisition bandwidth
FOVy = 0.22;                %%recon FOV (m) in YZ (axial) plane
FOVz = 0.18;                %%recon FOV (m) in YZ (axial) plane
N_recony = 2*FOVy*1000;     %%recon matrix size in YZ (axial) plane
N_reconz = 2*FOVz*1000;     %%recon matrix size in YZ (axial) plane
N_reconx = 1;

niters = 5;        %% PCG iteration number
Igzmax = 9.2;    %% max current in Gz coil
gzoffset = 7.5e-5; %% offset in Gz field
Igxmax = 3;        %% max current in Gx coil (not used for 2d recon)
parthick = 0.007;  %% 7 mm partition thickness (assumed for 2D recon)
Xoffset = -0.001;  %% offset between fft center slice and measured center slice
B0scale = 1.08;     %% scale factor applied to B0 map
B0_vec = [-1500];  %% freq offset of field map (may differ for each data set - temperature dependant)
usemask = 1;       %% mask recon FOV


%% SET UP VARIABLES FOR ENCODING MATRIX
%define data dimension
N_ro = size(data_all,1);  %% # of readout pts used for recon
N_pe =  size(data_all,2);   %% # of phase encode pts in Z (in plane) used for recon
N_pex = 1;  %% set to 1 for 2D recons (using FFT along x)
% calc readout time
readout_time = N_ro/(BW);                                   %%total readout time in seconds
time = linspace(-readout_time/2,readout_time/2,N_ro);       %%read out time vector of data in seconds
% calc x paritition positions    
xslicevec = [parthick*11:-parthick:-parthick*11] - Xoffset;  
xslicevec(5:8) = [0.056,0.048,0.04, 0.032];  %% the xgradient efficiency is lower near the neck, so partitions are thicker


%% loop recon for 2D partitions
count = 1;

for slicecount = 5:16
    
   %% use fft partitioned data for each 2d recon
   FOVx = xslicevec(slicecount);
   data_fftX = fftshift(fft(fftshift(data_all),[],3));
   data_partition = data_fftX(:,:,slicecount);
   datarecon = data_partition;
    
    %% SET UP FIELDMAPS FOR RECON
    disp('generating field maps');
    %%load interpolate maps to recon size
    [B0, Gx, Gz] = interpmap_1x_20191011( FOVx, FOVy, FOVz, N_reconx, N_recony, N_reconz, B0name, Gxname, Gzname ); 
    
    %%set up mask
    load('mask_small.mat'); bw = smart_interp2d(bw, N_recony,N_reconz);
    if usemask == 0; bw = ones(size(bw)); end
    [I] = find(bw);   %indices of non-zeros in mask
      
    %%setup B0 map (Gy permanet readout fieldmap)
    %scale B0 map and substract center field value
    field_maps_recon = (B0)*42.58e6*B0scale;  %% Hz
    field_maps_recon  = field_maps_recon(:,:,:,3) - field_maps_recon(ceil(end/2),ceil(end/2),1,3);
    field_maps_recon_masked = field_maps_recon(I);     %%% permanent magnet readout field map (Y)  
    
    %%setup Gz maps (forshot-to-shot PE in Z)
    %substract Gz offset value
    Gz_unit_field_3d = Gz(:,:,:,3)-gzoffset;  %use z component only  
    %calculate multiplier for G fields for gradient moments in sequence
    Gzscale_recon = -Igzmax*(1/2)*linspace(-1,1,N_pe)*42.576e6*0.666e-3;
    %scale the "unit" Gz map by multiplier for each PE
    G_maps_recon_masked = Gz_unit_field_3d(I)*Gzscale_recon;  %% Gz field map for each PEz

    
    %% CREATE RECONSTRUCTION STRUCTURE
    reconStruct = [];
    reconStruct.recondata(:,:,:,1) = datarecon;
    reconStruct.SEM_mat(:,:,:,:) = repmat(field_maps_recon_masked, 1, 1, 1, N_pe)  + B0_vec;
    reconStruct.grad_mat = 2*pi*G_maps_recon_masked ;
    reconStruct.numsamples = N_ro;
    reconStruct.traj_mat = squeeze(single([2*pi*time]')); %Switch sign on 2pi of mrsolutions
    reconStruct.nT = N_ro;
    reconStruct.nC = 1;
    reconStruct.TE = 2*readout_time;
    reconStruct.time = time;
    reconStruct.reconSize = size(field_maps_recon_masked);
    reconStruct.ngrad = 1;
    reconStruct.b0map=0;
    reconStruct.maxIters=niters;
    [XX,YY] = ndgrid(linspace(-FOVy/2,FOVy/2,N_recony),linspace(-FOVz/2,FOVz/2,N_reconz)); ZZ = zeros(size(XX));
    reconStruct.XX = XX;  reconStruct.YY = YY;  reconStruct.ZZ = ZZ;
    reconStruct.b1_mat = ones(numel(I),1,1, 1, N_pe);
    reconStruct.b1_plus = ones(numel(I),1,1, 1, N_pe);
    reconStruct.b1_plus_echo_exp = [-1     3    -5     7    -9    11];
    reconStruct.b1_plus_echo_exp = 1*reconStruct.b1_plus_echo_exp;
    reconStruct.reconSize(3) = numel(reconStruct.recondata(1,1,1,:));
    reconStruct.numPE = 1;
    reconStruct.nrots=1;
    reconStruct.lambda=0.1;
    reconStruct.I = I;
    reconStruct.useGPU = 0;
    
    %% RUN PCG
    tol = 1e-6; display = 1; freqs = 0; freq_weights = 1;
    figure(100); pause(0.1);
    [img] = RECON(reconStruct,niters,tol,display,freqs,freq_weights);
    recon_all(:,:,count) = img;
    
    count = count+1;
end


%% reshape images
image_all = zeros([size(field_maps_recon),count-1]);
for ii = 1:count-1
    temp = zeros(size(field_maps_recon));
    temp(I) = recon_all(:,1,ii);
    image_all(:,:,ii) = (temp);
end


%% display images
figure; mosaic1(abs(image_all),3,5); colormap gray;
caxis([0 5]); axis equal

save('test_T2_images.mat', 'image_all');
