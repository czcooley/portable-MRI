# portable-MRI
 reconstruction and processing code for portable MRI data

 - [Overview](#overview)
 - [System Requirements](#system-requirements)
 - [Installation Guide](#installation-guide)
 - [License](#license)


 # Overview
 ``portable-MRI`` contains MATLAB code for generalized reconstruction and processing of MRI data generated with a portable MRI scanner. The portable MRI scanner has a built-in (always on) readout gradient. 3D RARE sequences used with phase encoding in the phase-encoding in 2 dimensions. The input data to the reconstruction code is a .mat file containing raw phase-encoded spin-echo train data.

 <b>Image Reconstruction</b> <br />
`HALBACH_reconstruct_simp.m` is the top level script for running the image reconstruction code. An exemplary in vivo dataset (T2 weighted) is provided (T2_data_example.mat).

<b>Image Processing</b><br />
Further image processing is performed with `intensity_correction.m`. This script masks each 2D image with a custom mask and normalizes the image intensity using alow-pass filtered verison of the image. Reconstructed images are provided for loading into this script (T2_images_example.mat). Custom masks can be created using the `drawfreehand` and `createMask` built-in MATLAB functions.

 # System Requirements
 ## Hardware requirements
 `portable-MRI` package requires a standard computer with enough RAM to support the in-memory operations. Accelerated reconstruction requires GPUs.

 ## Software requirements
 This package is supported for *MATLAB V2018b*.

 ### MATLAB Dependencies
 `portable-MRI` mainly depends on the following MATLAB Add-on Toolboxes.

 ```
Image Processing Toolbox
Signal Processing Toolbox
Parallel Computing Toolbox
 ```

 # Guide:

 - git clone https://github.com/czcooley/portable-MRI
 - add all folders to path in MATLAB
 - run `Halbach_reconstruct_simp.m` from top level folder to generate images from raw data
 - run `intensity_correction.m` to normalize the image intensity of the partitions




 # License

 This project is covered under the **Apache 2.0 License**.
