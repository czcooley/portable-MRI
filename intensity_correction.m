%% intensity correction across images.  Compensates for B1- variation of the RF coil.

clear all;

load('T2_images_example.mat');

thresh = 0.02;
disksize = 50;
jj=1;

for ii = 2:size(image_all,3)-1
    
    load(['mask',num2str(ii)]);
    bw = smart_interp2d(bw, size(image_all,1),size(image_all,2));

    
    im = abs(image_all(:,:,ii));
    im_mask = im.*bw;
    [I1,I2,Y] =  find(im_mask < thresh);
    
    im_thresh = im_mask;
    for iii = 1: numel(I1)
        im_thresh(I1(iii),I2(iii)) = 1;
    end
    
    h = fspecial('disk', disksize);
    im_LP = imfilter(im_thresh,h,'replicate');
    
    im_corr(:,:,jj) = im_mask./im_LP;
    jj = jj+1;
    
end


figure;
subplot(2,1,1);
mosaic1(abs(image_all(:,:,2:end)),2,5);
caxis([0, 5]);
axis equal;
title('original recon');
subplot(2,1,2);
title(['theshold = ',num2str(thresh),', disksize = ',num2str(disksize)]);
mosaic1(im_corr,2,5);
caxis([0, 2.5])
axis equal
title('intensity corrected');



