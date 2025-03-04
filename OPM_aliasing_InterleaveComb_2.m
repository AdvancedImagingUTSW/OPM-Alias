clear;
% close all;

%%
mip = @(x,dim) squeeze(max(rescale(x),[],dim))';
pwr = @(x,dim) squeeze(log10(max(abs(fftshift(x)),[],dim)))';
p_rl = @(x,dim) squeeze(log10(max(real(fftshift(x)),[],dim)))';
p_im = @(x,dim) squeeze(log10(max(imag(fftshift(x)),[],dim)))';
apply = @(mask,x) fftshift(mask) .* x;
% dataPath = '/archive/bioinformatics/Danuser_lab/Fiolka/MicroscopeDevelopment/OPM/ZoomSystem/SiliconeOil/100nmBeads/green/241105';

% % % meso
dataPath='/archive/bioinformatics/Danuser_lab/Fiolka/MicroscopeDevelopment/Aliasing_OPM/meso/Beads/Stretched';
imPath = fullfile(dataPath, 'Cell2', '1_CH00_000000.tif');
% 

%% deskewing

%meso
dsFactor = 3;
xyPixelSize = 1.15;
dz = 1.6;
skewAngle = 45.0;

interpMethod = 'linear';

%% load, deskew and downsample image

im_full = permute(readtiff(imPath), [2,1,3]);

im_dsp = im_full(:, :, 1:dsFactor:end); %Undersampled data

im_dsk = zeros(size(im_full));

im_dsk(:,:,1:dsFactor:end) = im_full(:,:,1:dsFactor:end); %Undersampled, but with zero infill

%%
% deskew
im_dsk = deskewRotateFrame3D( ...
    im_dsk, ...
    skewAngle, ...
    dz , ...
    xyPixelSize, ...
    'reverse', true, ...
    'interpMethod', interpMethod...
    );

im_full = deskewRotateFrame3D( ...
    im_full, ...
    skewAngle, ...
    dz, ...
    xyPixelSize, ...
    'reverse', true, ...
    'interpMethod', interpMethod...
    );

im_interp = deskewRotateFrame3D( ...
    im_dsp, ...
    skewAngle, ...
    dz*dsFactor, ...
    xyPixelSize, ...
    'reverse', true, ...
    'interpMethod', interpMethod...
    );

figure(1); clf;
set(gcf, 'color', [1,1,1]);
imagesc(mip(im_full, 1));
axis image;
% imagesc(rescale(squeeze(im_full(256, :, :))));

figure(2); clf;
set(gcf, 'color', [1,1,1]);
imagesc(mip(im_dsk, 1));
axis image;
% title('stack (deskewed, downsampled)');

%% 
G = fftn(im_dsk);

figure(3); clf;
set(gcf, 'color', [1,1,1]);
imagesc(pwr(G, 1));
colormap hot;
axis image;

pbaspect([1,1,1]);

%% crop-and-stretch method
mid_z = round(size(G,3)/2);
dz = floor(mid_z/dsFactor) - 2*dsFactor;

G_shift = fftshift(G);
G_crop = G_shift(:,:,mid_z-dz-1:mid_z+dz);   %hard coded truncation in Fourier space
G_crop = ifftshift(G_crop);

figure(5); 
imagesc(pwr(G_crop, 1));
colormap hot;
axis image;

im_recon_cns = abs(ifftn(G_crop));
% im_recon_cns = norm_u16(im_recon_cns);
%writetiff(im_recon_dsp, fullfile(dataPath, "im_recon_dsp.tif"));

figure(6); clf;
imagesc(mip(im_recon_cns, 1));
axis image;

z1 = size(im_full, 3);
z2 = size(im_recon_cns, 3);

S = [1 0 0 0
    0 1 0 0
    0 0 z1/z2 0
    0 0 0 1];

im_crop_stretch = imwarp(im_recon_cns, affine3d(S), interpMethod); %interpolate in z. So far not that great.

figure(7); clf;
imagesc(mip(im_crop_stretch, 1));
axis image;

%% tophat method
G_band = G_shift;
G_band(:, :, 1:(mid_z-dz)) = 0;
G_band(:, :, (mid_z+dz):end) = 0;
G_band = ifftshift(G_band);

figure(11); clf;
subplot(2,1,1);
imagesc(pwr(G_band, 1));
colormap(gca, 'hot');
axis image;

im_tophat = abs(ifftn(G_band));
% im_tophat = norm_u16(im_tophat);

subplot(2,1,2);
imagesc(mip(im_tophat, 1));
colormap(gca, 'parula');
axis image;

%%
figure(8); clf;
imagesc(mip(im_interp, 1));
axis image;

%writetiff(im_interp_rot, fullfile(dataPath, "im_interp_rot.tif"));

%%

% errAlias = abs_err(im_full, im_recon_interp);
% errIntrp = abs_err(im_full, im_interp);
% 
% figure(9); clf;
% subplot(2,1,1);
% imagesc(mip(errAlias, 1));
% colormap jet;
% axis image;
% set(gca, "clim", [0, 2]);
% 
% subplot(2,1,2);
% imagesc(mip(errIntrp, 1));
% colormap jet;
% axis image;
% set(gca, "clim", [0, 2]);

%%
figure(10); clf;
set(gcf, 'color', [1,1,1]);

t = tiledlayout(2, 4, 'TileSpacing', 'compact', 'Padding', 'compact');

pct = 0.25;
crop_full = crop3d(im_full, pct);
crop_cns = crop3d(im_crop_stretch, pct);
crop_interp = crop3d(im_interp, pct);
crop_tophat = crop3d(im_tophat, pct);

nexttile(t);
imagesc(mip(crop_full, 1));
colormap(gca, 'parula');
% axis image;
title('ground truth');

nexttile(t);
colormap(gca, 'parula');
imagesc(mip(crop_cns, 1));
% axis image;
title('crop and stretch');

nexttile(t);
colormap(gca, 'parula');
imagesc(mip(crop_tophat, 1));
% axis image;
title('tophat');

nexttile(t);
colormap(gca, 'parula');
imagesc(mip(crop_interp, 1));
% axis image;
title('interpolation alone');

nexttile(t);
imagesc(pwr(fftn(crop_full), 1));
colormap(gca, 'hot');
% axis image;
nexttile(t);
imagesc(pwr(fftn(crop_cns), 1));
colormap(gca, 'hot');
% axis image;
nexttile(t);
imagesc(pwr(fftn(crop_tophat), 1));
colormap(gca, 'hot');
% axis image;
nexttile(t);
imagesc(pwr(fftn(crop_interp), 1));
colormap(gca, 'hot');
% axis image;

%% functions
function [ err ] = sqr_err( a, b )
    a = single(a);
    b = single(b);
    err = (a - b).^2;
end

function [ err ] = abs_err( a, b )
    a = single(a);
    b = single(b);
    err = abs(a - b);
end

function [ crop ] = crop3d( im, pct )
    if pct == 1.0
        crop = im;
        return;
    end

    c = floor(size(im)/2);
    dr = floor(c * pct);

    crop = im((c(1)-dr(1):c(1)+dr(1)), (c(2)-dr(2):c(2)+dr(2)), (c(3)-dr(3):c(3)+dr(3)));
end

function [ out ] = norm_u16( in )
    out = double(in);
    out = out - median(out(:));
    out(out < 0) = 0;
    out = out ./ max(out(:));
    out = uint16(65535 .* out);
end