% clear all;
D = gpuDevice;

%% helper functions
mip = @(x,dim) squeeze(max(rescale(x),[],dim))';
pwr = @(x,dim) squeeze(log10(max(abs(fftshift(x)),[],dim)))';

%% load data

if ~exist('imPath', 'var')
    
    dataPath = '/archive/bioinformatics/Danuser_lab/Fiolka/Manuscripts/OPM-ALIAS/DataToShare/mesoOPM';
    
    experimentName = 'mesoOPM';
    
    % CASE 1: You have critically sampled ground truth, simulate the
    % undersampling and reconstruct ---------------------------------------
        
        imPath = fullfile(dataPath, 'Cell1', '1_CH00_000000.tif');
        fullySampled = true;
    
    % ---------------------------------------------------------------------
    
    % CASE 2: You've acquired real downsampled data (no ground truth)
    % ---------------------------------------------------------------------
    
        % imPath = fullfile(dataPath, '', 'Cell1_dsp4.tif'); 
        % fullySampled = false;
    
    % ---------------------------------------------------------------------
end

%% parameter setup

% defaults
baseThick = 0;
zCrop = 0;

switch experimentName
    case 'highOPM'
        osFactor = 1;
        dsFactor = 4;
        xyPixelSize = 0.147;
        dz = 0.207;
        skewAngle = 45.0;
        fillMethod = 'median';
        angleBump = 3.0;
        outSize = [16,256,64];
    
    case 'highOPM2'
        osFactor = 2;
        dsFactor = 2;
        xyPixelSize = 0.147;
        dz = 0.207;
        skewAngle = 45.0;
        fillMethod = 'crop';
        angleBump = 3.0;
        outSize = [16,256,96];
        baseThick = 0.5;
        zCrop = 56;
    
    case 'mesoOPM'
        osFactor = 1;
        dsFactor = 4;
        xyPixelSize = 1.15;
        dz = 1.60;
        skewAngle = 45.0;
        fillMethod = 'median';
        angleBump = 1.0;
        outSize = [128,256,64];

    otherwise
        osFactor = 1;
        dsFactor = 1;
        xyPixelSize = 0.147;
        dz = 0.210;
        skewAngle = 45.0;
        fillMethod = 'median';
        angleBump = 0.0;   
        outSize = [256,256,64];
end

dz = dz * osFactor;

%% data preparation

im = permute(readtiff(imPath), [2,1,3]);
im = gpuArray(single(im));

im_dsp = im;

im_full = zeros(size(im_dsp) .* [1 1 dsFactor], 'like', im); % dummy image

% if the input is critically sampled (or better) then we treat it as
% ground truth and simulate the undersampling
if fullySampled
    % oversample -> critical
    im_full = im(:, :, 1:osFactor:end);
    
    % critical -> downsample
    im_dsp = im_full(:, :, 1:dsFactor:end);
end

% comb upsampling operation
im_comb = zeros(size(im_full), 'like', im);
im_comb(:, :, 1:dsFactor:end) = im_dsp;

%%
fillVal = median(im_dsp(:));

% deskewing

im_dsp = deskewFrame3D( ...
    im_dsp, ...
    skewAngle, ...
    dz * dsFactor, ...
    xyPixelSize, ...
    'reverse', true ...
    );

im_comb = deskewFrame3D( ...
    im_comb, ...
    skewAngle, ...
    dz, ...
    xyPixelSize, ...
    'reverse', true ...
    );

im_full = deskewFrame3D( ...
    im_full, ...
    skewAngle, ...
    dz, ...
    xyPixelSize, ...
    'reverse', true ...
    );

% try to reduce edge effects
if strcmp(fillMethod, 'crop')
    x1 = find(squeeze(im_dsp(1,:,1)) > 0);
    x1 = x1(1);
    
    x2 = find(squeeze(im_dsp(1,:,end)) > 0);
    x2 = x2(end);
    
    im_dsp = im_dsp(:, (x1+1):(x2-1), :);
    im_full = im_full(:, (x1+1):(x2-1), :);
elseif strcmp(fillMethod, 'median')
    im_dsp(im_dsp(:) == 0) = fillVal; 
end

%% Fourier space masking and plotting

nz_ds = size(im_dsp,3);

G_dsp  = fftn(im_dsp);
G_comb = fftn(im_comb);
G_full = fftn(im_full);

% handle even downsample case;
z_ds = 0;
if ~mod(dsFactor, 2)
    z_ds = uint16(size(G_comb,3)/dsFactor) + mod(size(G_comb,3),2);
end

% deskew and reciprocal space figures...
figure(1); clf;
set(gcf, 'color', [1,1,1]);
t = tiledlayout(2, 3, 'TileSpacing', 'none', 'Padding', 'compact');

nexttile(t);
imagesc(mip(im_full, 1));
colormap(gca, 'parula');
axis image;

if fullySampled
    title('Critical Nyquist sampling', 'FontSize', 16);
end

%%
nexttile(t);
imagesc(mip(im_dsp, 1));
colormap(gca, 'parula');
axis image;
set(gca, 'YTick', []);
title(sprintf('%dX downsampled', dsFactor), 'FontSize', 16);

nexttile(t);
imagesc(mip(im_comb, 1));
colormap(gca, 'parula');
axis image;
set(gca, 'YAxisLocation', 'right');
title('Comb upsampling', 'FontSize', 16);

%%
nexttile(t);
imagesc(pwr(G_full, 1));
colormap(gca, 'hot');
axis image;

nexttile(t);
imagesc(pwr(G_dsp, 1));
colormap(gca, 'hot');
axis image;
set(gca, 'YTick', []);

nexttile(t);
imagesc(pwr(G_comb, 1));
colormap(gca, 'hot');
axis image;

%%
set(gca, 'YAxisLocation', 'right');
[sy, sx, sz] = size(G_comb);

[x, y, z] = meshgrid(1:sx, 1:sy, 1:sz);
x = x - mean(x(:));
y = y - mean(y(:));
z = z - mean(z(:));

blurSize = 0.0;

mask = (z > -x .* (sz/sx) ./ tand(skewAngle + angleBump) - sz/dsFactor/2) ...
     & (z < -x .* (sz/sx) ./ tand(skewAngle + angleBump) + sz/dsFactor/2);

hold on;

notMask = mip(~mask,1);
tri1 = bwtraceboundary(notMask, [1, 1], "E");
tri2 = bwtraceboundary(notMask, [size(notMask,1), size(notMask,2)], "W");

fill(tri1(:,2), tri1(:,1), [1,1,1], 'FaceAlpha', 0.2, 'EdgeColor', 'w');
fill(tri2(:,2), tri2(:,1), [1,1,1], 'FaceAlpha', 0.2, 'EdgeColor', 'w');

for j = 1:(dsFactor - mod(dsFactor,2))
    yline((j - ~mod(dsFactor,2) * 0.5) * nz_ds, '--', 'color', 'cyan');
end

hold off;

if blurSize > 0
    mask = imgaussfilt3(double(mask), blurSize/dsFactor);
end

%%
G_mask = fftshift(mask) .* G_comb;

% masked G
figure(2); clf; 
set(gcf, 'color', [1,1,1]);
imagesc(pwr(G_mask, 1));
colormap hot;
axis image;
title('G masked');

%%
g_recon = ifftn(G_mask);
g_recon = real(g_recon);

figure(3); clf; 
set(gcf, 'color', [1,1,1]);
imagesc(mip(g_recon, 1));
colormap parula;
axis image;
title('Reconstructed (deskewed)');

%%
rz = size(im_full,3) / size(im_dsp,3);

S = [1 0 0 0
    0 1 0 0
    0 0 rz 0
    0 0 0 1];

im_interp = imwarp(im_dsp, affine3d(S));

figure(4); clf;
imagesc(mip(im_interp, 1));
axis image;
title('Interpolated (deskewed)');

%% rotations

im_full_rot = rotateFrame3D( ...
    im_full, ...
    skewAngle, ...
    osFactor, ...
    'reverse', true, ...
    'Crop', true, ...
    'outSize', outSize ...
    );
im_full_rot = norm_u16(im_full_rot);

if zCrop
    im_full_rot = im_full_rot(:, :, 1:zCrop);
end

% rotated results
figure(5); clf;
set(gcf, 'color', [1,1,1]);
subplot(3,1,1);
imagesc(mip(im_full_rot, 1));
set(gca, 'XTick', []);
ylabel('fully sampled');
axis image; 

im_recon_rot = rotateFrame3D( ...
    g_recon, ...
    skewAngle, ...
    osFactor, ...
    'reverse', true, ...
    'Crop', true, ...
    'outSize', outSize ...
    );
im_recon_rot = norm_u16(im_recon_rot);

if zCrop
    im_recon_rot = im_recon_rot(:, :, 1:zCrop);
end

subplot(3,1,2);
imagesc(mip(im_recon_rot, 1));
set(gca, 'XTick', []);
ylabel('reconstructed');
axis image;

im_interp_rot = rotateFrame3D( ...
    im_interp, ...
    skewAngle, ...
    osFactor, ...
    'reverse', true, ...
    'Crop', true, ...
    'outSize', outSize ...
    );
im_interp_rot = norm_u16(im_interp_rot);

if zCrop
    im_interp_rot = im_interp_rot(:, :, 1:zCrop);
end

subplot(3,1,3);
imagesc(mip(im_interp_rot, 1));
ylabel('interpolated');
axis image;

%% clear GPU cache

reset(D);

%% saving:
% writetiff(im_full_rot, fullfile(dataPath, "im_full_rot.tif"));
% writetiff(im_recon_rot, fullfile(dataPath, "im_recon_rot.tif"));
% writetiff(im_interp_rot, fullfile(dataPath, "im_interp_rot.tif"));

%% functions

function [ out ] = norm_u16( in )
    out = double(in);
    out = out - median(out(:));
    out(out < 0) = 0;
    out = out ./ max(out(:));
    out = uint16(65535 .* out);
end